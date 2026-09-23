program test_drive
   !! Standalone driver (vilma_drive): generate a synthetic reference state and an
   !! ice-thickness forcing on the model Gauss grid, write a sparse &fe3d config
   !! overlaid on the shipped defaults, run vilma_run, and check the output:
   !!   (1) one output slice per forcing slice,
   !!   (2) the bed subsides under the growing ice cap (z_bed < z_bed_eq),
   !!   (3) the ocean draws down (rsl < 0 in the far field).
   use vilma_precision,       only: wp
   use vilma_constants,       only: rad2deg, pi
   use vilma_sht,             only: sht_grid, sht_grid_init, sht_grid_destroy
   use vilma_drive,           only: vilma_run
   use vilma_radial_fe,       only: radial_fe_finalize
   use ncio,               only: nc_create, nc_write_dim, nc_write, nc_read, nc_size
   implicit none

   integer, parameter :: LMAX = 8, NLAT = 16, NPHI = 32, NT = 3
   character(len=*), parameter :: REF   = "obj/test_drive_ref.nc"
   character(len=*), parameter :: FORCE = "obj/test_drive_force.nc"
   character(len=*), parameter :: OUT   = "obj/test_drive_out.nc"
   character(len=*), parameter :: HOR   = "obj/test_drive_hor.nc"
   character(len=*), parameter :: CFG   = "obj/test_drive.nml"
   character(len=*), parameter :: DEFS  = "input/fastearth3d_defaults.nml"

   type(sht_grid), target :: sht
   real(wp), allocatable  :: lon_deg(:), lat_deg(:), z_bed_eq(:,:), h_ice_eq(:,:)
   real(wp), allocatable  :: h_ice(:,:,:), tyr(:), zb(:,:), rsl(:,:)
   real(wp), allocatable  :: ue(:,:), un(:,:), uet(:,:), unt(:,:)
   real(wp) :: thd, sub
   integer  :: i, j, k, u, jice, jocean, nout
   logical  :: ok

   ok = .true.
   call sht_grid_init(sht, LMAX, nlat=NLAT, nphi=NPHI)

   allocate(lon_deg(sht%nphi), lat_deg(sht%nlat))
   lon_deg = sht%lon*rad2deg
   lat_deg = 90.0_wp - sht%colat*rad2deg

   ! --- synthetic reference state (polar land cap over deep ocean) -------------
   allocate(z_bed_eq(sht%nphi,sht%nlat), h_ice_eq(sht%nphi,sht%nlat))
   do j = 1, sht%nlat
      thd = sht%colat(j)*rad2deg
      do i = 1, sht%nphi
         z_bed_eq(i,j) = merge(500.0_wp, -4000.0_wp, thd < 50.0_wp)
      end do
   end do
   h_ice_eq = 0.0_wp

   call nc_create(REF, overwrite=.true.)
   call nc_write_dim(REF, "lon", x=lon_deg, units="degrees_east")
   call nc_write_dim(REF, "lat", x=lat_deg, units="degrees_north")
   call nc_write(REF, "z_bed_eq",  z_bed_eq,  dim1="lon", dim2="lat")
   call nc_write(REF, "h_ice_eq", h_ice_eq, dim1="lon", dim2="lat")

   ! --- forcing: a grounded cap (colat<30) growing 0 -> 1 -> 2 km --------------
   allocate(h_ice(sht%nphi,sht%nlat,NT), tyr(NT))
   tyr = [0.0_wp, 1000.0_wp, 2000.0_wp]            ! years
   do k = 1, NT
      do j = 1, sht%nlat
         thd = sht%colat(j)*rad2deg
         do i = 1, sht%nphi
            h_ice(i,j,k) = merge(1000.0_wp*real(k-1,wp), 0.0_wp, thd < 30.0_wp)
         end do
      end do
   end do

   call nc_create(FORCE, overwrite=.true.)
   call nc_write_dim(FORCE, "lon",  x=lon_deg, units="degrees_east")
   call nc_write_dim(FORCE, "lat",  x=lat_deg, units="degrees_north")
   call nc_write_dim(FORCE, "time", x=tyr,     units="years", unlimited=.true.)
   call nc_write(FORCE, "h_ice", h_ice, dim1="lon", dim2="lat", dim3="time")

   ! --- run config: sparse &fe3d (overlaid on the physics defaults) + a COMPLETE
   ! --- &ctl group (the run config owns &ctl in full; it has no defaults file).
   open(newunit=u, file=CFG, status="replace", action="write")
   write(u,'(a)')    "&fe3d"
   write(u,'(a,i0)') "    lmax = ", LMAX
   write(u,'(a,i0)') "    nlat = ", NLAT
   write(u,'(a,i0)') "    nphi = ", NPHI
   write(u,'(a)')    "/"
   write(u,'(a)')    "&ctl"
   write(u,'(a)')    '    file_forcing = "'//FORCE//'"'
   write(u,'(a)')    '    name_ice     = "h_ice"'
   write(u,'(a)')    '    name_time    = "time"'
   write(u,'(a)')    '    file_ref     = "'//REF//'"'
   write(u,'(a)')    '    name_zbed_eq = "z_bed_eq"'
   write(u,'(a)')    '    name_hice_ref = "h_ice_eq"'
   write(u,'(a)')    '    file_out     = "'//OUT//'"'
   write(u,'(a)')    '    file_hor     = "'//HOR//'"'
   write(u,'(a)')    "    time_init    = -1.0e30"
   write(u,'(a)')    "    time_end     =  1.0e30"
   write(u,'(a)')    "    remap_input  = .false."   ! input already on the Gauss grid
   write(u,'(a)')    '    name_lon     = "lon"'
   write(u,'(a)')    '    name_lat     = "lat"'
   write(u,'(a)')    "    i_eq         = 0"         ! start slice (ice-free) is the reference
   write(u,'(a)')    '    z_bed_ref_file   = ""'
   write(u,'(a)')    '    h_ice_ref_file   = ""'
   write(u,'(a)')    '    z_bed_eq_file    = ""'
   write(u,'(a)')    '    h_ice_eq_file    = ""'
   write(u,'(a)')    '    rsl_restart_file = ""'
   write(u,'(a)')    '    name_z_bed_ref = "bedrock_topography"'
   write(u,'(a)')    '    name_h_ice_ref = "ice_thickness"'
   write(u,'(a)')    '    name_rsl       = "rsl"'
   write(u,'(a)')    '    restart_in_file = ""'
   write(u,'(a)')    "/"
   close(u)

   call vilma_run(CFG, defaults_file=DEFS)

   ! --- checks -----------------------------------------------------------------
   jice   = nearest_row(15.0_wp)
   jocean = nearest_row(70.0_wp)

   nout = nc_size(OUT, "time")
   write(*,'(a,i0)') '   output slices: ', nout
   if (nout /= NT) then
      write(*,'(a)') '   FAIL: expected one output slice per forcing slice'; ok = .false.
   end if

   allocate(zb(sht%nphi,sht%nlat), rsl(sht%nphi,sht%nlat))
   call nc_read(OUT, "z_bed", zb,  start=[1,1,nout], count=[sht%nphi,sht%nlat,1])
   call nc_read(OUT, "rsl",   rsl, start=[1,1,nout], count=[sht%nphi,sht%nlat,1])

   sub = z_bed_eq(1,jice) - zb(1,jice)
   write(*,'(a,f10.3,a)') '   bed subsidence under ice =', sub, ' m'
   write(*,'(a,f10.3,a)') '   far-field ocean rsl      =', rsl(1,jocean), ' m'
   if (sub < 10.0_wp) then
      write(*,'(a)') '   FAIL: no appreciable subsidence under the growing ice cap'; ok = .false.
   end if
   if (rsl(1,jocean) >= 0.0_wp) then
      write(*,'(a)') '   FAIL: building land ice should draw the ocean down'; ok = .false.
   end if

   ! --- the horizontal-displacement file (&ctl file_hor) --------------------------
   ! An axisymmetric cap on a radially symmetric Earth: the displacement is purely
   ! meridional (u_east ≡ 0) and there is no toroidal field at all, so the
   ! toroidal components must be exactly zero, not merely small.
   write(*,'(a,i0)') '   horizontal slices: ', nc_size(HOR, "time")
   if (nc_size(HOR, "time") /= NT) then
      write(*,'(a)') '   FAIL: expected one horizontal slice per forcing slice'; ok = .false.
   end if
   allocate(ue(sht%nphi,sht%nlat), un(sht%nphi,sht%nlat), uet(sht%nphi,sht%nlat), unt(sht%nphi,sht%nlat))
   call nc_read(HOR, "u_east",      ue,  start=[1,1,nout], count=[sht%nphi,sht%nlat,1])
   call nc_read(HOR, "u_north",     un,  start=[1,1,nout], count=[sht%nphi,sht%nlat,1])
   call nc_read(HOR, "u_east_tor",  uet, start=[1,1,nout], count=[sht%nphi,sht%nlat,1])
   call nc_read(HOR, "u_north_tor", unt, start=[1,1,nout], count=[sht%nphi,sht%nlat,1])
   write(*,'(a,f10.4,a,es10.2,a)') '   max|u_north| =', maxval(abs(un)), ' m,  max|u_east| =', &
        maxval(abs(ue)), ' m'
   write(*,'(a,2es10.2)') '   max|toroidal| (east, north) =', maxval(abs(uet)), maxval(abs(unt))
   if (maxval(abs(un)) < 0.01_wp) then
      write(*,'(a)') '   FAIL: no horizontal displacement under the growing cap'; ok = .false.
   end if
   if (maxval(abs(ue)) > 1.0e-12_wp*maxval(abs(un)) .or. maxval(abs(uet)) /= 0.0_wp &
       .or. maxval(abs(unt)) /= 0.0_wp) then
      write(*,'(a)') '   FAIL: an axisymmetric load on a radial Earth moved east or toroidally'; ok = .false.
   end if

   call sht_grid_destroy(sht)

   write(*,'(a)') ''
   if (ok) then
      write(*,'(a)') ' PASS: vilma_run reads ref+forcing, marches the model,'
      write(*,'(a)') '       and writes a sensible subsiding/draw-down output'
   else
      write(*,'(a)') ' FAIL: standalone driver did not all pass'
      call radial_fe_finalize()
      error stop 1
   end if
   call radial_fe_finalize()

contains

   integer function nearest_row(colat_deg) result(jbest)
      real(wp), intent(in) :: colat_deg
      integer  :: j
      real(wp) :: d, dbest
      jbest = 1;  dbest = huge(1.0_wp)
      do j = 1, sht%nlat
         d = abs(sht%colat(j)*rad2deg - colat_deg)
         if (d < dbest) then;  dbest = d;  jbest = j;  end if
      end do
   end function nearest_row

end program test_drive
