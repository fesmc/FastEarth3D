program test_couple_remap
   !! Host-grid coupling path (the CLIMBER-X contract). Init with a regular lon-lat host
   !! grid that differs from the model Gauss grid, drive a host-grid ice load, and check
   !! that the model: builds the remap, returns rsl/z_bed ON THE HOST GRID, conserves
   !! ocean mass through its internal (Gauss-grid) solve, and relaxes sensibly (bed
   !! subsides under the ice, ocean draws down). This exercises the conservative
   !! host->Gauss (ice in) and bilinear Gauss->host (rsl out) legs of fe_remap end to end.
   use fe_precision, only: wp
   use fe_constants, only: pi, sec_per_year
   use fe_params,    only: fe_param_class
   use fe_radial_fe, only: radial_fe_finalize
   use fe_coupling,  only: solid_earth, solid_earth_init, solid_earth_update, solid_earth_finalize
   use coords,       only: grid_class, grid_init
   implicit none

   integer, parameter :: LMAX = 16
   integer, parameter :: NLON = 96, NLAT = 48          ! host grid (regular lon-lat, != Gauss)
   integer, parameter :: NSTEP = 12
   type(fe_param_class) :: p
   type(solid_earth)    :: se
   type(grid_class)     :: hgrid
   real(wp), allocatable :: lon(:), lat(:), z_bed_eq(:,:), h_ice_eq(:,:), h_ice(:,:)
   real(wp) :: net_subs, ocean_rsl, worst_mass, bed_eq_ice
   integer  :: i, j, jice, joc, step
   logical  :: ok

   ok = .true.
   allocate(lon(NLON), lat(NLAT), z_bed_eq(NLON,NLAT), h_ice_eq(NLON,NLAT), h_ice(NLON,NLAT))
   do i = 1, NLON;  lon(i) = -180.0_wp + (real(i,wp)-0.5_wp)*(360.0_wp/NLON);  end do
   do j = 1, NLAT;  lat(j) =  -90.0_wp + (real(j,wp)-0.5_wp)*(180.0_wp/NLAT);  end do
   call grid_init(hgrid, name="host", mtype="latlon", units="degrees", lon180=.true., x=lon, y=lat)

   ! reference: a northern land cap (lat>40°, +500 m) over deep ocean (−4000 m); no ref ice
   do j = 1, NLAT
      do i = 1, NLON
         z_bed_eq(i,j) = merge(500.0_wp, -4000.0_wp, lat(j) > 40.0_wp)
      end do
   end do
   h_ice_eq = 0.0_wp

   jice = nearest(70.0_wp);  joc = nearest(-50.0_wp)     ! a row under the ice; an ocean row
   bed_eq_ice = z_bed_eq(1,jice)

   p%lmax = LMAX;  p%nlat = 2*LMAX;  p%nphi = 4*LMAX
   se%par = p
   call solid_earth_init(se, z_bed_eq, h_ice_eq, grid=hgrid)

   write(*,'(a,i0,a,i0,a,i0,a)') ' couple+remap: host ', NLON, 'x', NLAT, &
        ' -> Gauss (lmax=', LMAX, ')'

   ! the model must have built a remap and exposed host-grid output fields
   if (.not. se%remap) then
      write(*,'(a)') '   FAIL: no remap built for a non-Gauss host grid'; ok = .false.
   end if
   if (size(se%rsl,1) /= NLON .or. size(se%rsl,2) /= NLAT .or. &
       size(se%z_bed,1) /= NLON .or. size(se%z_bed,2) /= NLAT) then
      write(*,'(a)') '   FAIL: output fields not on the host grid'; ok = .false.
   end if

   ! a 2 km grounded ice cap on the land interior (lat>55°), held fixed
   h_ice = 0.0_wp
   do j = 1, NLAT
      do i = 1, NLON
         if (lat(j) > 55.0_wp) h_ice(i,j) = 2000.0_wp
      end do
   end do

   worst_mass = 0.0_wp
   do step = 1, NSTEP
      call solid_earth_update(se, h_ice, 2.0e3_wp)        ! 2 kyr/step [years]
      worst_mass = max(worst_mass, se%worst_mass_resid)
   end do

   net_subs  = bed_eq_ice - se%z_bed(1,jice)
   ocean_rsl = se%rsl(1,joc)
   write(*,'(a,f12.4,a)') '   net bed subsidence under ice (host) =', net_subs, ' m'
   write(*,'(a,f12.4,a)') '   ocean rsl (drawdown, should be < 0)  =', ocean_rsl, ' m'
   write(*,'(a,es11.2)')  '   worst ocean-mass residual            =', worst_mass

   ! checks: mass conserved on the Gauss solve, appreciable subsidence under the ice
   ! (z_bed reconstructed on the host bed), and ocean drawdown from building land ice.
   if (worst_mass > 1.0e-10_wp) then
      write(*,'(a)') '   FAIL: ocean mass not conserved during the Gauss-grid solve'; ok = .false.
   end if
   if (net_subs < 10.0_wp) then
      write(*,'(a)') '   FAIL: no appreciable subsidence under the ice'; ok = .false.
   end if
   if (ocean_rsl >= 0.0_wp) then
      write(*,'(a)') '   FAIL: building land ice should draw the ocean down'; ok = .false.
   end if

   call solid_earth_finalize(se)

   write(*,'(a)') ''
   if (ok) then
      write(*,'(a)') ' PASS: host-grid coupling remaps ice in / rsl out and relaxes sensibly'
   else
      write(*,'(a)') ' FAIL: host-grid coupling path'
      call radial_fe_finalize()
      error stop 1
   end if
   call radial_fe_finalize()

contains

   integer function nearest(lat_deg) result(jbest)
      real(wp), intent(in) :: lat_deg
      integer  :: j
      real(wp) :: d, dbest
      jbest = 1;  dbest = huge(1.0_wp)
      do j = 1, NLAT
         d = abs(lat(j) - lat_deg)
         if (d < dbest) then;  dbest = d;  jbest = j;  end if
      end do
   end function nearest

end program test_couple_remap
