program vilma_remap_main
   !! Offline lon-lat -> Gauss-Legendre remapper. Reads a &vilma config (the grid
   !! knobs lmax/nlat/nphi, the source file_forcing + name_ice/name_lon/name_lat/
   !! name_time, and file_out), conservatively remaps the ice variable over every
   !! time slice onto the model Gauss grid, and writes a Gauss-grid forcing file the
   !! standalone driver can consume with remap_input=.false.
   !!
   !!   ./bin/vilma_remap.x config.nml [defaults.nml]
   !!
   !! The online path (vilma_drive, remap_input=.true.) is the default and needs no
   !! preprocessing; this tool is for pre-baking a remapped forcing when that is
   !! preferred (e.g. reusing one remap across many runs). Same conservative engine
   !! (vilma_remap), so results are identical.
   use vilma_precision, only: wp
   use vilma_params,    only: vilma_param_class, vilma_par_load
   use vilma_control,   only: vilma_ctl_class, vilma_ctl_load, DEFAULTS_FILE
   use vilma_sht,       only: sht_grid, sht_grid_init, sht_grid_destroy
   use vilma_remap,     only: remap_ll_gauss, remap_init, remap_to_gauss
   use ncio,         only: nc_read, nc_size, nc_create, nc_write_dim, nc_write
   implicit none

   real(wp), parameter :: RAD2DEG = 57.295779513082323_wp
   character(len=512)  :: cfg
   type(vilma_param_class) :: p
   type(vilma_ctl_class)   :: c
   type(sht_grid)       :: sht
   type(remap_ll_gauss)   :: rmap
   real(wp), allocatable :: lon_s(:), lat_s(:), src(:,:), gauss(:,:)
   real(wp), allocatable :: lon_g(:), lat_g(:), tyr(:)
   integer :: nlon, nls, np, nl, nt, k, nlat, nphi

   ! --- config (grid from &vilma over the physics defaults; I/O from &ctl) -----
   if (command_argument_count() >= 1) then
      call get_command_argument(1, cfg)
   else
      cfg = "vilma.nml"
   end if
   call vilma_par_load(p, cfg, defaults_file=DEFAULTS_FILE)
   call vilma_ctl_load(c, cfg)
   if (len_trim(c%file_forcing) == 0) error stop 'vilma_remap: file_forcing not set'

   ! --- Gauss grid (same defaulting as the driver) ---------------------------
   nlat = p%nlat;  if (nlat <= 0) nlat = 2*p%lmax + 2
   nphi = p%nphi;  if (nphi <= 0) nphi = 4*p%lmax
   call sht_grid_init(sht, p%lmax, nlat=nlat, nphi=nphi)
   np = sht%nphi;  nl = sht%nlat

   ! --- conservative map from the source axes --------------------------------
   nlon = nc_size(c%file_forcing, trim(c%name_lon))
   nls  = nc_size(c%file_forcing, trim(c%name_lat))
   allocate(lon_s(nlon), lat_s(nls))
   call nc_read(c%file_forcing, trim(c%name_lon), lon_s)
   call nc_read(c%file_forcing, trim(c%name_lat), lat_s)
   call remap_init(rmap, sht, lon_s, lat_s)

   nt = nc_size(c%file_forcing, trim(c%name_time))
   allocate(tyr(nt));  call nc_read(c%file_forcing, trim(c%name_time), tyr)
   write(*,'(a,i0,a,i0,a,i0,a,i0,a,i0,a)') ' vilma_remap: ', nlon, 'x', nls, &
        ' lon-lat -> ', np, 'x', nl, ' Gauss, ', nt, ' slices'

   ! --- Gauss output axes ----------------------------------------------------
   allocate(lon_g(np), lat_g(nl))
   lon_g = sht%lon*RAD2DEG
   do k = 1, nl;  lat_g(k) = 90.0_wp - sht%colat(k)*RAD2DEG;  end do

   call nc_create(c%file_out)
   call nc_write_dim(c%file_out, "lon",  x=lon_g, units="degrees_east")
   call nc_write_dim(c%file_out, "lat",  x=lat_g, units="degrees_north")
   call nc_write_dim(c%file_out, "time", x=tyr,   units="years", unlimited=.true.)

   allocate(src(nlon,nls), gauss(np,nl))
   do k = 1, nt
      call nc_read(c%file_forcing, trim(c%name_ice), src, &
                   start=[1,1,k], count=[nlon, nls, 1])
      call remap_to_gauss(rmap, sht, src, gauss, conserve_mass=.true.)
      call nc_write(c%file_out, trim(c%name_ice), gauss, dim1="lon", dim2="lat", &
                    dim3="time", start=[1,1,k], count=[np, nl, 1])
   end do

   call sht_grid_destroy(sht)
   write(*,'(a,a)') ' vilma_remap: wrote ', trim(c%file_out)
end program vilma_remap_main
