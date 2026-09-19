program fastearth_mkref
   !! Offline reference generation: write the 2D reference state (bedrock + ice) onto
   !! a chosen Gauss grid. Reads a &fe3d config (grid knobs lmax/nlat/nphi, the
   !! reference source z_bed_ref_file / h_ice_ref_file + name_z_bed_ref /
   !! name_h_ice_ref / name_lon / name_lat, and file_out), conservatively remaps bed
   !! (as-is) and ice (mass-conserving) once, and writes a Gauss-grid reference file.
   !!
   !!   ./bin/fastearth_mkref.x mkref_l128.nml [defaults.nml]
   !!
   !! Used to generate the canonical reference data/reference/rtopo_gauss_l128.nc from
   !! the 0.5-deg source. Runs at other resolutions remap that canonical reference
   !! online (cached), so regenerating per-resolution files is no longer needed. Both
   !! reference vars are assumed on a common source grid (RTopo).
   use fe_precision, only: wp
   use fe_params,    only: fe_param_class, fe_par_load
   use fe_control,   only: fe_ctl_class, fe_ctl_load, DEFAULTS_FILE
   use fe_sht,       only: sht_grid, sht_grid_init, sht_grid_destroy
   use fe_remap,     only: remap_ll_gauss, remap_init, remap_to_gauss
   use ncio,         only: nc_read, nc_size, nc_create, nc_write_dim, nc_write
   implicit none

   real(wp), parameter :: RAD2DEG = 57.295779513082323_wp
   character(len=512)  :: cfg
   type(fe_param_class) :: p
   type(fe_ctl_class)   :: c
   type(sht_grid)       :: sht
   type(remap_ll_gauss)   :: rmap
   real(wp), allocatable :: lon_s(:), lat_s(:), buf(:,:), bed(:,:), ice(:,:)
   real(wp), allocatable :: lon_g(:), lat_g(:)
   integer :: nlon, nls, np, nl, k, nlat, nphi

   ! --- config (grid from &fe3d over the physics defaults; I/O from &ctl) -----
   if (command_argument_count() >= 1) then
      call get_command_argument(1, cfg)
   else
      cfg = "fastearth.nml"
   end if
   call fe_par_load(p, cfg, defaults_file=DEFAULTS_FILE)
   call fe_ctl_load(c, cfg)
   if (len_trim(c%z_bed_ref_file) == 0) error stop 'fastearth_mkref: z_bed_ref_file not set'
   if (len_trim(c%h_ice_ref_file) == 0) error stop 'fastearth_mkref: h_ice_ref_file not set'

   ! --- Gauss grid (same defaulting as the driver) ---------------------------
   nlat = p%nlat;  if (nlat <= 0) nlat = 2*p%lmax + 2
   nphi = p%nphi;  if (nphi <= 0) nphi = 4*p%lmax
   call sht_grid_init(sht, p%lmax, nlat=nlat, nphi=nphi)
   np = sht%nphi;  nl = sht%nlat

   ! --- conservative map from the reference source axes ----------------------
   nlon = nc_size(c%z_bed_ref_file, trim(c%name_lon))
   nls  = nc_size(c%z_bed_ref_file, trim(c%name_lat))
   allocate(lon_s(nlon), lat_s(nls), buf(nlon,nls), bed(np,nl), ice(np,nl))
   call nc_read(c%z_bed_ref_file, trim(c%name_lon), lon_s)
   call nc_read(c%z_bed_ref_file, trim(c%name_lat), lat_s)
   write(*,'(a,i0,a,i0,a,i0,a,i0,a)') ' fastearth_mkref: ', nlon, 'x', nls, &
        ' lon-lat -> ', np, 'x', nl, ' Gauss (building conservative map)'
   call remap_init(rmap, sht, lon_s, lat_s)

   call nc_read(c%z_bed_ref_file, trim(c%name_z_bed_ref), buf)
   call remap_to_gauss(rmap, sht, buf, bed, conserve_mass=.false.)
   call nc_read(c%h_ice_ref_file, trim(c%name_h_ice_ref), buf)
   call remap_to_gauss(rmap, sht, buf, ice, conserve_mass=.true.)

   ! --- Gauss output ---------------------------------------------------------
   allocate(lon_g(np), lat_g(nl))
   lon_g = sht%lon*RAD2DEG
   do k = 1, nl;  lat_g(k) = 90.0_wp - sht%colat(k)*RAD2DEG;  end do

   call nc_create(c%file_out)
   call nc_write_dim(c%file_out, "lon", x=lon_g, units="degrees_east")
   call nc_write_dim(c%file_out, "lat", x=lat_g, units="degrees_north")
   call nc_write(c%file_out, trim(c%name_z_bed_ref), bed, dim1="lon", dim2="lat", &
                 units="m", long_name="reference bedrock elevation (Gauss grid)")
   call nc_write(c%file_out, trim(c%name_h_ice_ref), ice, dim1="lon", dim2="lat", &
                 units="m", long_name="reference ice thickness (Gauss grid)")

   call sht_grid_destroy(sht)
   write(*,'(a,a)') ' fastearth_mkref: wrote ', trim(c%file_out)
end program fastearth_mkref
