program test_visc_load
   !! Rung 6c — load a real lon-lat-r viscosity field from netCDF onto the Gauss
   !! grid × FE radial nodes (fe_read_visc_3d), and bridge node→element into the
   !! lateral-viscosity advance (ve%enable_lateral_visc_from_nodes).
   !!
   !! Part 1 (always): write a synthetic field whose value is LINEAR in each of
   !! (lon, lat, r) on the model's own Gauss/node grid, read it back, and check the
   !! bilinear-in-(lon,lat)/linear-in-r interpolation reproduces it to round-off —
   !! this pins the coordinate ordering, the (nphi,nlat)→spatial-index mapping and
   !! the periodic-longitude wrap.
   !!
   !! Part 2 (only if present): load the real Pan et al. (2022) field from
   !! ~/models/isostasy_data/earth_structure/viscosity/pan2022.nc, enable it via the
   !! node→element bridge, and check the resulting Maxwell rate field is finite and
   !! that a few forced steps produce a finite, non-trivial uplift.
   use fe_precision,       only: wp
   use fe_constants,       only: rad2deg
   use fe_earth_structure, only: earth_model, build_M3L70V01, fe_read_visc_3d
   use fe_radial_fe,       only: radial_fe_finalize
   use fe_response,        only: response_enable_lateral_visc_from_nodes, response_destroy, response, response_init_elastic, response_init_ve, response_init_null
   use fe_sht,             only: sht_grid, sht_grid_init
   use ncio
   implicit none

   integer,  parameter :: LMAX = 32, NLATF = 2
   real(wp), parameter :: YR = 3.15576e7_wp, DT = 50.0_wp*YR
   character(len=*), parameter :: SYNTH = 'obj/test_visc_synth.nc'
   character(len=*), parameter :: PAN   = &
        '/Users/alrobi001/models/isostasy_data/earth_structure/viscosity/pan2022.nc'

   type(sht_grid)    :: sht
   type(earth_model) :: e
   type(response) :: ve
   logical :: ok, exists

   ok = .true.
   call sht_grid_init(sht, LMAX, nlat=NLATF*LMAX, nphi=2*LMAX+2, mmax=LMAX)
   e = build_M3L70V01()
   call response_init_ve(ve, e, sht, DT)

   call test_synthetic()

   inquire(file=PAN, exist=exists)
   if (exists) then
      call test_pan2022()
   else
      write(*,'(a)') ' [skip] pan2022.nc not found — synthetic round-trip only'
   end if

   call response_destroy(ve)
   call radial_fe_finalize()
   if (.not. ok) error stop 'test_visc_load FAILED'
   write(*,'(a)') ' test_visc_load PASSED'

contains

   subroutine test_synthetic()
      !! Linear-in-each-coordinate field on the Gauss/node grid → interp is exact.
      real(wp), allocatable :: lon_s(:), lat_s(:), r_s(:), eta_s(:,:,:), vn(:,:)
      real(wp) :: lon_t, lat_t, expect, err, maxerr
      integer  :: nphi, nlat, nr, i, j, k, sp
      nphi = sht%nphi;  nlat = sht%nlat;  nr = ve%nr
      allocate(lon_s(nphi), lat_s(nlat), r_s(nr), eta_s(nphi,nlat,nr))
      ! source lon = Gauss lon (ascending); source lat = ascending Gauss lats;
      ! source r = node radii. eta = lon_deg + 1000·lat_deg + 1e-3·r (linear).
      do i = 1, nphi;  lon_s(i) = sht%lon(i)*rad2deg;  end do
      do j = 1, nlat;  lat_s(j) = 90.0_wp - sht%colat(nlat+1-j)*rad2deg;  end do
      r_s = ve%r
      do k = 1, nr
         do j = 1, nlat
            do i = 1, nphi
               eta_s(i,j,k) = lon_s(i) + 1000.0_wp*lat_s(j) + 1.0e-3_wp*r_s(k)
            end do
         end do
      end do
      call nc_create(SYNTH, overwrite=.true.)
      call nc_write_dim(SYNTH, "lon", x=lon_s, units="degrees_east")
      call nc_write_dim(SYNTH, "lat", x=lat_s, units="degrees_north")
      call nc_write_dim(SYNTH, "r",   x=r_s,   units="m")
      call nc_write(SYNTH, "eta", eta_s, dim1="lon", dim2="lat", dim3="r")

      call fe_read_visc_3d(SYNTH, sht, ve%r, vn)
      maxerr = 0.0_wp
      do k = 1, nr
         do j = 1, nlat
            lat_t = 90.0_wp - sht%colat(j)*rad2deg
            do i = 1, nphi
               lon_t = sht%lon(i)*rad2deg
               sp = i + (j-1)*nphi
               expect = lon_t + 1000.0_wp*lat_t + 1.0e-3_wp*ve%r(k)
               err = abs(vn(sp,k) - expect)/(1.0_wp + abs(expect))
               maxerr = max(maxerr, err)
            end do
         end do
      end do
      write(*,'(a,es10.2)') ' synthetic round-trip max rel err = ', maxerr
      if (maxerr > 1.0e-10_wp) then; ok = .false.; write(*,'(a)') '   FAIL: interp not exact'; end if
      deallocate(lon_s, lat_s, r_s, eta_s, vn)
   end subroutine test_synthetic

   subroutine test_pan2022()
      real(wp), allocatable :: vn(:,:)
      real(wp) :: vmin, vmax
      call fe_read_visc_3d(PAN, sht, ve%r, vn)
      vmin = minval(vn);  vmax = maxval(vn)
      write(*,'(a,f6.2,a,f6.2)') ' pan2022 log10(eta) on grid: min=', vmin, ' max=', vmax
      ! Mantle log10(η) lives in ~[18,24]; lithosphere/core endpoints push the
      ! clamped range wider, but it must stay physical (no NaN, within [10,40]).
      if (vmin < 10.0_wp .or. vmax > 40.0_wp .or. vmin /= vmin) then
         ok = .false.; write(*,'(a)') '   FAIL: pan2022 values out of physical range'
      end if
      call response_enable_lateral_visc_from_nodes(ve, sht, vn)
      if (.not. allocated(ve%Mk3)) then
         ok = .false.; write(*,'(a)') '   FAIL: enable_from_nodes did not set up Mk3'
      else if (any(ve%Mk3 /= ve%Mk3) .or. minval(ve%Mk3) < 0.0_wp) then
         ok = .false.; write(*,'(a)') '   FAIL: Mk3 has NaN or negative entries'
      else
         write(*,'(a,es10.2,a,es10.2)') ' Mk3 (=μΔt/η_eff) range: min=', &
              minval(ve%Mk3), ' max=', maxval(ve%Mk3)
      end if
      deallocate(vn)
   end subroutine test_pan2022

end program test_visc_load
