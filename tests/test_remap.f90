program test_remap
   !! Conservative lon-lat -> Gauss remap (fe_remap): constant preservation, latitude
   !! orientation (the ascending<->SHTns-row flip), zonal-field accuracy, and the
   !! optional global mass-rescale (SHTns surface integral == source area integral).
   use fe_precision, only: wp
   use fe_sht,       only: sht_grid, sht_grid_init, sht_grid_surface_integral, sht_grid_destroy
   use fe_remap,     only: remap_ll_gauss, remap_init, remap_to_gauss, remap_to_ll
   implicit none

   integer,  parameter :: LMAX = 32
   integer,  parameter :: NLON = 72, NLAT_S = 36       ! 5-degree source grid
   real(wp), parameter :: DEG = 3.141592653589793_wp/180.0_wp
   real(wp), parameter :: RAD2DEG = 57.295779513082323_wp

   type(sht_grid)      :: sht
   type(remap_ll_gauss)  :: rmap
   real(wp), allocatable :: lon_s(:), lat_s(:), fsrc(:,:), fdst(:,:)
   real(wp) :: latg, emax, e, src_mean, gauss_mean, isrc, idst
   integer  :: i, j
   logical  :: ok

   ok = .true.
   call sht_grid_init(sht, LMAX, nlat=2*LMAX, nphi=4*LMAX)

   allocate(lon_s(NLON), lat_s(NLAT_S))
   do i = 1, NLON;   lon_s(i) = -180.0_wp + (real(i,wp)-0.5_wp)*(360.0_wp/NLON);  end do
   do j = 1, NLAT_S; lat_s(j) =  -90.0_wp + (real(j,wp)-0.5_wp)*(180.0_wp/NLAT_S); end do
   call remap_init(rmap, sht, lon_s, lat_s)

   allocate(fsrc(NLON,NLAT_S), fdst(sht%nphi,sht%nlat))

   ! (1) constant field preserved exactly
   fsrc = 7.0_wp
   call remap_to_gauss(rmap, sht, fsrc, fdst)
   emax = maxval(abs(fdst - 7.0_wp))
   write(*,'(a,es10.2)') '   constant-field max error      = ', emax
   if (emax > 1.0e-10_wp) then; write(*,*) 'FAIL: constant not preserved'; ok = .false.; end if

   ! (2) zonal field f = sin(lat): accuracy + orientation (north row positive)
   do j = 1, NLAT_S
      do i = 1, NLON;  fsrc(i,j) = sin(lat_s(j)*DEG);  end do
   end do
   call remap_to_gauss(rmap, sht, fsrc, fdst)
   emax = 0.0_wp
   do j = 1, sht%nlat
      latg = 90.0_wp - sht%colat(j)*RAD2DEG
      e = maxval(abs(fdst(:,j) - sin(latg*DEG)))
      emax = max(emax, e)
   end do
   write(*,'(a,es10.2)') '   zonal sin(lat) max error      = ', emax
   if (emax > 2.0e-2_wp) then; write(*,*) 'FAIL: zonal field beyond resolution tol'; ok = .false.; end if
   ! orientation: SHTns row 1 is the northern pole-most cell (lat>0 -> sin>0)
   if (fdst(1,1) <= 0.0_wp .or. fdst(1,sht%nlat) >= 0.0_wp) then
      write(*,*) 'FAIL: latitude orientation flipped'; ok = .false.
   end if
   write(*,'(a,f8.3,a,f8.3)') '   north row=', fdst(1,1), '   south row=', fdst(1,sht%nlat)

   ! (3) mass-rescale: a positive bump, conserve_mass=.true. -> SHTns integral matches
   !     the source area-integral (mean) to machine precision.
   do j = 1, NLAT_S
      do i = 1, NLON
         fsrc(i,j) = 1.0_wp + 0.5_wp*cos(lat_s(j)*DEG)*cos(lon_s(i)*DEG)
      end do
   end do
   call remap_to_gauss(rmap, sht, fsrc, fdst, conserve_mass=.true.)
   ! source area-weighted mean and the model's gauss-weighted mean (integral/4pi)
   isrc = sum(fsrc*rmap%ll%area)/sum(rmap%ll%area)
   idst = sht_grid_surface_integral(sht, fdst)/(16.0_wp*atan(1.0_wp))
   write(*,'(a,f12.8,a,f12.8,a,es10.2)') '   src mean=', isrc, '  gauss mean=', idst, &
        '  rel=', abs(isrc-idst)/abs(isrc)
   if (abs(isrc-idst) > 1.0e-10_wp*abs(isrc)) then
      write(*,*) 'FAIL: mass not conserved with conserve_mass'; ok = .false.
   end if

   ! (4) reverse leg: bilinear Gauss -> host lon-lat reproduces a smooth zonal field
   block
      real(wp), allocatable :: gsrc(:,:), fll(:,:)
      real(wp) :: latg2
      allocate(gsrc(sht%nphi,sht%nlat), fll(NLON,NLAT_S))
      do j = 1, sht%nlat
         latg2 = 90.0_wp - sht%colat(j)*RAD2DEG
         gsrc(:,j) = sin(latg2*DEG)            ! smooth field on the Gauss grid (SHTns rows)
      end do
      call remap_to_ll(rmap, gsrc, fll)
      emax = 0.0_wp
      do j = 1, NLAT_S
         e = maxval(abs(fll(:,j) - sin(lat_s(j)*DEG)))
         emax = max(emax, e)
      end do
      write(*,'(a,es10.2)') '   reverse bilinear sin(lat) err = ', emax
      if (emax > 5.0e-2_wp) then; write(*,*) 'FAIL: reverse bilinear beyond tol'; ok = .false.; end if
   end block

   call sht_grid_destroy(sht)
   if (ok) then
      write(*,'(a)') ' PASS: fe_remap conservative lon-lat -> Gauss + bilinear Gauss -> lon-lat'
   else
      write(*,'(a)') ' FAIL: fe_remap'
      error stop 1
   end if
end program test_remap
