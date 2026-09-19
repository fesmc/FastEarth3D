program test_sht
   !! Round-trip test for the SHTns transform wrapper (fe_sht): synthesize a
   !! known band-limited spectrum to the spatial grid, analyze it back, and check
   !! the coefficients are recovered to machine precision. This proves the build
   !! system, the SHTns linkage, and the spectral<->spatial kernel before any
   !! physics is added.
   use, intrinsic :: iso_c_binding, only: c_double
   use fe_precision, only: wp
   use fe_sht,       only: sht_grid, sht_grid_init, sht_grid_lmidx, sht_grid_synthesis, sht_grid_analysis, sht_grid_surface_integral, sht_grid_destroy
   implicit none

   type(sht_grid) :: g
   integer  :: lmax, l, m, lm, j
   real(wp),    allocatable :: sh(:,:)
   complex(wp), allocatable :: slm(:), slm0(:)
   real(wp) :: err, tol, fourpi, area, integ, y10err
   logical  :: ok

   lmax = 24
   call sht_grid_init(g, lmax, nlat=48, nphi=96)
   print '(a,i0,a,i0,a,i0,a,i0)', ' SHTns grid: lmax=', g%lmax, &
        '  nlat=', g%nlat, '  nphi=', g%nphi, '  nlm=', g%nlm

   allocate(sh(g%nphi, g%nlat), slm(g%nlm), slm0(g%nlm))

   ! Deterministic band-limited spectrum. m=0 coefficients must be real for a
   ! real-valued spatial field; m>0 may be complex.
   slm0 = (0.0_wp, 0.0_wp)
   do l = 0, lmax
      do m = 0, l
         lm = sht_grid_lmidx(g, l, m)
         if (m == 0) then
            slm0(lm) = cmplx(1.0_wp/real(l+1, wp), 0.0_wp, wp)
         else
            slm0(lm) = cmplx(1.0_wp/real(l+1, wp), &
                             0.1_wp*real(m, wp)/real(l+1, wp), wp)
         end if
      end do
   end do

   slm = slm0
   call sht_grid_synthesis(g, slm, sh)    ! spectral -> spatial
   call sht_grid_analysis(g, sh, slm)     ! spatial  -> spectral (round trip)

   err = maxval(abs(slm - slm0))
   tol = 1.0e-11_wp
   fourpi = 4.0_wp*acos(-1.0_wp)
   ok = .true.
   print '(a,es12.4,a,es12.4)', ' round-trip max error = ', err, '   tol = ', tol
   if (err >= tol) ok = .false.

   ! --- Quadrature: surface integral of a uniform field must give 4*pi --------
   sh   = 1.0_wp
   area = sht_grid_surface_integral(g, sh)
   print '(a,es22.14,a,es12.4)', ' area  ∫1 dΩ          = ', area, &
        '   err = ', abs(area - fourpi)
   if (abs(area - fourpi) >= 1.0e-11_wp) ok = .false.

   ! --- Normalization: a unit Y00 coefficient is orthonormal, ∫ Y00^2 dΩ = 1 --
   slm = (0.0_wp, 0.0_wp)
   slm(sht_grid_lmidx(g, 0,0)) = (1.0_wp, 0.0_wp)
   call sht_grid_synthesis(g, slm, sh)
   integ = sht_grid_surface_integral(g, sh*sh)
   print '(a,es22.14,a,es12.4)', ' norm  ∫ Y00^2 dΩ     = ', integ, &
        '   err = ', abs(integ - 1.0_wp)
   if (abs(integ - 1.0_wp) >= 1.0e-11_wp) ok = .false.

   ! --- Geometry + convention: Y(1,0) field = sqrt(3/4pi) * cos(colat) --------
   ! Confirms orthonormal norm, no Condon-Shortley phase, and that g%colat
   ! matches the spatial grid row ordering.
   slm = (0.0_wp, 0.0_wp)
   slm(sht_grid_lmidx(g, 1,0)) = (1.0_wp, 0.0_wp)
   call sht_grid_synthesis(g, slm, sh)
   y10err = 0.0_wp
   do j = 1, g%nlat
      y10err = max(y10err, abs(sh(1,j) - sqrt(3.0_wp/fourpi)*cos(g%colat(j))))
   end do
   print '(a,es12.4)', ' Y(1,0) vs sqrt(3/4pi)cosθ max err = ', y10err
   if (y10err >= 1.0e-12_wp) ok = .false.

   call sht_grid_destroy(g)

   if (ok) then
      print '(a)', ' PASS: SHTns transform, quadrature, normalization, geometry'
   else
      print '(a)', ' FAIL: one or more SHT checks exceeded tolerance'
      error stop 1
   end if
end program test_sht
