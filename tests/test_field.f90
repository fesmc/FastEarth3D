program test_field
   !! Validate the arbitrary-point evaluator (sht_grid_eval_point) and the analytic
   !! field generators (fe_field), the infrastructure the Martinec-2018 SLE
   !! benchmark profiles are built on.
   use fe_precision, only: wp
   use fe_sht,       only: sht_grid, sht_grid_init, sht_grid_lmidx, sht_grid_eval_point, sht_grid_synthesis, sht_grid_eval_point_horiz, sht_grid_destroy, sht_grid_surface_integral
   use fe_field,     only: spherical_cap, exp_basin, angular_distance
   implicit none
   real(wp), parameter :: pi = acos(-1.0_wp)
   integer,  parameter :: LMAX = 32
   type(sht_grid) :: sht
   logical :: ok
   integer :: i, j, lm
   real(wp) :: v, ref, err, c0
   real(wp), allocatable :: f(:,:), g(:,:)
   complex(wp), allocatable :: flm(:)

   ok = .true.
   call sht_grid_init(sht, LMAX, nlat=2*LMAX, nphi=4*LMAX)
   allocate(f(sht%nphi,sht%nlat), g(sht%nphi,sht%nlat), flm(sht%nlm))

   ! --- (1) eval_point normalization: pure Y00 is constant 1/sqrt(4pi) ----------
   flm = (0.0_wp, 0.0_wp);  flm(sht_grid_lmidx(sht, 0,0)) = (1.0_wp, 0.0_wp)
   c0  = 1.0_wp/sqrt(4.0_wp*pi)
   call sht_grid_eval_point(sht, flm, 0.7_wp, 1.3_wp, v)
   write(*,'(a,es12.4,a,es12.4)') ' (1) Y00 eval = ', v, '  expected ', c0
   if (abs(v - c0) > 1.0e-12_wp) then
      write(*,'(a)') '     FAIL: eval_point normalization off';  ok = .false.
   end if

   ! --- (2) eval_point at grid nodes reproduces the synthesized field -----------
   flm = (0.0_wp, 0.0_wp)
   flm(sht_grid_lmidx(sht, 2,0)) = ( 1.5_wp,  0.0_wp)
   flm(sht_grid_lmidx(sht, 3,1)) = (-0.7_wp,  0.4_wp)
   flm(sht_grid_lmidx(sht, 5,4)) = ( 0.3_wp, -0.9_wp)
   call sht_grid_synthesis(sht, flm, f)
   err = 0.0_wp
   do j = 1, sht%nlat
      do i = 1, sht%nphi
         call sht_grid_eval_point(sht, flm, sht%colat(j), sht%lon(i), v)
         err = max(err, abs(v - f(i,j)))
      end do
   end do
   write(*,'(a,es12.4)') ' (2) eval_point vs synthesis at grid nodes, max err = ', err
   if (err > 1.0e-11_wp) then
      write(*,'(a)') '     FAIL: eval_point disagrees with synthesis';  ok = .false.
   end if

   ! --- (3) spherical cap: pointwise-exact shape, zero outside, mass vs analytic -
   call spherical_cap(sht, 0.0_wp, 0.0_wp, 10.0_wp*pi/180.0_wp, 1500.0_wp, f)
   ! pointwise: the generated field must equal the analytic profile at every node
   ! (no grid node sits exactly at the pole, so this is the right correctness test)
   err = 0.0_wp
   do j = 1, sht%nlat
      do i = 1, sht%nphi
         err = max(err, abs(f(i,j) - cap_ref(sht%colat(j), sht%lon(i), &
                                              10.0_wp*pi/180.0_wp, 1500.0_wp)))
      end do
   end do
   write(*,'(a,es12.4)') ' (3) cap pointwise vs analytic shape, max err = ', err
   if (err > 1.0e-9_wp) then
      write(*,'(a)') '     FAIL: cap shape wrong';  ok = .false.
   end if
   if (maxval(f(:, sht%nlat/2:)) > 1.0e-9_wp) then
      write(*,'(a)') '     FAIL: cap nonzero far from centre (>90deg colat)';  ok = .false.
   end if
   ! cap mass converges to analytic with resolution (Gauss quadrature of a sharp
   ! feature): check on a finer grid to ~1%.
   call cap_mass_error(err)
   write(*,'(a,f6.2,a)') '     cap area-integral error at LMAX=160 = ', 100.0_wp*err, ' %'
   if (err > 1.5e-2_wp) then
      write(*,'(a)') '     FAIL: cap mass not converging to analytic';  ok = .false.
   end if

   ! --- (4) exponential basin: centre depth and far plateau ---------------------
   call exp_basin(sht, 35.0_wp*pi/180.0_wp, 25.0_wp*pi/180.0_wp, &
                  760.0_wp, 1200.0_wp, 26.0_wp*pi/180.0_wp, g)
   ! at the centre B = bmax - b0; far away (delta>>sigma) B -> bmax
   call eval_basin_centre(g, sht, 35.0_wp*pi/180.0_wp, 25.0_wp*pi/180.0_wp, v)
   write(*,'(a,f10.3,a)') ' (4) basin centre B = ', v, '  (expected -440)'
   if (abs(v - (760.0_wp-1200.0_wp)) > 5.0_wp) then
      write(*,'(a)') '     FAIL: basin centre wrong';  ok = .false.
   end if
   ! antipode of the basin centre -> essentially bmax
   ref = maxval(g)
   if (abs(ref - 760.0_wp) > 1.0_wp) then
      write(*,'(a)') '     FAIL: basin far-field plateau wrong';  ok = .false.
   end if

   ! --- (5) horizontal eval = surface gradient of the scalar field --------------
   ! eval_point_horiz(s_lm) must return (∂_θ f, (1/sinθ)∂_φ f) for f = Σ s_lm Y_lm.
   ! Validate convention-agnostically against a central finite-difference of the
   ! model's OWN scalar eval_point (no dependence on the real-SH normalization).
   block
      real(wp) :: vth, vph, fp, fm, dth, dph, hstep, colat0, lon0, sgrad_t, sgrad_p, eh
      flm = (0.0_wp, 0.0_wp)
      flm(sht_grid_lmidx(sht, 2,0)) = ( 1.5_wp,  0.0_wp)
      flm(sht_grid_lmidx(sht, 3,1)) = (-0.7_wp,  0.4_wp)
      flm(sht_grid_lmidx(sht, 5,4)) = ( 0.3_wp, -0.9_wp)
      colat0 = 1.0_wp;  lon0 = 2.0_wp;  hstep = 1.0e-4_wp   ! interior point, away from poles
      call sht_grid_eval_point_horiz(sht, flm, colat0, lon0, vth, vph)
      call sht_grid_eval_point(sht, flm, colat0+hstep, lon0, fp)
      call sht_grid_eval_point(sht, flm, colat0-hstep, lon0, fm)
      sgrad_t = (fp - fm)/(2.0_wp*hstep)                    ! ∂_θ f
      call sht_grid_eval_point(sht, flm, colat0, lon0+hstep, fp)
      call sht_grid_eval_point(sht, flm, colat0, lon0-hstep, fm)
      sgrad_p = (fp - fm)/(2.0_wp*hstep)/sin(colat0)        ! (1/sinθ) ∂_φ f
      dth = abs(vth - sgrad_t);  dph = abs(vph - sgrad_p)
      eh  = max(dth, dph)
      write(*,'(a,es12.4,a,es12.4)') ' (5) horiz vs FD-gradient: |dvth|=', dth, '  |dvph|=', dph
      if (eh > 1.0e-4_wp) then
         write(*,'(a)') '     FAIL: eval_point_horiz != surface gradient';  ok = .false.
      end if
   end block

   write(*,'(a)') ''
   if (ok) then
      write(*,'(a)') ' PASS: eval_point + analytic field generators validated'
   else
      write(*,'(a)') ' FAIL: field-infrastructure validation did not all pass'
      error stop 1
   end if
   call sht_grid_destroy(sht)

contains

   real(wp) function cap_ref(colat, lon, alpha, h0) result(h)
      !! Analytic spherical-cap height at a point (centre at the pole).
      real(wp), intent(in) :: colat, lon, alpha, h0
      real(wp) :: ca, d
      ca = cos(alpha);  d = angular_distance(colat, lon, 0.0_wp, 0.0_wp)
      if (d <= alpha) then
         h = h0*sqrt(max((cos(d)-ca)/(1.0_wp-ca), 0.0_wp))
      else
         h = 0.0_wp
      end if
   end function cap_ref

   subroutine cap_mass_error(err)
      !! Cap area-integral on a fine grid vs the analytic value (relative).
      real(wp), intent(out) :: err
      type(sht_grid) :: s2
      real(wp), allocatable :: ff(:,:)
      real(wp) :: gi, an
      call sht_grid_init(s2, 160, nlat=320, nphi=640)
      allocate(ff(s2%nphi, s2%nlat))
      call spherical_cap(s2, 0.0_wp, 0.0_wp, 10.0_wp*pi/180.0_wp, 1500.0_wp, ff)
      gi = sht_grid_surface_integral(s2, ff)
      an = cap_volume_over_R2(10.0_wp*pi/180.0_wp, 1500.0_wp)
      err = abs(gi - an)/an
      call sht_grid_destroy(s2)
   end subroutine cap_mass_error

   real(wp) function cap_volume_over_R2(alpha, h0) result(v)
      !! \int_cap h dOmega = 2*pi*h0 \int_0^alpha sqrt((cos d-ca)/(1-ca)) sin d dd.
      real(wp), intent(in) :: alpha, h0
      integer,  parameter  :: NQ = 200000
      real(wp) :: ca, dd, d, s, acc
      integer  :: k
      ca = cos(alpha);  dd = alpha/real(NQ,wp);  acc = 0.0_wp
      do k = 1, NQ
         d = (real(k,wp)-0.5_wp)*dd
         s = sqrt(max((cos(d)-ca)/(1.0_wp-ca), 0.0_wp))
         acc = acc + s*sin(d)*dd
      end do
      v = 2.0_wp*pi*h0*acc
   end function cap_volume_over_R2

   subroutine eval_basin_centre(g, sht, colat_c, lon_c, v)
      !! Nearest grid node to the basin centre (the analytic min is at the centre).
      real(wp),       intent(in)  :: g(:,:)
      type(sht_grid), intent(in)  :: sht
      real(wp),       intent(in)  :: colat_c, lon_c
      real(wp),       intent(out) :: v
      real(wp) :: dmin, d
      integer  :: i, j, ic, jc
      dmin = huge(1.0_wp);  ic = 1;  jc = 1
      do j = 1, sht%nlat
         do i = 1, sht%nphi
            d = angular_distance(sht%colat(j), sht%lon(i), colat_c, lon_c)
            if (d < dmin) then;  dmin = d;  ic = i;  jc = j;  end if
         end do
      end do
      v = g(ic,jc)
   end subroutine eval_basin_centre

end program test_field
