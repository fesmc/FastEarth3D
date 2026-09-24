program test_toroidal
   !! The toroidal degree of freedom under laterally varying viscosity
   !! (design-toroidal.md §4, V3 and V4). No external data.
   !!
   !! A held Y₂₀ load over M3-L70-V01 whose Maxwell layers carry a log10-viscosity
   !! perturbation p(θ,φ) = δ·sin²θ·cos 2φ (a Re Y₂₂ pattern):
   !!
   !!   (1) ALIVE — the toroidal field W is driven, far above round-off. A radially
   !!       symmetric Earth never forces W; this pattern must.
   !!   (2) SELECTION RULE — the configuration is symmetric under φ → −φ. The
   !!       spheroidal potentials are then even (cos mφ: real coefficients), the
   !!       toroidal potential ODD (sin mφ: imaginary coefficients), since e_r×∇T
   !!       picks up a sign under reflection. Load m = 0 times viscosity m = 0, ±2, ±4…
   !!       (10^p is not linear) populates even m only. All four statements are exact.
   !!   (3) SCALING — W appears at first order in the viscosity contrast, so halving
   !!       δ halves it (slope 1). Its feedback on the spheroidal field closes a loop
   !!       spheroidal → toroidal → spheroidal, second order: the uplift change it
   !!       causes falls 4× (slope 2). The slopes are perturbation reasoning, not
   !!       Martinec; their RATIO is the falsifiable part.
   !!   (4) ZERO (V3) — an axisymmetric perturbation p = δ·cos²θ under the same
   !!       axisymmetric load has every vertical plane as a mirror, which leaves the
   !!       toroidal potential nothing to be odd under: W ≡ 0 to round-off.
   use vilma_precision,       only: wp
   use vilma_constants,       only: kyr
   use vilma_earth_structure, only: earth_model, build_M3L70V01
   use vilma_radial_fe,       only: radial_fe_finalize
   use vilma_response,        only: response, response_init_ve, response_enable_lateral_visc, &
                                 response_begin_step, response_apply, response_commit_step, &
                                 response_horizontal_toroidal, response_destroy
   use vilma_sht,             only: sht_grid, sht_grid_init, sht_grid_destroy, sht_grid_lmidx
   implicit none

   integer,  parameter :: LMAX  = 12
   integer,  parameter :: NSTEP = 40
   real(wp), parameter :: DELTA = 0.2_wp          ! dex
   type(sht_grid)    :: sht
   type(earth_model) :: e
   real(wp) :: dt
   complex(wp), allocatable :: u_on(:), t_on(:), u_off(:), t_off(:), u_h(:), t_h(:), u_hoff(:), t_hoff(:)
   real(wp) :: wmax, umax, err, slope_w, slope_du, du1, du2
   integer  :: l, m, lm
   logical  :: ok

   ok = .true.
   dt = 0.04_wp*kyr                 ! 40 yr: M ≤ 0.1 in M3-L70's 1e21 mantle
   call sht_grid_init(sht, LMAX, nlat=3*LMAX, nphi=3*LMAX, mmax=LMAX)
   e = build_M3L70V01()
   allocate(u_on(sht%nlm), t_on(sht%nlm), u_off(sht%nlm), t_off(sht%nlm))
   allocate(u_h(sht%nlm), t_h(sht%nlm), u_hoff(sht%nlm), t_hoff(sht%nlm))

   ! --- (1) alive, (2) selection rule --------------------------------------------
   call run(DELTA, .false., .true., u_on, t_on)
   wmax = maxval(abs(t_on));  umax = maxval(abs(u_on))
   write(*,'(a,es11.3,a,es11.3,a)') ' (1) max|W(a)| =', wmax, ' m,  max|U(a)| =', umax, ' m'
   if (wmax < 1.0e-6_wp*umax) then
      write(*,'(a)') '     FAIL: the toroidal field is not driven'
      ok = .false.
   end if
   err = 0.0_wp
   do l = 1, LMAX
      do m = 0, l
         lm = sht_grid_lmidx(sht, l, m)
         if (mod(m,2) == 1) then                      ! odd m: nothing at all
            err = max(err, abs(u_on(lm))/umax, abs(t_on(lm))/wmax)
         else                                         ! even m: U real (cos), W imaginary (sin)
            err = max(err, abs(aimag(u_on(lm)))/umax, abs(real(t_on(lm), wp))/wmax)
         end if
      end do
   end do
   write(*,'(a,es10.2)') ' (2) reflection selection rule (U even/cos, W odd/sin, even m), worst =', err
   if (err > 1.0e-10_wp) then
      write(*,'(a)') '     FAIL: the toroidal field breaks the configuration''s symmetry'
      ok = .false.
   end if

   ! --- (3) scaling in δ -----------------------------------------------------------
   call run(DELTA,        .false., .false., u_off,  t_off)
   call run(0.5_wp*DELTA, .false., .true.,  u_h,    t_h)
   call run(0.5_wp*DELTA, .false., .false., u_hoff, t_hoff)
   slope_w  = log(maxval(abs(t_on))/maxval(abs(t_h)))/log(2.0_wp)
   du1 = maxval(abs(u_on - u_off));  du2 = maxval(abs(u_h - u_hoff))
   slope_du = log(du1/du2)/log(2.0_wp)
   write(*,'(a,f7.3,a)') ' (3) W(a) vs δ:                     slope', slope_w,  '  (expect 1)'
   write(*,'(a,f7.3,a,es10.2,a)') '     uplift change from W vs δ:     slope', slope_du, &
        '  (expect 2; |ΔU| =', du1, ' m)'
   if (maxval(abs(t_off)) /= 0.0_wp) then
      write(*,'(a)') '     FAIL: toroidal switched off, yet W was reported'
      ok = .false.
   end if
   if (abs(slope_w - 1.0_wp) > 0.1_wp .or. abs(slope_du - 2.0_wp) > 0.2_wp) then
      write(*,'(a)') '     FAIL: the toroidal field does not scale as first/second order in δ'
      ok = .false.
   end if

   ! --- (4) V3: axisymmetric perturbation, W ≡ 0 ----------------------------------
   call run(DELTA, .true., .true., u_on, t_on)
   err = maxval(abs(t_on))/maxval(abs(u_on))
   write(*,'(a,es10.2)') ' (4) axisymmetric perturbation: max|W(a)|/max|U(a)| =', err
   if (err > 1.0e-12_wp) then
      write(*,'(a)') '     FAIL: an axisymmetric configuration drove toroidal flow'
      ok = .false.
   end if

   write(*,'(a)') ''
   if (ok) then
      write(*,'(a)') ' PASS: lateral viscosity drives the toroidal field, with the'
      write(*,'(a)') '       symmetry, scaling and axisymmetric null it must have'
   else
      write(*,'(a)') ' FAIL: toroidal validation did not all pass'
      call sht_grid_destroy(sht);  call radial_fe_finalize()
      error stop 1
   end if
   call sht_grid_destroy(sht);  call radial_fe_finalize()

contains

   subroutine run(delta, axisym, tor, u_lm, t_lm)
      !! Hold a unit-ish Y₂₀ load for NSTEP steps over the perturbed Earth; return
      !! the surface uplift U(a) and toroidal W(a) coefficients after the last step.
      real(wp),    intent(in)  :: delta
      logical,     intent(in)  :: axisym    !! p = δ cos²θ instead of δ sin²θ cos 2φ
      logical,     intent(in)  :: tor       !! carry the toroidal degree of freedom
      complex(wp), intent(out) :: u_lm(:), t_lm(:)
      type(response) :: ve
      real(wp), allocatable :: pert(:,:,:)
      complex(wp), allocatable :: slm(:), n_lm(:)
      integer :: i, ie, ip, it
      real(wp) :: st, ct
      call response_init_ve(ve, e, sht, dt)
      ve%toroidal = tor
      allocate(pert(sht%nphi, sht%nlat, ve%ne));  pert = 0.0_wp
      do ie = 1, ve%ne
         if (ve%MkPerDt(ie) == 0.0_wp) cycle          ! Maxwell layers only
         do it = 1, sht%nlat
            st = sin(sht%colat(it));  ct = cos(sht%colat(it))
            do ip = 1, sht%nphi
               if (axisym) then
                  pert(ip,it,ie) = delta*ct*ct
               else
                  pert(ip,it,ie) = delta*st*st*cos(2.0_wp*sht%lon(ip))
               end if
            end do
         end do
      end do
      call response_enable_lateral_visc(ve, sht, pert)
      allocate(slm(sht%nlm), n_lm(sht%nlm));  slm = (0.0_wp, 0.0_wp)
      slm(sht_grid_lmidx(sht, 2, 0)) = (1.0e3_wp, 0.0_wp)
      do i = 1, NSTEP
         call response_begin_step(ve, sht)
         call response_apply(ve, sht, slm, u_lm, n_lm)
         call response_commit_step(ve, sht, slm)
      end do
      call response_begin_step(ve, sht)
      call response_apply(ve, sht, slm, u_lm, n_lm)
      call response_horizontal_toroidal(ve, sht, t_lm)
      call response_destroy(ve)
   end subroutine run

end program test_toroidal
