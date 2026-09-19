program test_couple_order
   !! Convergence-order characterization of the strain<->memory COUPLING (the §3c
   !! lever) in the 1-D ve_degree stepper. Companion to test_etd1, which rejected
   !! the exponential memory integrators and concluded that the order bottleneck is
   !! the EXPLICIT coupling, not the memory rule. This test measures that directly.
   !!
   !! Same homogeneous Maxwell sphere + held degree-2 load as test_relax/test_etd1
   !! (eta=1e21, t_efold~1.1 kyr). For this sphere M = mu*dt*yr/eta ~ 4.42e-3 * dt[yr],
   !! so the dt sweep below (5..80 yr) stays in the resolved regime M < 0.4 where the
   !! asymptotic order is clean.
   !!
   !! FOUR variants that separate "iterate the coupling" from "raise the order":
   !!   A  FE   explicit          -- the historical scheme (1st-order memory rule)
   !!   B  BE   iterated          -- implicit 1st-order rule, coupling iterated to
   !!                                consistency: the control for "iteration alone"
   !!   C  TRAP, 1 endpoint iter  -- trapezoid rule, single predictor (no fixed point)
   !!   D  TRAP, iterated         -- trapezoid rule, coupling iterated to consistency
   !!
   !! WHAT IT ESTABLISHES (correcting the §3c hypothesis that "iterate the coupling"
   !! is itself the order-lifting lever):
   !!   * The observable's order equals the order of the MEMORY time-integration --
   !!     the balance solve is algebraically exact given tau. So iterating the
   !!     coupling to consistency with a 1st-order rule (B) stays 1st-order: the
   !!     iteration is NOT an independent order lever.
   !!   * 2nd order needs a 2nd-order memory rule (trapezoidal). That rule is IMPLICIT
   !!     in the endpoint strain, so it REQUIRES the coupling iteration to solve (D).
   !!     The two §3c bullets are one lever: an implicit 2nd-order advance solved by
   !!     coupling iteration. The fixed point converges in a handful of iterations.
   !!   * The report is always the balance against the ENTERING memory tau_n (time-
   !!     aligned with the sample time); only the carried-forward memory is iterated.
   !!
   !! The reference is the most accurate scheme (D) at a tiny dt, cross-checked
   !! against FE single-pass at the same tiny dt to confirm a shared dt->0 limit
   !! (so the measured "error" is a genuine truncation error, not a scheme artefact).
   use fe_precision,       only: wp
   use fe_constants,       only: pi, grav_G, sec_per_year
   use fe_earth_structure, only: earth_model, earth_layer, RHEOL_MAXWELL
   use fe_radial_fe,       only: radial_mesh_build, radial_mesh, radial_operator, radial_fe_finalize
   use fe_viscoelastic,    only: ve_step_double, ve_destroy, ve_step, ve_init, ve_degree, SCHEME_FE, SCHEME_TRAP, SCHEME_BE
   implicit none

   real(wp), parameter :: km = 1.0e3_wp, yr = sec_per_year
   real(wp), parameter :: a = 6371.0_wp*km, rho = 5511.0_wp, mu = 1.4e11_wp
   real(wp), parameter :: eta = 1.0e21_wp                 ! t_efold ~ 1.1 kyr
   integer,  parameter :: j = 2
   integer,  parameter :: ND = 4                          ! diagnostic times
   real(wp), parameter :: tdiag(ND) = [1.0_wp, 2.0_wp, 5.0_wp, 10.0_wp]  ! kyr
   real(wp), parameter :: T_END = 10.0_wp                 ! kyr (run length)
   integer,  parameter :: MAXIT = 50                      ! coupling-iteration cap

   ! dt sweep (yr), each half the previous so order p = log2(e(dt)/e(dt/2)). All
   ! keep M < 0.4 (resolved regime where the asymptotic order is clean).
   integer,  parameter :: NS = 5
   real(wp), parameter :: dtsweep(NS) = [80.0_wp, 40.0_wp, 20.0_wp, 10.0_wp, 5.0_wp]
   real(wp), parameter :: DT_REF = 0.25_wp                ! reference step (yr); M~1.1e-3
   real(wp), parameter :: TOL_TIGHT = 1.0e-12_wp          ! coupling tol for clean order measurement

   real(wp) :: g, phiL, hfluid
   real(wp) :: href(ND), hfe_ref(ND), consist
   real(wp) :: hA(ND), hB(ND), hC(ND), hD(ND)
   real(wp) :: eA(NS), eB(NS), eC(NS), eD(NS)
   integer  :: itB(NS), itC(NS), itD(NS)
   real(wp) :: pA, pB, pC, pD
   integer  :: i, dummy
   logical  :: ok

   ok = .true.
   g      = grav_G*(4.0_wp/3.0_wp)*pi*rho*a
   phiL   = 4.0_wp*pi*grav_G*a/real(2*j+1, wp)
   hfluid = -real(2*j+1, wp)/3.0_wp

   ! dt->0 reference: D (most accurate) at a tiny step; cross-check vs FE single-pass.
   call run(SCHEME_TRAP, MAXIT, TOL_TIGHT, DT_REF, href,    dummy)
   call run(SCHEME_FE,   1,     TOL_TIGHT, DT_REF, hfe_ref, dummy)
   consist = maxval(abs(hfe_ref - href))
   write(*,'(a)')       ' Coupling-order characterization -- held degree-2 load, homogeneous Maxwell sphere'
   write(*,'(a,f8.4,a,f8.4,a,f8.4,a)') '   reference D (dt=0.25 yr): h(1kyr)=', href(1), &
        '  h(10kyr)=', href(ND), '   (fluid h_inf=', hfluid, ')'
   write(*,'(a,es10.2)') '   reference cross-check |FE-D| at dt=0.25 yr (max) = ', consist
   if (consist > 5.0e-3_wp) then
      write(*,'(a)') '   FAIL: schemes do not share a dt->0 limit'; ok = .false.
   end if

   ! Sweep all four variants.
   do i = 1, NS
      call run(SCHEME_FE,   1,     TOL_TIGHT, dtsweep(i), hA, dummy)
      call run(SCHEME_BE,   MAXIT, TOL_TIGHT, dtsweep(i), hB, itB(i))
      call run(SCHEME_TRAP, 1,     TOL_TIGHT, dtsweep(i), hC, itC(i))
      call run(SCHEME_TRAP, MAXIT, TOL_TIGHT, dtsweep(i), hD, itD(i))
      eA(i) = maxval(abs(hA - href))
      eB(i) = maxval(abs(hB - href))
      eC(i) = maxval(abs(hC - href))
      eD(i) = maxval(abs(hD - href))
   end do

   write(*,'(a)') ''
   write(*,'(a)') '   variant: A=FE  B=BE/iter  C=TRAP/1it  D=TRAP/iter'
   write(*,'(a)') '   dt[yr]    M     err A      err B      err C      err D    (B,D iters)'
   write(*,'(a)') '   --------------------------------------------------------------------------'
   do i = 1, NS
      write(*,'(a,f6.1,f7.3,4es11.2,a,i3,a,i3,a)') '   ', dtsweep(i), &
           mu*dtsweep(i)*yr/eta, eA(i), eB(i), eC(i), eD(i), &
           '   (', itB(i), ',', itD(i), ')'
   end do

   ! Observed order between the two finest steps (cleanest asymptotic estimate).
   pA = order(eA(NS-1), eA(NS))
   pB = order(eB(NS-1), eB(NS))
   pC = order(eC(NS-1), eC(NS))
   pD = order(eD(NS-1), eD(NS))
   write(*,'(a)') ''
   write(*,'(a,4(f6.2))') '   observed order (finest pair)   A,B,C,D = ', pA, pB, pC, pD

   ! --- iteration-cost curve (TRAP iterated): iters & error vs couple_tol --------
   ! As the coupling tol tightens, iterations rise but the error plateaus once the
   ! fixed point is resolved below the trapezoidal truncation floor. The knee is the
   ! practical tol: tightening past it buys iterations, not accuracy.
   call cost_curve()

   ! --- step-doubling local-error estimate (§3c part ii) -------------------------
   ! The Richardson estimate of an order-p method's LOCAL error scales as dt^(p+1):
   ! ~dt^3 for trapezoidal (p=2), ~dt^2 for forward-Euler (p=1). Confirms the
   ! estimate is a faithful, order-aware accept/reject signal for an adaptive dt.
   call step_doubling_check()

   write(*,'(a)') ''
   write(*,'(a)') '   FINDING: iterating the coupling to consistency with a 1st-order rule (B)'
   write(*,'(a)') '            stays 1st order -- the iteration is not itself an order lever.'
   write(*,'(a)') '            2nd order comes from the trapezoidal RULE (D), which is implicit'
   write(*,'(a)') '            and so NEEDS the iteration to solve. One lever, not two (doc/§3c).'

   ! --- guards (the decision-relevant claims) ---------------------------------
   ! D reaches ~2nd order in the resolved regime.
   if (pD < 1.7_wp) then
      write(*,'(a,f5.2)') '   FAIL: TRAP iterated did not reach 2nd order, p=', pD; ok = .false.
   end if
   ! A (historical) is ~1st order.
   if (pA < 0.7_wp .or. pA > 1.4_wp) then
      write(*,'(a,f5.2)') '   FAIL: FE not ~1st order, p=', pA; ok = .false.
   end if
   ! Iterating the coupling with a 1st-order rule (B) stays 1st order -- the control.
   if (pB < 0.7_wp .or. pB > 1.4_wp) then
      write(*,'(a,f5.2)') '   FAIL: BE iterated not ~1st order (iteration is not an order lever), p=', pB
      ok = .false.
   end if
   ! D is strictly more accurate than A at the finest dt (the practical payoff).
   if (eD(NS) >= eA(NS)) then
      write(*,'(a)') '   FAIL: TRAP iterated not more accurate than FE'; ok = .false.
   end if
   ! The coupling fixed point is cheap (converges well inside the cap).
   if (maxval(itD) >= MAXIT) then
      write(*,'(a)') '   FAIL: coupling iteration hit its cap (did not converge)'; ok = .false.
   end if

   write(*,'(a)') ''
   if (ok) then
      write(*,'(a)') ' PASS: trapezoidal + iterated coupling is 2nd-order and cheap; iterating'
      write(*,'(a)') '       the coupling with a 1st-order rule stays 1st-order (one lever, not two).'
   else
      write(*,'(a)') ' FAIL: coupling-order characterization did not all pass'
      call radial_fe_finalize();  error stop 1
   end if
   call radial_fe_finalize()

contains

   subroutine run(scheme, max_iter, tol, dt_yr, hd, iters_max)
      !! Step a held unit load with the given memory scheme, coupling-iteration cap,
      !! and coupling tolerance; sample h at the diagnostic times and report the peak
      !! coupling-iteration count over the run.
      integer,  intent(in)  :: scheme, max_iter
      real(wp), intent(in)  :: tol, dt_yr
      real(wp), intent(out) :: hd(ND)
      integer,  intent(out) :: iters_max
      type(earth_model) :: e
      type(radial_mesh) :: m
      type(ve_degree)   :: ve
      real(wp) :: dt, t, ua, va, fa, h, tk
      integer  :: istep, nstep, kd
      call mk_earth(e)
      call radial_mesh_build(m, e)
      dt = dt_yr*yr
      call ve_init(ve, e, m, j, dt)
      ve%scheme          = scheme
      ve%max_couple_iter = max_iter
      ve%couple_tol      = tol
      nstep = nint(T_END*1.0e3_wp*yr/dt)
      hd = 0.0_wp;  kd = 1;  iters_max = 0
      do istep = 1, nstep
         call ve_step(ve, 1.0_wp, t, ua, va, fa)
         h  = g*ua/phiL
         tk = t/(1.0e3_wp*yr)
         iters_max = max(iters_max, ve%couple_iters_last)
         do while (kd <= ND .and. tk >= tdiag(kd) - 0.5_wp*dt_yr/1.0e3_wp)
            hd(kd) = h;  kd = kd + 1
         end do
      end do
      do while (kd <= ND)        ! fill any unreached diagnostics with the last value
         hd(kd) = h;  kd = kd + 1
      end do
      call ve_destroy(ve)
   end subroutine run

   subroutine cost_curve()
      !! TRAP iterated: peak coupling iterations and error at a few (dt, couple_tol)
      !! pairs. Shows the accuracy plateau (trapezoidal truncation floor) and how few
      !! iterations a practical 1e-6 tol needs vs the 1e-12 used for order measurement.
      integer,  parameter :: NDT = 3, NTOL = 4
      real(wp), parameter :: dtc(NDT)  = [80.0_wp, 20.0_wp, 5.0_wp]
      real(wp), parameter :: tolc(NTOL) = [1.0e-3_wp, 1.0e-6_wp, 1.0e-9_wp, 1.0e-12_wp]
      real(wp) :: hh(ND), err
      integer  :: id, it, iters
      write(*,'(a)') ''
      write(*,'(a)') '   Iteration-cost curve (TRAP iterated): peak iters / error vs couple_tol'
      write(*,'(a)') '   dt[yr]   M     tol=1e-3      1e-6        1e-9        1e-12'
      write(*,'(a)') '   ------------------------------------------------------------------------'
      do id = 1, NDT
         write(*,'(a,f6.1,f7.3,a)',advance='no') '   ', dtc(id), mu*dtc(id)*yr/eta, '  '
         do it = 1, NTOL
            call run(SCHEME_TRAP, MAXIT, tolc(it), dtc(id), hh, iters)
            err = maxval(abs(hh - href))
            write(*,'(a,i2,a,es8.1)',advance='no') ' ', iters, '/', err
         end do
         write(*,'(a)') ''
      end do
      write(*,'(a)') '   (cells are peak-iters/error; error stops improving once tol clears the floor)'
   end subroutine cost_curve

   subroutine step_doubling_check()
      !! Validate ve%step_double: the local-error estimate of the first step from the
      !! relaxed (tau=0) start, swept over dt. For an order-p method it must scale as
      !! dt^(p+1) -- order ~3 for trapezoidal, ~2 for forward-Euler. Reported AND
      !! guarded (this is the §3c part-ii deliverable).
      real(wp) :: esdT(NS), esdF(NS), pT, pF
      integer  :: i
      write(*,'(a)') ''
      write(*,'(a)') '   Step-doubling local-error estimate vs dt (first step from tau=0)'
      write(*,'(a)') '   dt[yr]   M      est FE       est TRAP'
      write(*,'(a)') '   ---------------------------------------------'
      do i = 1, NS
         esdF(i) = one_estimate(SCHEME_FE,   1,     dtsweep(i))
         esdT(i) = one_estimate(SCHEME_TRAP, MAXIT, dtsweep(i))
         write(*,'(a,f6.1,f7.3,2es13.3)') '   ', dtsweep(i), &
              mu*dtsweep(i)*yr/eta, esdF(i), esdT(i)
      end do
      pF = order(esdF(NS-1), esdF(NS))     ! expect ~2 (FE local error ~ dt^2)
      pT = order(esdT(NS-1), esdT(NS))     ! expect ~3 (TRAP local error ~ dt^3)
      write(*,'(a,2f7.2)') '   est scaling order (p+1)  FE, TRAP = ', pF, pT
      if (pF < 1.6_wp .or. pF > 2.4_wp) then
         write(*,'(a,f5.2)') '   FAIL: FE step-doubling estimate not ~dt^2, p+1=', pF; ok = .false.
      end if
      if (pT < 2.6_wp .or. pT > 3.4_wp) then
         write(*,'(a,f5.2)') '   FAIL: TRAP step-doubling estimate not ~dt^3, p+1=', pT; ok = .false.
      end if
   end subroutine step_doubling_check

   function one_estimate(scheme, max_iter, dt_yr) result(est)
      !! One step-doubling estimate of the first step from a fresh (relaxed) stepper.
      integer,  intent(in) :: scheme, max_iter
      real(wp), intent(in) :: dt_yr
      real(wp) :: est
      type(earth_model) :: e
      type(radial_mesh) :: m
      type(ve_degree)   :: ve
      call mk_earth(e)
      call radial_mesh_build(m, e)
      call ve_init(ve, e, m, j, dt_yr*yr)
      ve%scheme = scheme;  ve%max_couple_iter = max_iter;  ve%couple_tol = TOL_TIGHT
      call ve_step_double(ve, 1.0_wp, est)
      call ve_destroy(ve)
   end function one_estimate

   pure function order(e_coarse, e_fine) result(p)
      !! Observed convergence order from a dt -> dt/2 error pair: p = log2(e/e_half).
      real(wp), intent(in) :: e_coarse, e_fine
      real(wp) :: p
      p = log(e_coarse/max(e_fine, tiny(1.0_wp)))/log(2.0_wp)
   end function order

   subroutine mk_earth(e)
      type(earth_model), intent(out) :: e
      e%name = "maxwell";  e%r_earth = a;  e%r_core = 0.0_wp
      allocate(e%layers(1))
      e%layers(1) = earth_layer(0.0_wp, a, rho, mu, eta, RHEOL_MAXWELL)
   end subroutine mk_earth

end program test_couple_order
