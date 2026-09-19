program test_etd1
   !! Characterization of the ETD1 exponential memory integrator vs forward-Euler
   !! (fe_viscoelastic, the 1-D stepper). Same homogeneous Maxwell sphere + held
   !! degree-2 load as test_relax. This documents a DELIBERATE NEGATIVE RESULT: ETD1
   !! was implemented and measured as a candidate larger-dt / adaptive-stepping core,
   !! and found NOT to help for this model class. Kept as the reproducible evidence
   !! for that conclusion (so neither ETD0 nor ETD1 is re-attempted blindly), and as
   !! the scheme-pluggable infrastructure for the levers that DO help (see below).
   !!
   !! WHAT WAS MEASURED (held degree-2 load, eta=1e21, t_relax~1.1 kyr; full table
   !! printed below; reference = each scheme at dt=0.1 yr, which agree to 2e-4 so the
   !! dt->0 limit is real and scheme-independent):
   !!
   !!   dt[yr]   M     FE err    ETD1 err   FE/ETD1
   !!     25   0.11   5.8e-3    4.6e-2      0.12     <- resolved regime: FE ~8x better
   !!    100   0.44   2.3e-2    1.6e-1      0.14
   !!    400   1.8    1.9e-1    4.4e-1      0.43
   !!   1000   4.4    4.1e-1    7.1e-1      0.58
   !!   2000   8.8    1.3e0     8.7e-1      1.50     <- FE now diverging; ETD1 bounded
   !!   4000  18      9.2e0     1.1e0       8.3
   !!
   !! TWO CONCLUSIONS (both decision-relevant; see doc/performance-assessment.md):
   !!  (1) ETD1 is LESS accurate than FE in the usable (resolved, M<1) regime. Both
   !!      are 1st-order in dt here; ETD1 just has a ~8x worse error constant. The
   !!      cause is the EXPLICIT strain<->memory coupling: the strain fed to the
   !!      memory update is computed from the PREVIOUS step's memory (lagged), so it
   !!      is only 1st-order accurate -> the coupling, not the memory integrator, is
   !!      the order bottleneck. ETD1's higher-order memory treatment is wasted, and
   !!      its exponential under-relaxes per step (forcing weight 2mu*M*phi1 < 2mu*M),
   !!      the SAME "wrong direction" that sank ETD0 (fastearth3d-exp-memory-finding).
   !!  (2) FE is PRACTICALLY unconditionally stable for this model: it stays finite to
   !!      M~35 (the elastic/self-gravity feedback damps the naive M<2 scalar limit).
   !!      It produces garbage above M~2, but resolving the kyr relaxation needs M<1
   !!      anyway, so FE's stability is not a practical constraint here. ETD1's only
   !!      genuine edge -- staying bounded past M~8 -- lands in the under-resolved
   !!      regime where the answer is meaningless regardless. ETD's unconditional
   !!      stability would matter only with a genuinely weak (low-eta) layer that
   !!      destabilizes FE without needing fine resolution; no benchmark model has one.
   !!
   !! IMPLICATION: for the adaptive-stepping goal the lever is NOT the memory
   !! integrator. It is (a) iterating the strain<->memory coupling to consistency
   !! (lifts the 1st-order lag), and (b) FE step-doubling for the local-error estimate
   !! (FE is already stable, so no exponential integrator is needed).
   !!
   !! This test PASSES when the kernel stays correct (FE/ETD1 share the dt->0 limit)
   !! and ETD1 stays bounded at every dt; it REPORTS the FE-vs-ETD1 accuracy ranking
   !! rather than asserting it, so an eventual coupling fix that lets ETD1 win shows
   !! up in the printout instead of breaking the build.
   use fe_precision,       only: wp
   use fe_constants,       only: pi, grav_G, sec_per_year
   use fe_earth_structure, only: earth_model, earth_layer, RHEOL_MAXWELL
   use fe_radial_fe,       only: radial_mesh_build, radial_mesh, radial_operator, radial_fe_finalize
   use fe_viscoelastic,    only: ve_destroy, ve_step, ve_init, ve_degree, SCHEME_FE, SCHEME_ETD1
   implicit none

   real(wp), parameter :: km = 1.0e3_wp, yr = sec_per_year
   real(wp), parameter :: a = 6371.0_wp*km, rho = 5511.0_wp, mu = 1.4e11_wp
   real(wp), parameter :: eta = 1.0e21_wp                 ! t_efold ~ 1.1 kyr
   integer,  parameter :: j = 2
   integer,  parameter :: ND = 4                          ! diagnostic times
   real(wp), parameter :: tdiag(ND) = [1.0_wp, 2.0_wp, 5.0_wp, 10.0_wp]  ! kyr
   real(wp), parameter :: T_END = 10.0_wp                 ! kyr (run length)

   ! dt sweep (yr). M = mu*dt/eta ~ 4.42e-3 * dt[yr] for this sphere, so FE goes
   ! unstable (M > 2) near dt ~ 452 yr; the last two entries straddle that wall.
   integer,  parameter :: NS = 10
   real(wp), parameter :: dtsweep(NS) = [25.0_wp, 50.0_wp, 100.0_wp, 200.0_wp, &
                                         400.0_wp, 600.0_wp, 1000.0_wp, &
                                         2000.0_wp, 4000.0_wp, 8000.0_wp]

   real(wp), parameter :: DT_REF = 0.1_wp                 ! reference step (yr); M~4.4e-4
   real(wp) :: g, phiL, hfluid
   real(wp) :: href(ND), hfe_ref(ND)                      ! dt->0 reference h(tdiag)
   real(wp) :: hfe(ND), het(ND), efe, eet, eest, consist
   logical  :: fe_ok, et_ok, ok
   integer  :: i
   real(wp) :: worst_ratio

   ok = .true.
   g      = grav_G*(4.0_wp/3.0_wp)*pi*rho*a
   phiL   = 4.0_wp*pi*grav_G*a/real(2*j+1, wp)
   hfluid = -real(2*j+1, wp)/3.0_wp

   ! dt->0 truth. Both schemes converge to the same time-exact (space-discrete)
   ! solution as dt->0, so run BOTH at a tiny step: ETD1@DT_REF is the reference
   ! (2nd order, ~1e-7), and FE@DT_REF must agree with it (cross-check that the
   ! reference is scheme-independent, i.e. a genuine limit, not an ETD1 artefact).
   ! NB measuring scheme error against FE@1yr -- as a naive test would -- is biased:
   ! FE's own truncation cancels for FE and inflates for ETD1.
   call run(SCHEME_ETD1, DT_REF, href,    et_ok, eest)
   call run(SCHEME_FE,   DT_REF, hfe_ref, fe_ok, eest)
   consist = maxval(abs(hfe_ref - href))
   write(*,'(a)')       ' ETD1 vs forward-Euler -- held degree-2 load, homogeneous Maxwell sphere'
   write(*,'(a,f8.4,a,f8.4,a,f8.4,a)') '   reference (dt=0.1 yr):  h(1kyr)=', href(1), &
        '  h(10kyr)=', href(ND), '   (fluid h_inf=', hfluid, ')'
   write(*,'(a,es10.2)') '   reference cross-check |FE-ETD1| at dt=0.1 yr (max) = ', consist
   if (consist > 5.0e-3_wp) then
      write(*,'(a)') '   FAIL: FE and ETD1 do not share a dt->0 limit'; ok = .false.
   end if

   ! Elastic->fluid swing sets the scale; a bounded scheme cannot err by more than
   ! roughly this. Use it as the ETD1 stability bound.
   write(*,'(a)') ''
   write(*,'(a)') '   dt[yr]      M     FE err    ETD1 err   FE/ETD1   ETD1 est   stable'
   write(*,'(a)') '   --------------------------------------------------------------------'
   worst_ratio = 0.0_wp        ! worst (largest) ETD1/FE accuracy ratio in the resolved regime
   do i = 1, NS
      call run(SCHEME_FE,   dtsweep(i), hfe, fe_ok, eest)
      call run(SCHEME_ETD1, dtsweep(i), het, et_ok, eest)
      efe = errmax(hfe, href, fe_ok)
      eet = errmax(het, href, et_ok)
      write(*,'(a,f7.1,f7.3,2es11.2,f9.2,es11.2,a,2a4)') '   ', dtsweep(i), &
           mu*dtsweep(i)*yr/eta, efe, eet, ratio(efe,eet,fe_ok), eest, '   ', &
           trim(yn(fe_ok)), trim(yn(et_ok))
      ! HARD GUARD 1 -- ETD1 must stay BOUNDED at every dt (unconditional stability):
      ! a finite response no larger than ~the full elastic->fluid swing.
      if (.not. et_ok .or. eet > 2.0_wp*abs(hfluid)) then
         write(*,'(a)') '   FAIL: ETD1 not bounded (lost unconditional stability)'; ok = .false.
      end if
      ! DOCUMENTED FINDING (reported, not asserted): in the resolved regime (M<1)
      ! FE is the more accurate scheme. Track the worst ratio for the summary.
      if (fe_ok .and. mu*dtsweep(i)*yr/eta < 1.0_wp) &
         worst_ratio = max(worst_ratio, eet/max(efe, tiny(1.0_wp)))
   end do
   write(*,'(a)') '   (FE/ETD1 < 1 means FE more accurate; "stable" = FE,ETD1 finite)'

   write(*,'(a)') ''
   write(*,'(a,f5.1,a)') '   FINDING: in the resolved regime (M<1) FE is ', worst_ratio, &
        'x more accurate than ETD1 -- the explicit strain<->memory coupling, not the'
   write(*,'(a)') '            memory integrator, is the 1st-order bottleneck (see header).'
   write(*,'(a)') '            ETD1 stays bounded for all dt, but only wins where the kyr'
   write(*,'(a)') '            relaxation is already under-resolved. Negative result: do not'
   write(*,'(a)') '            adopt ETD1; pursue coupling-iteration + FE step-doubling instead.'

   write(*,'(a)') ''
   if (ok) then
      write(*,'(a)') ' PASS: ETD1 kernel correct (shares dt->0 limit) and bounded at all dt;'
      write(*,'(a)') '       documents that FE is more accurate in the resolved regime.'
   else
      write(*,'(a)') ' FAIL: ETD1 characterization checks did not all pass'
      call radial_fe_finalize();  error stop 1
   end if
   call radial_fe_finalize()

contains

   subroutine run(scheme, dt_yr, hd, stable, err_est)
      !! Step a held unit load with the given scheme/dt; sample h at the diagnostic
      !! times and report the peak embedded error estimate over the run. `stable`
      !! is false if the response blows up (NaN or |h| >> fluid swing).
      integer,  intent(in)  :: scheme
      real(wp), intent(in)  :: dt_yr
      real(wp), intent(out) :: hd(ND)
      logical,  intent(out) :: stable
      real(wp), intent(out) :: err_est
      type(earth_model) :: e
      type(radial_mesh) :: m
      type(ve_degree)   :: ve
      real(wp) :: dt, t, ua, va, fa, h, tk
      integer  :: istep, nstep, kd
      call mk_earth(e)
      call radial_mesh_build(m, e)
      dt = dt_yr*yr
      call ve_init(ve, e, m, j, dt)
      ve%scheme = scheme
      nstep = nint(T_END*1.0e3_wp*yr/dt)
      hd = 0.0_wp;  stable = .true.;  err_est = 0.0_wp;  kd = 1
      do istep = 1, nstep
         call ve_step(ve, 1.0_wp, t, ua, va, fa)
         h  = g*ua/phiL
         tk = t/(1.0e3_wp*yr)
         err_est = max(err_est, ve%err_last)
         if (h /= h .or. abs(h) > 50.0_wp) stable = .false.   ! fluid swing ~ |hfluid| = 5/3
         ! sample at (or just past) each diagnostic time
         do while (kd <= ND .and. tk >= tdiag(kd) - 0.5_wp*dt_yr/1.0e3_wp)
            hd(kd) = h;  kd = kd + 1
         end do
      end do
      do while (kd <= ND)        ! fill any unreached diagnostics with the last value
         hd(kd) = h;  kd = kd + 1
      end do
      call ve_destroy(ve)
   end subroutine run

   pure function errmax(h, hr, stable) result(em)
      real(wp), intent(in) :: h(ND), hr(ND)
      logical,  intent(in) :: stable
      real(wp) :: em
      if (.not. stable) then
         em = huge(1.0_wp)
      else
         em = maxval(abs(h - hr))
      end if
   end function errmax

   pure function ratio(efe, eet, fe_ok) result(r)
      real(wp), intent(in) :: efe, eet
      logical,  intent(in) :: fe_ok
      real(wp) :: r
      if (.not. fe_ok) then
         r = -1.0_wp                 ! FE unstable: ratio undefined
      else
         r = efe/max(eet, tiny(1.0_wp))
      end if
   end function ratio

   pure function yn(ok) result(s)
      logical, intent(in) :: ok
      character(len=4) :: s
      if (ok) then;  s = '  Y'
      else;          s = '  N'
      end if
   end function yn

   subroutine mk_earth(e)
      type(earth_model), intent(out) :: e
      e%name = "maxwell";  e%r_earth = a;  e%r_core = 0.0_wp
      allocate(e%layers(1))
      e%layers(1) = earth_layer(0.0_wp, a, rho, mu, eta, RHEOL_MAXWELL)
   end subroutine mk_earth

end program test_etd1
