program test_modal
   !! Tests for fe_modal (RESP_MODAL).
   !!
   !! Part A — propagator foundation: modal_apply_p reproduces fe_viscoelastic's
   !! homogeneous (σ=0) ve_step relaxation bit-for-bit (same kernel, wrapped as a
   !! matrix-free linear map); the relaxation decays; the backward-Euler
   !! propagator (used by the eigensolve) is stable.
   !!
   !! Part B — modal solve: the extracted modes {τ_k, C^u_k} reconstruct the
   !! held-load step response. Pins the two analytic limits (elastic t=0, fluid
   !! t→∞) tightly and tracks the transient against ve_step.
   use fe_precision,       only: wp
   use fe_constants,       only: pi, grav_G, sec_per_year
   use fe_earth_structure, only: earth_model, earth_layer, RHEOL_MAXWELL, RHEOL_FLUID, &
                                 build_M3L70V01, earth_gravity_at
   use fe_radial_fe,       only: radial_mesh, radial_mesh_build, radial_fe_finalize, &
                                 radial_operator, radial_operator_assemble, &
                                 radial_operator_solve, radial_operator_destroy
   use fe_viscoelastic,    only: ve_degree, ve_init, ve_step, ve_destroy, NLAM, &
                                 SCHEME_FE, SCHEME_BE
   use fe_modal,           only: modal_degree, modal_degree_init, &
                                 modal_degree_destroy, modal_apply_p, &
                                 modal_spectrum, modal_solve, modal_spectrum_destroy
   implicit none

   real(wp), parameter :: km = 1.0e3_wp, yr = sec_per_year
   real(wp), parameter :: a = 6371.0_wp*km, rho = 5511.0_wp, mu = 1.4e11_wp
   real(wp), parameter :: eta = 1.0e21_wp
   integer,  parameter :: j = 2

   type(earth_model) :: e
   type(radial_mesh) :: m
   logical  :: ok

   ok = .true.
   call mk_earth(e)
   call radial_mesh_build(m, e)

   call part_a()
   call part_b()
   call part_c()

   call radial_fe_finalize()
   write(*,'(a)') ''
   if (ok) then
      write(*,'(a)') ' PASS: fe_modal propagator + modal solve checks'
   else
      write(*,'(a)') ' FAIL: fe_modal checks failed'
      error stop 1
   end if

contains

   subroutine part_a()
      type(ve_degree)   :: ref
      type(modal_degree):: md
      real(wp), allocatable :: t0A(:,:), t0B(:,:), t0C(:,:)
      real(wp), allocatable :: cA(:,:), cB(:,:), cC(:,:), oA(:,:), oB(:,:), oC(:,:)
      real(wp) :: dt, t, ua, va, fa, n0, maxdiff, mdone, mbe
      integer  :: ne, istep, nstep
      logical  :: finite
      dt = 10.0_wp*yr;  ne = m%ne
      allocate(t0A(NLAM,ne), t0B(NLAM,ne), t0C(NLAM,ne))
      allocate(cA(NLAM,ne), cB(NLAM,ne), cC(NLAM,ne))
      allocate(oA(NLAM,ne), oB(NLAM,ne), oC(NLAM,ne))

      call ve_init(ref, e, m, j, dt)
      do istep = 1, 50
         call ve_step(ref, 1.0_wp, t, ua, va, fa)
      end do
      t0A = ref%Am;  t0B = ref%Bm;  t0C = ref%Cm
      call ve_destroy(ref)
      n0 = norm3(t0A, t0B, t0C)

      call ve_init(ref, e, m, j, dt)
      ref%Am = t0A;  ref%Bm = t0B;  ref%Cm = t0C
      ref%Un_prev = 0.0_wp;  ref%Vn_prev = 0.0_wp;  ref%time = 0.0_wp
      call modal_degree_init(md, e, m, j, dt, SCHEME_FE)
      cA = t0A;  cB = t0B;  cC = t0C
      nstep = 100;  maxdiff = 0.0_wp;  finite = .true.
      do istep = 1, nstep
         call ve_step(ref, 0.0_wp, t, ua, va, fa)
         call modal_apply_p(md, cA, cB, cC, oA, oB, oC)
         cA = oA;  cB = oB;  cC = oC
         maxdiff = max(maxdiff, norm3(cA-ref%Am, cB-ref%Bm, cC-ref%Cm))
         if (any(cA /= cA)) finite = .false.
      end do
      mdone = norm3(cA, cB, cC)
      call ve_destroy(ref);  call modal_degree_destroy(md)

      write(*,'(a)')        ' Part A: propagator foundation (homogeneous sphere, degree 2)'
      write(*,'(a,es10.3)') '   FE prop vs ve_step max|diff|     = ', maxdiff
      write(*,'(a,es10.3,a,es10.3)') '   memory |tau0|=', n0, '  -> after 100 steps ', mdone
      if (.not. finite)               then; write(*,'(a)') '   FAIL: FE prop not finite';        ok=.false.; end if
      if (maxdiff > 1.0e-12_wp*n0)    then; write(*,'(a)') '   FAIL: FE prop /= ve_step(0)';      ok=.false.; end if
      if (mdone >= n0)                then; write(*,'(a)') '   FAIL: relaxation did not decay';   ok=.false.; end if

      call modal_degree_init(md, e, m, j, 1.0e3_wp*yr, SCHEME_BE)
      call modal_apply_p(md, t0A, t0B, t0C, oA, oB, oC)
      mbe = norm3(oA, oB, oC)
      finite = .not. any(oA /= oA)
      call modal_degree_destroy(md)
      write(*,'(a,es10.3)') '   BE prop (dt=1kyr) |P.tau0|       = ', mbe
      if (.not. finite)  then; write(*,'(a)') '   FAIL: BE prop not finite';        ok=.false.; end if
      if (mbe >= n0)     then; write(*,'(a)') '   FAIL: BE prop did not relax';      ok=.false.; end if
   end subroutine part_a

   subroutine part_b()
      type(modal_spectrum) :: spec, spec1
      type(radial_operator):: op
      type(ve_degree)      :: ref
      real(wp) :: g, phiL, hfluid, u_el, u_fl, ua, va, fa
      real(wp) :: dt, t, umod, uve, swing, maxrel, tdom, cdom
      real(wp) :: gu_err, fl_err, fl_err1
      integer  :: istep, nstep, k

      g      = grav_G*(4.0_wp/3.0_wp)*pi*rho*a
      phiL   = 4.0_wp*pi*grav_G*a/real(2*j+1, wp)
      hfluid = -real(2*j+1, wp)/3.0_wp
      u_fl   = hfluid*phiL/g                         ! analytic fluid surface uplift

      ! independent elastic reference
      call radial_operator_assemble(op, e, m, j)
      call radial_operator_solve(op, 1.0_wp, ua, va, fa)
      u_el = ua
      call radial_operator_destroy(op)

      ! --- modal solve, all significant modes ---
      call modal_solve(spec, e, m, j, n_modes=-1)
      ! --- modal solve, single dominant mode (convergence demo) ---
      call modal_solve(spec1, e, m, j, n_modes=1)

      gu_err = abs(spec%gu - u_el)/abs(u_el)
      fl_err  = abs((spec%gu  + sum(spec%Cu )) - u_fl)/abs(u_fl)
      fl_err1 = abs((spec1%gu + sum(spec1%Cu)) - u_fl)/abs(u_fl)

      ! dominant mode (largest |Cu|)
      tdom = 0.0_wp;  cdom = 0.0_wp
      do k = 1, spec%nmode
         if (abs(spec%Cu(k)) > cdom) then;  cdom = abs(spec%Cu(k));  tdom = spec%tau(k);  end if
      end do

      ! --- transient vs ve_step (small dt to limit the reference's time error) ---
      dt = 5.0_wp*yr
      call ve_init(ref, e, m, j, dt)
      nstep = nint(60.0e3_wp*yr/dt)
      swing = abs(u_fl - u_el);  maxrel = 0.0_wp
      do istep = 0, nstep
         call ve_step(ref, 1.0_wp, t, ua, va, fa)
         uve = ua
         umod = spec%gu
         do k = 1, spec%nmode
            umod = umod + spec%Cu(k)*(1.0_wp - exp(-t/spec%tau(k)))
         end do
         maxrel = max(maxrel, abs(umod - uve)/swing)
      end do
      call ve_destroy(ref)

      write(*,'(a)')        ''
      write(*,'(a)')        ' Part B: modal solve step response (homogeneous sphere, degree 2)'
      write(*,'(a,i0)')     '   modes kept (K=all)               = ', spec%nmode
      write(*,'(a,f8.3,a)') '   dominant relaxation time tau     = ', tdom/(1000*yr), ' kyr'
      write(*,'(a,es10.3)') '   elastic gain error |gu-u_el|/u_el= ', gu_err
      write(*,'(a,es10.3)') '   fluid-limit error  (K=all)       = ', fl_err
      write(*,'(a,es10.3)') '   fluid-limit error  (K=1)         = ', fl_err1
      write(*,'(a,es10.3)') '   transient max rel err vs ve_step = ', maxrel

      if (gu_err > 1.0e-6_wp)  then; write(*,'(a)') '   FAIL: elastic gain /= elastic solve'; ok=.false.; end if
      if (fl_err  > 2.0e-3_wp) then; write(*,'(a)') '   FAIL: modal sum misses fluid limit';  ok=.false.; end if
      if (maxrel  > 3.0e-2_wp) then; write(*,'(a)') '   FAIL: transient off vs ve_step';      ok=.false.; end if

      call modal_spectrum_destroy(spec);  call modal_spectrum_destroy(spec1)
   end subroutine part_b

   subroutine part_c()
      !! Multilayer M3-L70-V01 (elastic litho + 3 Maxwell mantle layers + fluid
      !! core): exercises the Maxwell masking and the across-layer Galerkin
      !! reduction. Compare the modal step response to ve_step at degree 2.
      type(earth_model)    :: em, em_fl
      type(radial_mesh)    :: mm
      type(modal_spectrum) :: spec
      type(radial_operator):: op
      type(ve_degree)      :: ref
      real(wp) :: dt, t, ua, va, fa, u_el, u_fl, umod, uve, swing, maxrel
      real(wp) :: gu_err, fl_err, tdom, cdom, wmin, wsum_err, w_off
      integer  :: jj, istep, nstep, k, ee
      jj = 2
      em = build_M3L70V01()
      call radial_mesh_build(mm, em)

      ! elastic (t=0) reference
      call radial_operator_assemble(op, em, mm, jj)
      call radial_operator_solve(op, 1.0_wp, ua, va, fa)
      u_el = ua
      call radial_operator_destroy(op)

      ! fluid (t→∞) reference: same earth with Maxwell layers relaxed to fluid
      ! (μ=0). The modal asymptote gu + Σ Cu must equal this exactly.
      em_fl = em
      do k = 1, size(em_fl%layers)
         if (em_fl%layers(k)%rheology == RHEOL_MAXWELL) then
            em_fl%layers(k)%mu = 0.0_wp
            em_fl%layers(k)%rheology = RHEOL_FLUID
         end if
      end do
      call radial_operator_assemble(op, em_fl, mm, jj)
      call radial_operator_solve(op, 1.0_wp, ua, va, fa)
      u_fl = ua
      call radial_operator_destroy(op)

      ! larger BE step for the slow mantle modes (conditioning only; τ exact)
      call modal_solve(spec, em, mm, jj, n_modes=-1, dt_be=5.0e3_wp*yr)

      gu_err = abs(spec%gu - u_el)/abs(u_el)
      fl_err = abs((spec%gu + sum(spec%Cu)) - u_fl)/abs(u_fl)
      tdom = 0.0_wp;  cdom = 0.0_wp
      do k = 1, spec%nmode
         if (abs(spec%Cu(k)) > cdom) then;  cdom = abs(spec%Cu(k));  tdom = spec%tau(k);  end if
      end do

      ! transient: modal reconstruction vs ve_step over a long window (dt small
      ! enough for ve_step accuracy; covers the dominant relaxation)
      dt = 25.0_wp*yr
      nstep = nint(400.0e3_wp*yr/dt)
      call ve_init(ref, em, mm, jj, dt)
      swing = abs(u_fl - u_el);  maxrel = 0.0_wp
      do istep = 0, nstep
         call ve_step(ref, 1.0_wp, t, ua, va, fa)
         uve = ua;  umod = spec%gu
         do k = 1, spec%nmode
            umod = umod + spec%Cu(k)*(1.0_wp - exp(-t/spec%tau(k)))
         end do
         maxrel = max(maxrel, abs(umod - uve)/swing)
      end do
      call ve_destroy(ref)

      write(*,'(a)')        ''
      write(*,'(a)')        ' Part C: modal solve step response (M3-L70-V01, degree 2)'
      write(*,'(a,i0)')     '   modes kept (K=all)               = ', spec%nmode
      write(*,'(a,f9.3,a)') '   dominant relaxation time tau     = ', tdom/(1000*yr), ' kyr'
      write(*,'(a,es10.3)') '   elastic gain error               = ', gu_err
      write(*,'(a,es10.3)') '   fluid-limit error (K=all)        = ', fl_err
      write(*,'(a,es10.3)') '   transient max rel err vs ve_step = ', maxrel

      if (gu_err > 1.0e-6_wp)  then; write(*,'(a)') '   FAIL: elastic gain /= elastic solve'; ok=.false.; end if
      if (fl_err  > 1.0e-2_wp) then; write(*,'(a)') '   FAIL: modal sum misses fluid limit';  ok=.false.; end if
      if (maxrel  > 3.0e-2_wp) then; write(*,'(a)') '   FAIL: transient off vs ve_step';      ok=.false.; end if

      ! per-mode radial strain-energy weights (depth profile for lateral η): each kept
      ! mode's column is normalized (Σ_e w=1), non-negative, and lives only in the Maxwell
      ! mantle (zero in the elastic lithosphere and fluid core, which carry no memory).
      wmin = minval(spec%w);  wsum_err = 0.0_wp;  w_off = 0.0_wp
      do k = 1, spec%nmode
         wsum_err = max(wsum_err, abs(sum(spec%w(:,k)) - 1.0_wp))
         do ee = 1, mm%ne
            if (em%layers(mm%elem_layer(ee))%rheology /= RHEOL_MAXWELL) &
               w_off = max(w_off, abs(spec%w(ee,k)))
         end do
      end do
      write(*,'(a,es10.3,a,es10.3,a,es10.3)') '   mode weights: min=', wmin, &
           '  |Σ-1|max=', wsum_err, '  off-mantle max=', w_off
      if (wmin < -1.0e-12_wp)    then; write(*,'(a)') '   FAIL: negative mode weight';       ok=.false.; end if
      if (wsum_err > 1.0e-10_wp) then; write(*,'(a)') '   FAIL: mode weight not normalized'; ok=.false.; end if
      if (w_off > 1.0e-12_wp)    then; write(*,'(a)') '   FAIL: mode weight leaks outside the Maxwell mantle'; ok=.false.; end if

      call modal_spectrum_destroy(spec)
   end subroutine part_c

   pure real(wp) function norm3(A, B, C) result(nrm)
      real(wp), intent(in) :: A(:,:), B(:,:), C(:,:)
      nrm = max(maxval(abs(A)), maxval(abs(B)), maxval(abs(C)))
   end function norm3

   subroutine mk_earth(em)
      type(earth_model), intent(out) :: em
      em%name = "maxwell";  em%r_earth = a;  em%r_core = 0.0_wp
      allocate(em%layers(1))
      em%layers(1) = earth_layer(0.0_wp, a, rho, mu, eta, RHEOL_MAXWELL)
   end subroutine mk_earth

end program test_modal
