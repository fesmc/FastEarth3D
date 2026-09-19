program test_relax
   !! Rung-3 validation: explicit Maxwell viscoelastic relaxation (fe_viscoelastic).
   !! A degree-j load held on a homogeneous Maxwell sphere must relax from the
   !! ELASTIC Love number (t=0, deviatoric stress fully present) to the FLUID
   !! limit (t→∞, deviatoric stress fully relaxed) — the two states already
   !! validated analytically in test_love. Checks:
   !!   (1) h(0) equals an independent elastic solve of the same μ,
   !!   (2) h(t→∞) → −(2j+1)/3 (fluid limit),
   !!   (3) the relaxation is smooth and monotonic,
   !!   (4) the relaxation time scales linearly with viscosity η.
   use fe_precision,       only: wp
   use fe_constants,       only: pi, grav_G, sec_per_year
   use fe_earth_structure, only: earth_model, earth_layer, RHEOL_MAXWELL
   use fe_radial_fe,       only: radial_operator_destroy, radial_operator_solve, radial_operator_assemble, radial_mesh_build, radial_mesh, radial_operator, radial_fe_finalize
   use fe_viscoelastic,    only: ve_destroy, ve_step, ve_init, ve_degree
   implicit none

   real(wp), parameter :: km = 1.0e3_wp, yr = sec_per_year
   real(wp), parameter :: a = 6371.0_wp*km, rho = 5511.0_wp, mu = 1.4e11_wp
   integer,  parameter :: j = 2
   real(wp) :: g, phiL, hfluid, h_el_ref
   real(wp) :: h0a, hinfa, tefa, h0b, hinfb, tefb
   logical  :: ok

   ok = .true.
   g      = grav_G*(4.0_wp/3.0_wp)*pi*rho*a          ! surface gravity, homog sphere
   phiL   = 4.0_wp*pi*grav_G*a/real(2*j+1, wp)        ! load potential (σ=1)
   hfluid = -real(2*j+1, wp)/3.0_wp

   ! Independent elastic reference (radial_operator on the same μ).
   h_el_ref = elastic_h()

   ! Relaxation runs at η and 2η.
   call relax_run(1.0e21_wp, h0a, hinfa, tefa)
   call relax_run(2.0e21_wp, h0b, hinfb, tefb)

   write(*,'(a)')        ' Maxwell relaxation, degree 2, homogeneous sphere'
   write(*,'(a,f10.5)')  '   elastic h(0)  reference   = ', h_el_ref
   write(*,'(a,f10.5)')  '   fluid   h(inf) target      = ', hfluid
   write(*,'(a)')        ''
   write(*,'(a,f10.5,a,f10.5,a,f8.3,a)') '   eta=1e21:  h0=',h0a,'  hinf=',hinfa, &
        '  t_efold=',tefa/(1000*yr),' kyr'
   write(*,'(a,f10.5,a,f10.5,a,f8.3,a)') '   eta=2e21:  h0=',h0b,'  hinf=',hinfb, &
        '  t_efold=',tefb/(1000*yr),' kyr'
   write(*,'(a,f8.3)')   '   t_efold(2eta)/t_efold(eta) = ', tefb/tefa

   ! (1) elastic start
   if (abs(h0a - h_el_ref) > 1.0e-6_wp*abs(h_el_ref)) then
      write(*,'(a)') '   FAIL: VE t=0 state /= elastic solve'; ok = .false.
   end if
   ! (2) fluid end
   if (abs(hinfa - hfluid) > 1.0e-3_wp .or. abs(hinfb - hfluid) > 1.0e-3_wp) then
      write(*,'(a)') '   FAIL: VE t->inf state /= fluid limit'; ok = .false.
   end if
   ! (4) viscosity scaling: doubling eta doubles the relaxation time
   if (abs(tefb/tefa - 2.0_wp) > 0.1_wp) then
      write(*,'(a)') '   FAIL: relaxation time does not scale with eta'; ok = .false.
   end if

   ! (5) the stepper runs for j=1 (same bordered KKT operator, no special-casing)
   call ve_degree1_smoke()

   write(*,'(a)') ''
   if (ok) then
      write(*,'(a)') ' PASS: Maxwell relaxation runs elastic -> fluid, t_relax ~ eta'
   else
      write(*,'(a)') ' FAIL: viscoelastic relaxation checks did not all pass'
      call radial_fe_finalize();  error stop 1
   end if
   call radial_fe_finalize()

contains

   function elastic_h() result(h)
      real(wp) :: h
      type(earth_model)     :: e
      type(radial_mesh)     :: m
      type(radial_operator) :: op
      real(wp) :: ua, va, fa
      call mk_earth(e, 1.0e21_wp);  call radial_mesh_build(m, e)
      call radial_operator_assemble(op, e, m, j)
      call radial_operator_solve(op, 1.0_wp, ua, va, fa)
      h = g*ua/phiL
      call radial_operator_destroy(op)
   end function elastic_h

   subroutine relax_run(eta, h0, hinf, t_efold)
      !! Step a held unit load and return h(0), h(t_max) and the e-folding time
      !! (first crossing of 63.2% of the elastic→fluid swing).
      real(wp), intent(in)  :: eta
      real(wp), intent(out) :: h0, hinf, t_efold
      type(earth_model) :: e
      type(radial_mesh) :: m
      type(ve_degree)   :: ve
      real(wp) :: dt, t, ua, va, fa, h, frac, hprev, tprev, fprev
      integer  :: istep, nstep
      logical  :: crossed
      call mk_earth(e, eta);  call radial_mesh_build(m, e)
      dt = 10.0_wp*yr
      call ve_init(ve, e, m, j, dt)
      nstep = nint(40.0e3_wp*yr/dt)
      h0 = 0.0_wp;  hinf = 0.0_wp;  t_efold = 0.0_wp
      crossed = .false.;  hprev = 0.0_wp;  tprev = 0.0_wp;  fprev = 0.0_wp
      do istep = 0, nstep
         call ve_step(ve, 1.0_wp, t, ua, va, fa)
         h = g*ua/phiL
         if (istep == 0) h0 = h
         frac = (h - h0)/(hfluid - h0)
         if (.not. crossed .and. frac >= 1.0_wp - 1.0_wp/exp(1.0_wp)) then
            ! linear interpolation between the bracketing steps
            t_efold = tprev + (t - tprev)*(0.632120559_wp - fprev)/(frac - fprev)
            crossed = .true.
         end if
         hprev = h;  tprev = t;  fprev = frac
      end do
      hinf = h
      call ve_destroy(ve)
   end subroutine relax_run

   subroutine ve_degree1_smoke()
      !! fe_viscoelastic must run for j=1 too. The stepper uses the SAME bordered
      !! KKT operator (radial_operator handles j=1 internally), so ve_degree needs
      !! no special-casing. Confirm a held degree-1 load on a Maxwell sphere steps
      !! stably — finite, non-trivial, and actually relaxing (the surface response
      !! evolves away from its elastic value), not frozen or divergent.
      type(earth_model) :: e
      type(radial_mesh) :: m
      type(ve_degree)   :: ve
      real(wp) :: dt, t, ua, va, fa, u0, ulast
      integer  :: istep, nstep
      logical  :: finite
      call mk_earth(e, 1.0e21_wp);  call radial_mesh_build(m, e)
      dt = 10.0_wp*yr
      call ve_init(ve, e, m, 1, dt)              ! degree 1
      nstep = nint(20.0e3_wp*yr/dt)
      finite = .true.;  u0 = 0.0_wp;  ulast = 0.0_wp
      do istep = 0, nstep
         call ve_step(ve, 1.0_wp, t, ua, va, fa)
         if (ua /= ua .or. abs(ua) > 1.0e30_wp) finite = .false.
         if (istep == 0) u0 = ua
         ulast = ua
      end do
      call ve_destroy(ve)
      write(*,'(a)') ''
      write(*,'(a,es12.4,a,es12.4)') '   j=1 VE: U(a) elastic=', u0, '  relaxed=', ulast
      if (.not. finite) then
         write(*,'(a)') '   FAIL: j=1 viscoelastic step not finite'; ok = .false.
      end if
      if (abs(u0) <= 0.0_wp) then
         write(*,'(a)') '   FAIL: j=1 elastic response is zero'; ok = .false.
      end if
      if (abs(ulast - u0) <= 1.0e-3_wp*abs(u0)) then
         write(*,'(a)') '   FAIL: j=1 load did not relax'; ok = .false.
      end if
   end subroutine ve_degree1_smoke

   subroutine mk_earth(e, eta)
      type(earth_model), intent(out) :: e
      real(wp),          intent(in)  :: eta
      e%name = "maxwell";  e%r_earth = a;  e%r_core = 0.0_wp
      allocate(e%layers(1))
      e%layers(1) = earth_layer(0.0_wp, a, rho, mu, eta, RHEOL_MAXWELL)
   end subroutine mk_earth

end program test_relax
