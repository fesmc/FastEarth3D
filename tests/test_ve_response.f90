program test_ve_response
   !! Rung-4 viscoelastic field driver (ve_response). Two checks pin it to
   !! already-validated code with no new reference data:
   !!
   !!   (1) elastic limit — with zero memory (the first step) the affine drift is
   !!       zero, so the per-degree gains gu(l), gn(l) must equal the standalone
   !!       elastic_response gains to machine precision;
   !!   (2) 1-D stepper agreement — driving a single held coefficient through
   !!       begin/apply/commit must reproduce the validated 1-D ve_degree
   !!       displacement history U(a,t) at every degree, and the geoid history
   !!       F(a,t) for degrees >=2, confirming the per-(l,m) memory path is the
   !!       same Maxwell kernel. Degree 1 is special: its geoid is referenced to
   !!       the CM frame (N₁≡0), so the field driver's degree-1 geoid is zero by
   !!       design and does not track the 1-D stepper's raw F — checked as N₁=0.
   use fe_precision,       only: wp
   use fe_constants,       only: kyr
   use fe_earth_structure, only: earth_model, build_M3L70V01
   use fe_radial_fe,       only: radial_mesh_build, radial_mesh, radial_fe_finalize
   use fe_viscoelastic,    only: ve_destroy, ve_step, ve_init, ve_degree, SCHEME_FE, SCHEME_TRAP
   use fe_response,        only: response_commit_step, response_apply, response_begin_step, response_destroy, response, response_init_elastic, response_init_ve, response_init_null
   use fe_sht,             only: sht_grid, sht_grid_init, sht_grid_destroy, sht_grid_lmidx
   implicit none

   integer, parameter :: LMAX = 8
   type(sht_grid)     :: sht
   type(earth_model)  :: e
   real(wp) :: dt
   logical  :: ok

   ok = .true.
   dt = 0.02_wp*kyr                 ! 20 yr explicit step (VEGA)
   call sht_grid_init(sht, LMAX, nlat=2*LMAX, nphi=4*LMAX)
   e = build_M3L70V01()

   write(*,'(a)') ' (1) elastic limit: ve_response gains == elastic_response'
   call elastic_limit(ok)

   write(*,'(a)') ''
   write(*,'(a)') ' (2) held single-degree load reproduces the 1-D ve_degree stepper (FE)'
   write(*,'(a)') '     j = 1 (geocenter, CM frame):'
   call stepper_agreement(1, SCHEME_FE, 1, ok)
   write(*,'(a)') ''
   write(*,'(a)') '     j = 2:'
   call stepper_agreement(2, SCHEME_FE, 1, ok)

   write(*,'(a)') ''
   write(*,'(a)') ' (3) TRAP field path reproduces the 1-D ve_degree TRAP (implicit coupling)'
   write(*,'(a)') '     j = 1 (geocenter, CM frame):'
   call stepper_agreement(1, SCHEME_TRAP, 50, ok)
   write(*,'(a)') ''
   write(*,'(a)') '     j = 2:'
   call stepper_agreement(2, SCHEME_TRAP, 50, ok)

   write(*,'(a)') ''
   if (ok) then
      write(*,'(a)') ' PASS: ve_response matches the elastic gains and the 1-D'
      write(*,'(a)') '       Maxwell stepper through the per-(l,m) memory path'
   else
      write(*,'(a)') ' FAIL: ve_response validation did not all pass'
      call radial_fe_finalize()
      error stop 1
   end if
   call sht_grid_destroy(sht)
   call radial_fe_finalize()

contains

   subroutine elastic_limit(ok)
      logical, intent(inout) :: ok
      type(response)      :: ve
      type(response) :: el
      integer  :: l
      real(wp) :: du, dn, dmax
      call response_init_ve(ve, e, sht, dt)
      call response_init_elastic(el, e, lmax=LMAX)
      ! The field driver carries degrees l>=1 (degree 1 in the CM frame via the
      ! sparse KKT border; degree 0 is the monopole geoid only), so compare the
      ! gains over the full deforming range.
      dmax = 0.0_wp
      do l = 1, LMAX
         du = abs(ve%gu(l) - el%ugain(l))
         dn = abs(ve%gn(l) - el%ngain(l))
         dmax = max(dmax, du, dn)
      end do
      write(*,'(a,es11.2)') '      max |gain difference| (l>=1) =', dmax
      if (dmax > 1.0e-12_wp) then
         write(*,'(a)') '      FAIL: ve gains differ from the elastic response'
         ok = .false.
      end if
      if (ve%gu(1) == 0.0_wp) then
         write(*,'(a)') '      FAIL: degree-1 deformation should be nonzero'
         ok = .false.
      end if
      call response_destroy(ve);  call response_destroy(el)
   end subroutine elastic_limit

   subroutine stepper_agreement(j, scheme, max_iter, ok)
      integer, intent(in)    :: j, scheme, max_iter
      logical, intent(inout) :: ok
      type(response)  :: ve
      type(ve_degree)    :: vd
      type(radial_mesh)  :: mesh
      complex(wp), allocatable :: slm(:), ulm(:), nlm(:)
      integer, parameter :: NSTEP = 25
      real(wp), parameter :: TOL = 1.0e-10_wp     ! coupling tol (both drivers; tight for the A/B)
      integer  :: lm, i
      real(wp) :: sigma, t1, ua1, va1, fa1, ua2, fa2, g
      real(wp) :: emax_u, emax_f, ref_u, nmax1

      sigma = 1.0_wp
      ! reference: validated 1-D stepper at degree j, same scheme/coupling settings
      call radial_mesh_build(mesh, e)
      call ve_init(vd, e, mesh, j=j, dt=dt)
      vd%scheme = scheme;  vd%max_couple_iter = max_iter;  vd%couple_tol = TOL

      ! field driver, single held (l=j,m=0) coefficient, same scheme/coupling settings
      call response_init_ve(ve, e, sht, dt)
      ve%scheme = scheme;  ve%max_couple_iter = max_iter;  ve%couple_tol = TOL
      g  = ve%g
      lm = sht_grid_lmidx(sht, j, 0)
      allocate(slm(sht%nlm), ulm(sht%nlm), nlm(sht%nlm))
      slm = (0.0_wp, 0.0_wp);  slm(lm) = cmplx(sigma, 0.0_wp, wp)

      emax_u = 0.0_wp;  emax_f = 0.0_wp;  ref_u = 0.0_wp;  nmax1 = 0.0_wp
      write(*,'(a)') '       step   U_a(1D)      U_a(field)     F_a(1D)      F_a(field)'
      do i = 1, NSTEP
         call ve_step(vd, sigma, t1, ua1, va1, fa1)          ! 1-D: state then advance
         call response_begin_step(ve, sht)
         call response_apply(ve, sht, slm, ulm, nlm)               ! field: state at this t
         ua2 = real(ulm(lm), wp)
         fa2 = -real(nlm(lm), wp)*g                       ! N = -F/g  ->  F = -N g
         call response_commit_step(ve, sht, slm)                    ! advance memory
         emax_u = max(emax_u, abs(ua1 - ua2))
         ref_u  = max(ref_u, abs(ua1))
         if (j >= 2) then
            emax_f = max(emax_f, abs(fa1 - fa2))          ! F comparison (deg>=2)
         else
            nmax1  = max(nmax1, abs(real(nlm(lm), wp)))   ! deg-1 geoid must stay 0
         end if
         if (mod(i,5) == 1 .or. i == NSTEP) &
            write(*,'(i9,4es14.5)') i, ua1, ua2, fa1, fa2
         ! imaginary channel must stay exactly zero
         if (abs(aimag(ulm(lm))) > 1.0e-14_wp) then
            write(*,'(a)') '      FAIL: imaginary channel leaked'; ok = .false.
         end if
      end do
      ! Displacement: the field driver must reproduce the 1-D ve_degree stepper at
      ! every degree (including the geocenter j=1, CE-like gauge).
      write(*,'(a,es11.2)')          '      relative U_a error =', emax_u/ref_u
      if (emax_u/ref_u > 1.0e-6_wp) then
         write(*,'(a)') '      FAIL: field driver displacement disagrees with the 1-D stepper'
         ok = .false.
      end if
      if (j >= 2) then
         ! Geoid: degrees >=2 must match the 1-D stepper's F_a.
         write(*,'(a,es11.2)')       '      max|ΔF_a| =', emax_f
         if (emax_f/abs(fa1) > 1.0e-6_wp) then
            write(*,'(a)') '      FAIL: field driver geoid disagrees with the 1-D stepper'
            ok = .false.
         end if
      else
         ! Geoid: degree-1 is referenced to the CM frame, so N₁ (hence F via the
         ! field driver) is identically zero by design — it deliberately does NOT
         ! track the 1-D stepper's raw F_a. Assert it stays zero.
         write(*,'(a,es11.2,a)')     '      max|N_1(field)| =', nmax1, '  (CM frame: must be 0)'
         if (nmax1 > 1.0e-30_wp) then
            write(*,'(a)') '      FAIL: degree-1 geoid not zero (CM-frame convention violated)'
            ok = .false.
         end if
      end if
      call ve_destroy(vd);  call response_destroy(ve)
   end subroutine stepper_agreement

end program test_ve_response
