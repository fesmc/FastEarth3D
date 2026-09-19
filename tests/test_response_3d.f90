program test_response_3d
   !! Rung 6a — laterally-varying viscosity, pseudo-spectral memory advance.
   !! Two self-checks pin the new 3D path to the validated 1-D field driver with
   !! no external data (a laterally-UNIFORM viscosity field must reduce the
   !! pseudo-spectral advance to the 1-D scalar advance, since the SHT round-trip
   !! is exact on a band-limited field):
   !!
   !!   (1) zero perturbation — enabling lateral viscosity with an all-zero log10
   !!       perturbation field must reproduce the ordinary 1-D ve_response memory
   !!       and uplift trajectory to SHT round-trip (machine) precision;
   !!   (2) uniform perturbation — a spatially-constant perturbation p must match a
   !!       1-D run on an Earth whose Maxwell-layer viscosities are scaled by 10^p,
   !!       validating the perturbation→Maxwell-rate mapping (η_eff = η·10^p).
   !!
   !! Both drive a multi-(l,m) load (degrees and orders, complex coefficients)
   !! through begin/apply/commit so the full spectral round-trip is exercised.
   use fe_precision,       only: wp
   use fe_constants,       only: kyr
   use fe_earth_structure, only: earth_n_layers, earth_model, build_M3L70V01, RHEOL_MAXWELL
   use fe_radial_fe,       only: radial_fe_finalize
   use fe_response,        only: response_apply, response_destroy, response_enable_lateral_visc, response_commit_step, response_begin_step, response, response_init_elastic, response_init_ve, response_init_null
   use fe_sht,             only: sht_grid, sht_grid_init, sht_grid_destroy, sht_grid_lmidx
   use fe_viscoelastic,    only: SCHEME_TRAP
   implicit none

   integer, parameter :: LMAX = 8
   type(sht_grid)     :: sht
   type(earth_model)  :: e
   real(wp) :: dt
   integer  :: nk1                  ! # degree-1 slots (contiguous at the start)
   logical  :: ok

   ok = .true.
   dt = 0.02_wp*kyr                 ! 20 yr explicit step (VEGA)
   ! Non-axisymmetric (mmax=lmax): the corrected lateral-viscosity advance reconstructs
   ! the strain/stress tensor on the grid via the general-order dyadic transforms. A
   ! uniform M still reduces to the 1-D advance per (l,m), so this exercises the full
   ! m>0 path. nlat/nphi = 3*lmax de-alias the spin-2 (G,H) dyadic channels.
   call sht_grid_init(sht, LMAX, nlat=3*LMAX, nphi=3*LMAX, mmax=LMAX)
   nk1 = min(1, sht%mmax) + 1        ! degree-1 orders m=0..min(1,mmax)
   e = build_M3L70V01()

   write(*,'(a)') ' (1) zero perturbation: lateral path == 1-D ve_response'
   call regression_zero(ok)

   write(*,'(a)') ''
   write(*,'(a)') ' (2) uniform perturbation p: lateral path == 1-D run with eta*10^p'
   call consistency_uniform(-0.4_wp, ok)

   write(*,'(a)') ''
   write(*,'(a)') ' (3) TRAP-3D: trapezoidal lateral path == 1-D trapezoidal run (eta*10^p)'
   call consistency_trap(-0.4_wp, ok)

   write(*,'(a)') ''
   if (ok) then
      write(*,'(a)') ' PASS: the pseudo-spectral memory advance reduces to the 1-D'
      write(*,'(a)') '       advance for a laterally-uniform viscosity field'
   else
      write(*,'(a)') ' FAIL: rung-6a lateral-viscosity validation did not all pass'
      call sht_grid_destroy(sht);  call radial_fe_finalize()
      error stop 1
   end if
   call sht_grid_destroy(sht);  call radial_fe_finalize()

contains

   subroutine build_load(slm)
      !! Deterministic NON-axisymmetric real-field load: every (l,m) up to LMAX
      !! populated (m=0 real; m>0 complex), so the general-order dyadic path and the
      !! m>0 advance are exercised.
      complex(wp), intent(out) :: slm(:)
      integer  :: l, m, lm
      slm = (0.0_wp, 0.0_wp)
      do l = 1, LMAX
         do m = 0, l
            lm = sht_grid_lmidx(sht, l, m)
            if (m == 0) then
               slm(lm) = cmplx(1000.0_wp/real(l*l, wp), 0.0_wp, wp)
            else
               slm(lm) = cmplx(700.0_wp/real(l*l, wp), 400.0_wp/real(l*(m+1), wp), wp)
            end if
         end do
      end do
   end subroutine build_load

   real(wp) function mem_diff(a, b) result(d)
      !! Max |a-b| over a memory array, masking the physically-inert λ=6 (local index
      !! 4) memory at the degree-1 slot (k=1 for mmax=0). Z⁶ has no degree-1 harmonic
      !! (its B13 norm is 0), so the dissipation never sees it; the 1-D scalar advance
      !! keeps a phantom value there while the tensor advance correctly zeros it. The
      !! difference is in the null space of the observable (the uplift check confirms).
      real(wp), intent(in) :: a(:,:,:), b(:,:,:)
      real(wp), allocatable :: da(:,:,:)
      da = a - b
      da(4,:,1:nk1) = 0.0_wp          ! λ=6 has no degree-1 harmonic (all its orders)
      d = maxval(abs(da))
   end function mem_diff

   real(wp) function mem_scale(a) result(s)
      real(wp), intent(in) :: a(:,:,:)
      s = maxval(abs(a))
   end function mem_scale

   subroutine drive_and_compare(ve3d, ve1d, label, ok)
      !! Drive both responses with the same load through NSTEP begin/apply/commit
      !! and report the max relative memory + uplift disagreement.
      type(response), intent(inout) :: ve3d, ve1d
      character(*),      intent(in)    :: label
      logical,           intent(inout) :: ok
      integer, parameter :: NSTEP = 25
      real(wp), parameter :: TOL = 1.0e-9_wp
      complex(wp), allocatable :: slm(:), u3(:), n3(:), u1(:), n1(:)
      integer  :: i
      real(wp) :: dmem, smem, dupl, supl

      allocate(slm(sht%nlm), u3(sht%nlm), n3(sht%nlm), u1(sht%nlm), n1(sht%nlm))
      call build_load(slm)
      dmem = 0.0_wp;  smem = 0.0_wp;  dupl = 0.0_wp;  supl = 0.0_wp
      do i = 1, NSTEP
         call response_begin_step(ve3d, sht);  call response_apply(ve3d, sht, slm, u3, n3)
         call response_begin_step(ve1d, sht);  call response_apply(ve1d, sht, slm, u1, n1)
         dupl = max(dupl, maxval(abs(u3 - u1)));  supl = max(supl, maxval(abs(u1)))
         call response_commit_step(ve3d, sht, slm)
         call response_commit_step(ve1d, sht, slm)
         dmem = max(dmem, mem_diff(ve3d%Are, ve1d%Are), mem_diff(ve3d%Bre, ve1d%Bre), &
                          mem_diff(ve3d%Cre, ve1d%Cre), mem_diff(ve3d%Aim, ve1d%Aim), &
                          mem_diff(ve3d%Bim, ve1d%Bim), mem_diff(ve3d%Cim, ve1d%Cim))
         smem = max(smem, mem_scale(ve1d%Are), mem_scale(ve1d%Bre), mem_scale(ve1d%Cre))
      end do
      write(*,'(a,a,a,es11.2,a,es11.2)') '      ', label, &
         ' rel mem err =', dmem/max(smem,tiny(1.0_wp)), &
         '   rel uplift err =', dupl/max(supl,tiny(1.0_wp))
      if (dmem/max(smem,tiny(1.0_wp)) > TOL .or. dupl/max(supl,tiny(1.0_wp)) > TOL) then
         write(*,'(a)') '      FAIL: lateral path diverges from the 1-D reference'
         ok = .false.
      end if
      deallocate(slm, u3, n3, u1, n1)
   end subroutine drive_and_compare

   subroutine regression_zero(ok)
      logical, intent(inout) :: ok
      type(response)     :: ve3d, ve1d
      real(wp), allocatable :: pert(:,:,:)
      call response_init_ve(ve1d, e, sht, dt)
      call response_init_ve(ve3d, e, sht, dt)
      allocate(pert(sht%nphi, sht%nlat, ve3d%ne));  pert = 0.0_wp
      ve3d%visc3d_tol = -1.0_wp     ! force every Maxwell element through the pseudo-spectral
      call response_enable_lateral_visc(ve3d, sht, pert)   ! kernel (else a uniform field collapses to 1-D)
      call drive_and_compare(ve3d, ve1d, 'pert=0 :', ok)
      deallocate(pert)
      call response_destroy(ve3d);  call response_destroy(ve1d)
   end subroutine regression_zero

   subroutine consistency_uniform(p, ok)
      real(wp), intent(in)    :: p
      logical,  intent(inout) :: ok
      type(response)     :: ve3d, ve1d
      type(earth_model)     :: es
      real(wp), allocatable :: pert(:,:,:)
      integer :: k
      ! 1-D reference Earth: Maxwell-layer viscosities scaled by 10^p (η_eff = η·10^p)
      es = build_M3L70V01()
      do k = 1, earth_n_layers(es)
         if (es%layers(k)%rheology == RHEOL_MAXWELL) &
            es%layers(k)%eta = es%layers(k)%eta * 10.0_wp**p
      end do
      call response_init_ve(ve1d, es, sht, dt)
      ! 3-D path: base Earth + a spatially-uniform perturbation p
      call response_init_ve(ve3d, e, sht, dt)
      allocate(pert(sht%nphi, sht%nlat, ve3d%ne));  pert = p
      ve3d%visc3d_tol = -1.0_wp     ! force the pseudo-spectral kernel (validate it reduces to 1-D)
      call response_enable_lateral_visc(ve3d, sht, pert)
      call drive_and_compare(ve3d, ve1d, 'uniform:', ok)
      deallocate(pert)
      call response_destroy(ve3d);  call response_destroy(ve1d)
   end subroutine consistency_uniform

   subroutine consistency_trap(p, ok)
      !! Same as consistency_uniform but with the implicit trapezoidal scheme
      !! (max_couple_iter>1): the TRAP-3D pseudo-spectral advance must reduce to the
      !! 1-D trapezoidal advance for a uniform field, to SHT round-trip precision.
      real(wp), intent(in)    :: p
      logical,  intent(inout) :: ok
      type(response)     :: ve3d, ve1d
      type(earth_model)     :: es
      real(wp), allocatable :: pert(:,:,:)
      integer :: k
      es = build_M3L70V01()
      do k = 1, earth_n_layers(es)
         if (es%layers(k)%rheology == RHEOL_MAXWELL) &
            es%layers(k)%eta = es%layers(k)%eta * 10.0_wp**p
      end do
      call response_init_ve(ve1d, es, sht, dt)
      call response_init_ve(ve3d, e, sht, dt)
      ve1d%scheme = SCHEME_TRAP;  ve1d%max_couple_iter = 4
      ve3d%scheme = SCHEME_TRAP;  ve3d%max_couple_iter = 4
      allocate(pert(sht%nphi, sht%nlat, ve3d%ne));  pert = p
      ve3d%visc3d_tol = -1.0_wp     ! force the pseudo-spectral kernel (validate it reduces to 1-D)
      call response_enable_lateral_visc(ve3d, sht, pert)
      call drive_and_compare(ve3d, ve1d, 'trap-uni:', ok)
      deallocate(pert)
      call response_destroy(ve3d);  call response_destroy(ve1d)
   end subroutine consistency_trap

end program test_response_3d
