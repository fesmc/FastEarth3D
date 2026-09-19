program test_benchmark_lvz
   !! Rung 6b — laterally-varying viscosity vs the ASPECT/Abaqus/TABOO 3D-code
   !! intercomparison of Weerdesteijn et al. (2023), G-cubed 24:e2022GC010813,
   !! Section 5.2 (short-timescale low-viscosity-region benchmark).
   !!
   !! Setup (their Fig 1 + Table 1 = Earth model M3-L70-V01, = build_M3L70V01):
   !!   - axisymmetric ice disc, radius 100 km, ρ_ice = 931 kg/m³, height H(t)
   !!     ramping 0 → 100 m over 100 yr, then held at 100 m to 200 yr;
   !!   - a cylindrical low-viscosity zone (LVZ) under the load: radius 100 km,
   !!     depth 70–170 km, η = 1e19 Pa·s (vs 1e21 upper mantle), μ unchanged.
   !!
   !! Published load-center maximum vertical surface deformation at t = 200 yr:
   !!   homogeneous (no LVZ) = −0.75 m;  with LVZ = −1.23 m (≈1.6× amplification
   !!   from the soft zone). ASPECT–Abaqus themselves agree to 1–3% (their Table 4),
   !!   and ASPECT (no self-gravity) vs TABOO (self-gravity, same M3 model) agree to
   !!   0.28% at the load center — so self-gravity is a sub-percent effect here and
   !!   the absolute numbers are a valid gate for our (self-gravitating, spherical)
   !!   model to a few percent.
   !!
   !! The load and the LVZ are both axisymmetric, so the whole problem is m=0:
   !! we run with mmax=0 (Legendre-only transforms) which keeps the pseudo-spectral
   !! memory advance cheap and lets us push lmax high enough to resolve the 100 km
   !! disc. We drive ve_response directly (pure ice-load deformation, no SLE/ocean).
   use fe_precision,       only: wp
   use fe_constants,       only: kyr
   use fe_earth_structure, only: earth_model, build_M3L70V01
   use fe_radial_fe,       only: radial_fe_finalize
   use fe_response,        only: response, response_init_elastic, response_init_ve, response_init_null, &
                                 response_apply, response_begin_step, response_commit_step, &
                                 response_enable_lateral_visc, response_destroy
   use fe_sht,             only: sht_grid, sht_grid_init, sht_grid_destroy, sht_grid_eval_point, sht_grid_lmidx
   implicit none

   integer,  parameter :: LMAX  = 512       ! axisymmetric; ~0.7° ≈ 78 km at the edge
   integer,  parameter :: NLATF = 3         ! nlat = NLATF*lmax (de-aliasing factor)
   real(wp), parameter :: YR    = 0.001_wp*kyr
   real(wp), parameter :: DT    = 2.5_wp*YR ! ASPECT's short-timescale step
   real(wp), parameter :: R_ICE = 100.0e3_wp
   real(wp), parameter :: RHO_ICE = 931.0_wp
   real(wp), parameter :: H_MAX = 100.0_wp
   real(wp), parameter :: T_RAMP = 100.0_wp*YR, T_END = 200.0_wp*YR
   real(wp), parameter :: LVZ_TOP = 70.0e3_wp, LVZ_BOT = 170.0e3_wp
   real(wp), parameter :: LVZ_PERT = -2.0_wp  ! log10(1e19/1e21): η 1e21 → 1e19
   real(wp), parameter :: U_HOMOG = -0.75_wp, U_LVZ = -1.23_wp  ! Weerdesteijn refs [m]
   real(wp), parameter :: TOL = 0.08_wp        ! 8% (codes themselves differ 1–3%)

   type(sht_grid)    :: sht
   type(earth_model) :: e
   real(wp) :: uc_homog, uc_lvz, re_homog, re_lvz
   logical  :: ok

   ok = .true.
   call sht_grid_init(sht, LMAX, nlat=NLATF*LMAX, nphi=2, mmax=0)
   e = build_M3L70V01()

   write(*,'(a,i0,a,f4.1,a)') ' Weerdesteijn 2023 LVZ benchmark (lmax=', LMAX, &
        ', dt=', DT/YR, ' yr, axisymmetric)'
   write(*,'(a)') ''
   call run_case(.false., uc_homog)
   call run_case(.true.,  uc_lvz)

   re_homog = abs(uc_homog - U_HOMOG)/abs(U_HOMOG)
   re_lvz   = abs(uc_lvz   - U_LVZ  )/abs(U_LVZ)
   write(*,'(a)') ''
   write(*,'(a)')               '   case          U_center(200 yr)   reference   rel.err'
   write(*,'(a,f12.4,a,f10.2,a,f8.1,a)') '   homogeneous ', uc_homog, ' m', U_HOMOG, ' m', 100*re_homog, ' %'
   write(*,'(a,f12.4,a,f10.2,a,f8.1,a)') '   with LVZ    ', uc_lvz,   ' m', U_LVZ,   ' m', 100*re_lvz,   ' %'
   write(*,'(a,f6.3,a)')        '   LVZ amplification U_lvz/U_homog =', uc_lvz/uc_homog, &
        '  (ref 1.64)'

   if (re_homog > TOL) then
      write(*,'(a)') '   FAIL: homogeneous disc uplift off the published -0.75 m'
      ok = .false.
   end if
   if (re_lvz > TOL) then
      write(*,'(a)') '   FAIL: LVZ disc uplift off the published -1.23 m'
      ok = .false.
   end if

   write(*,'(a)') ''
   if (ok) then
      write(*,'(a)') ' PASS: lateral-viscosity deformation matches Weerdesteijn 2023'
   else
      write(*,'(a)') ' FAIL: rung-6b LVZ benchmark out of tolerance'
      call sht_grid_destroy(sht);  call radial_fe_finalize()
      error stop 1
   end if
   call sht_grid_destroy(sht);  call radial_fe_finalize()

contains

   subroutine run_case(lvz_on, u_end)
      !! Drive the disc load 0→200 yr and return the load-center (pole) uplift at
      !! t = 200 yr. With lvz_on, the soft column is injected via enable_lateral_visc.
      logical,  intent(in)  :: lvz_on
      real(wp), intent(out) :: u_end
      type(response) :: ve
      complex(wp), allocatable :: disc_lm(:), slm(:), ulm(:), nlm(:)
      real(wp),    allocatable :: pert(:,:,:)
      real(wp) :: theta_c, t, H, rmid, depth, uval
      integer  :: i, ie, j, nstep

      call response_init_ve(ve, e, sht, DT)
      allocate(disc_lm(sht%nlm), slm(sht%nlm), ulm(sht%nlm), nlm(sht%nlm))

      ! Disc load → spectral via the EXACT spherical-cap coefficients (mass-correct,
      ! grid-independent), not by analysing a sampled step (whose quadrature error in
      ! the discontinuity makes the load mass — hence the uplift — depend on nlat).
      theta_c = R_ICE / e%r_earth
      call cap_coeffs(theta_c, disc_lm)

      if (lvz_on) then
         allocate(pert(sht%nphi, sht%nlat, ve%ne));  pert = 0.0_wp
         do ie = 1, ve%ne
            rmid  = 0.5_wp*(ve%r(ie) + ve%r(ie+1))
            depth = e%r_earth - rmid
            if (depth > LVZ_TOP .and. depth < LVZ_BOT) then
               do j = 1, sht%nlat
                  if (sht%colat(j) <= theta_c) pert(:,j,ie) = LVZ_PERT
               end do
            end if
         end do
         call response_enable_lateral_visc(ve, sht, pert)
         deallocate(pert)
      end if

      nstep = nint(T_END/DT)
      u_end = 0.0_wp
      do i = 1, nstep
         t = real(i,wp)*DT
         H = min(t/T_RAMP, 1.0_wp)*H_MAX           ! ramp then hold
         slm = RHO_ICE*H*disc_lm                   ! surface mass load [kg/m²]
         call response_begin_step(ve, sht)
         call response_apply(ve, sht, slm, ulm, nlm)
         call response_commit_step(ve, sht, slm)
         call sht_grid_eval_point(sht, ulm, 0.0_wp, 0.0_wp, uval)   ! pole = load center
         if (lvz_on .and. (mod(i,8)==0 .or. i==nstep)) &
            write(*,'(a,f7.1,a,f10.4)') '     t=', t/YR, ' yr  U_center=', uval
         if (i == nstep) u_end = uval
      end do

      deallocate(disc_lm, slm, ulm, nlm)
      call response_destroy(ve)
   end subroutine run_case

   subroutine cap_coeffs(tc, lm_coeffs)
      !! Exact m=0 spherical-harmonic coefficients (SHTns orthonormal) of the unit
      !! axisymmetric cap indicator (1 for θ<tc): with c=cos(tc),
      !!   a_0 = √π (1−c),   a_l = √(π/(2l+1)) [P_{l−1}(c) − P_{l+1}(c)],  l≥1,
      !! from ∫_c^1 P_l dx = [P_{l−1}−P_{l+1}]/(2l+1). Band-limited at lmax (the disc
      !! edge keeps Gibbs, but the load mass and the central response converge cleanly).
      real(wp),    intent(in)  :: tc
      complex(wp), intent(out) :: lm_coeffs(:)
      real(wp) :: c, p(0:LMAX+1)
      real(wp), parameter :: sqpi = 1.7724538509055160_wp   ! √π
      integer  :: l
      c = cos(tc)
      p(0) = 1.0_wp;  p(1) = c
      do l = 1, LMAX
         p(l+1) = (real(2*l+1,wp)*c*p(l) - real(l,wp)*p(l-1))/real(l+1,wp)
      end do
      lm_coeffs = (0.0_wp, 0.0_wp)
      lm_coeffs(sht_grid_lmidx(sht,0,0)) = cmplx(sqpi*(1.0_wp - c), 0.0_wp, wp)
      do l = 1, LMAX
         lm_coeffs(sht_grid_lmidx(sht,l,0)) = cmplx( &
            sqpi/sqrt(real(2*l+1,wp))*(p(l-1) - p(l+1)), 0.0_wp, wp)
      end do
   end subroutine cap_coeffs

end program test_benchmark_lvz
