program test_rotinv
   !! Rung 6c — off-pole rotational-invariance check (a strong end-to-end m>0 test
   !! of the laterally-varying-viscosity path, run with the implicit TRAP-3D advance).
   !!
   !! By rotational invariance of the (self-gravitating, spherical) physics, an ice
   !! cap sitting over a low-viscosity zone (LVZ) must produce the SAME peak vertical
   !! deformation at the cap centre regardless of where the cap+LVZ axis points. We
   !! drive an axisymmetric-about-n̂ load + LVZ for two axes — the north pole (m=0
   !! only) and an off-pole axis (all m up to lmax) — and require the cap-centre
   !! uplift to agree. The on-pole case is the validated axisymmetric path; matching
   !! it off-pole exercises the general-order tensor-SH dyadic advance with a
   !! genuinely non-axisymmetric, laterally-NON-uniform viscosity field.
   !!
   !! The cap height and the LVZ lateral extent use a smooth raised-cosine taper so
   !! the load and the viscosity are well resolved at the grid's lmax — the residual
   !! on/off-pole difference then measures the m>0 advance, not the cap's Gibbs ring.
   !!
   !! lmax defaults to 32 (in `make check`, a few seconds). Pass a larger lmax on the
   !! command line for the full-resolution standalone run, e.g.
   !!     make openmp=1 test_rotinv && bin/test_rotinv.x 192
   use fe_precision,       only: wp
   use fe_constants,       only: kyr
   use fe_earth_structure, only: earth_model, build_M3L70V01
   use fe_radial_fe,       only: radial_fe_finalize
   use fe_response,        only: response_destroy, response_commit_step, response_apply, response_begin_step, response_enable_lateral_visc, response, response_init_elastic, response_init_ve, response_init_null
   use fe_sht,             only: sht_grid, sht_grid_init, sht_grid_destroy, sht_grid_analysis, sht_grid_eval_point
   use fe_viscoelastic,    only: SCHEME_TRAP
   implicit none

   real(wp), parameter :: YR     = 0.001_wp*kyr
   real(wp), parameter :: DT     = 40.0_wp*YR        ! TRAP is A-stable: coarse step OK
   real(wp), parameter :: T_RAMP = 100.0_wp*YR, T_END = 200.0_wp*YR
   real(wp), parameter :: RHO_ICE = 931.0_wp, H_MAX = 200.0_wp
   real(wp), parameter :: CAP_IN = 0.20_wp, CAP_OUT = 0.40_wp  ! taper edges [rad] (~11.5°–23°)
   real(wp), parameter :: LVZ_TOP = 70.0e3_wp, LVZ_BOT = 170.0e3_wp
   real(wp), parameter :: LVZ_PERT = -2.0_wp         ! η 1e21 → 1e19 under the cap
   real(wp), parameter :: DEG = acos(-1.0_wp)/180.0_wp

   integer           :: lmax
   type(sht_grid)    :: sht
   type(earth_model) :: e
   real(wp)          :: u_pole, u_off, relerr, tol
   character(len=32) :: arg
   logical           :: ok

   lmax = 16
   if (command_argument_count() >= 1) then
      call get_command_argument(1, arg);  read(arg,*) lmax
   end if
   tol = max(0.01_wp, 48.0_wp/real(lmax,wp)/100.0_wp)   ! ~1.5% at lmax=32, 1% by lmax~48

   call sht_grid_init(sht, lmax, nlat=3*lmax, nphi=3*lmax, mmax=lmax)
   e = build_M3L70V01()

   write(*,'(a,i0,a,f5.2,a)') ' off-pole rotational-invariance (lmax=', lmax, &
        ', TRAP-3D, tol=', tol*100.0_wp, '%)'
   call run_case(0.0_wp,     0.0_wp,     u_pole)   ! axis = north pole (m=0)
   call run_case(55.0_wp*DEG, 40.0_wp*DEG, u_off)  ! off-pole axis (all m)

   relerr = abs(u_off - u_pole)/max(abs(u_pole), tiny(1.0_wp))
   write(*,'(a)') ''
   write(*,'(a,f10.4,a)') '   on-pole  cap-centre uplift = ', u_pole, ' m'
   write(*,'(a,f10.4,a)') '   off-pole cap-centre uplift = ', u_off,  ' m'
   write(*,'(a,es10.2)')  '   relative difference        = ', relerr
   ok = (relerr <= tol)

   call sht_grid_destroy(sht);  call radial_fe_finalize()
   write(*,'(a)') ''
   if (ok) then
      write(*,'(a)') ' PASS: off-pole LVZ matches on-pole peak uplift (rotational invariance)'
   else
      write(*,'(a)') ' FAIL: off-pole peak uplift differs from on-pole beyond tolerance'
      error stop 1
   end if

contains

   subroutine run_case(beta_c, lon_c, u_peak)
      !! Drive the cap+LVZ centred on axis n̂ = (colat beta_c, lon lon_c) to T_END and
      !! return the cap-centre uplift. Cap height and LVZ extent taper smoothly about n̂.
      real(wp), intent(in)  :: beta_c, lon_c   !! axis colatitude / longitude [rad]
      real(wp), intent(out) :: u_peak
      type(response) :: ve
      complex(wp), allocatable :: cap_lm(:), slm(:), ulm(:), nlm(:)
      real(wp),    allocatable :: cap(:,:), pert(:,:,:)
      real(wp) :: gam, rmid, depth, t, H, uval
      integer  :: i, j, ie, nstep, istep

      call response_init_ve(ve, e, sht, DT)
      ve%scheme = SCHEME_TRAP;  ve%max_couple_iter = 2
      allocate(cap_lm(sht%nlm), slm(sht%nlm), ulm(sht%nlm), nlm(sht%nlm))
      allocate(cap(sht%nphi, sht%nlat), pert(sht%nphi, sht%nlat, ve%ne))

      ! Unit cap shape on the grid (raised-cosine taper in angular distance from n̂).
      do j = 1, sht%nlat
         do i = 1, sht%nphi
            cap(i,j) = taper(ang_dist(sht%colat(j), sht%lon(i), beta_c, lon_c))
         end do
      end do
      call sht_grid_analysis(sht, cap, cap_lm)               ! NB: overwrites cap (unused after)

      ! LVZ: soft column under the cap, same lateral taper, in the depth band.
      pert = 0.0_wp
      do ie = 1, ve%ne
         rmid  = 0.5_wp*(ve%r(ie) + ve%r(ie+1))
         depth = e%r_earth - rmid
         if (depth > LVZ_TOP .and. depth < LVZ_BOT) then
            do j = 1, sht%nlat
               do i = 1, sht%nphi
                  pert(i,j,ie) = LVZ_PERT*taper(ang_dist(sht%colat(j), sht%lon(i), beta_c, lon_c))
               end do
            end do
         end if
      end do
      call response_enable_lateral_visc(ve, sht, pert)

      nstep = nint(T_END/DT)
      u_peak = 0.0_wp
      do istep = 1, nstep
         t   = real(istep,wp)*DT
         H   = min(t/T_RAMP, 1.0_wp)*H_MAX
         slm = RHO_ICE*H*cap_lm
         call response_begin_step(ve, sht)
         call response_apply(ve, sht, slm, ulm, nlm)
         call response_commit_step(ve, sht, slm)
         if (istep == nstep) then
            call sht_grid_eval_point(sht, ulm, beta_c, lon_c, uval)   ! uplift at the cap centre
            u_peak = uval
         end if
      end do

      deallocate(cap_lm, slm, ulm, nlm, cap, pert)
      call response_destroy(ve)
   end subroutine run_case

   pure real(wp) function ang_dist(th, ph, th0, ph0) result(g)
      !! Great-circle angle between grid point (th,ph) and axis (th0,ph0) [rad].
      real(wp), intent(in) :: th, ph, th0, ph0
      real(wp) :: c
      c = cos(th)*cos(th0) + sin(th)*sin(th0)*cos(ph - ph0)
      g = acos(max(-1.0_wp, min(1.0_wp, c)))
   end function ang_dist

   pure real(wp) function taper(g) result(s)
      !! Raised-cosine: 1 for g<=CAP_IN, 0 for g>=CAP_OUT, smooth between.
      real(wp), intent(in) :: g
      real(wp), parameter :: pi = acos(-1.0_wp)
      if (g <= CAP_IN) then
         s = 1.0_wp
      else if (g >= CAP_OUT) then
         s = 0.0_wp
      else
         s = 0.5_wp*(1.0_wp + cos(pi*(g - CAP_IN)/(CAP_OUT - CAP_IN)))
      end if
   end function taper

end program test_rotinv
