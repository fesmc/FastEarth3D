program test_modal_visc3d
   !! RESP_MODAL lateral viscosity (design-modal.md §4): depth-weighted split-operator
   !! rate modulation, validated by the uniform-η invariants (mirrors test_response_3d):
   !!
   !!   (A) uniform η perturbation = 0, the real-space anomaly path FORCED on
   !!       (visc3d_tol<0): the per-rank SHT round-trip must be the identity, so the
   !!       lateral response reproduces the 1-D modal advance to SHT precision.
   !!   (B) uniform η scaled by 10^p: the mean modulation makes the rate 10^(−p)/τ, i.e.
   !!       the 1-D advance at a step Δt·10^(−p). Comparing against a 1-D modal driven at
   !!       that scaled Δt isolates the mean factor + SHT identity (same eigensolve, so no
   !!       eigensolve-reproducibility noise), and must match to SHT precision.
   !!   (C) a genuinely non-uniform field must RUN and actually change the response vs the
   !!       1-D advance (the anomaly substep is wired and active).
   use fe_precision,       only: wp
   use fe_constants,       only: kyr, pi
   use fe_earth_structure, only: earth_model, build_M3L70V01
   use fe_radial_fe,       only: radial_fe_finalize
   use fe_response,        only: response, response_init_modal, response_set_dt, &
                                 response_begin_step, response_apply, response_commit_step, &
                                 response_destroy, response_enable_lateral_visc_modal, &
                                 LAT_LIE_CHAR, LAT_COUPLED
   use fe_sht,             only: sht_grid, sht_grid_init, sht_grid_destroy, sht_grid_lmidx
   implicit none

   integer, parameter :: LMAX = 8, NSTEP = 200
   type(sht_grid)    :: sht
   type(earth_model) :: e
   type(response)    :: refA, latA, refB, latB, refC, latC, refD, latD
   complex(wp), allocatable :: slm(:)
   real(wp), allocatable    :: pert(:,:,:)
   real(wp) :: dt, p, dzero, dscaled, dnon, dcoup
   logical  :: ok

   ok = .true.
   dt = 0.05_wp*kyr
   p  = -0.4_wp
   call sht_grid_init(sht, LMAX, nlat=2*LMAX, nphi=4*LMAX)
   e = build_M3L70V01()
   allocate(slm(sht%nlm));  call make_load(sht, slm)

   ! (A) uniform η perturbation 0, anomaly path forced -> SHT round-trip identity
   call response_init_modal(refA, e, sht, n_modes=-1, mode_rank=1, dt_be=5.0_wp*kyr)
   call response_set_dt(refA, dt)
   call response_init_modal(latA, e, sht, n_modes=-1, mode_rank=1, dt_be=5.0_wp*kyr)
   call response_set_dt(latA, dt)
   allocate(pert(sht%nphi, sht%nlat, latA%ne));  pert = 0.0_wp
   latA%lat_method = LAT_LIE_CHAR                               ! cases A–C exercise the LIE split path
   latA%visc3d_tol = -1.0_wp                                    ! force the anomaly SHT path
   call response_enable_lateral_visc_modal(latA, sht, pert)
   call drive_compare(refA, latA, dzero)

   ! (B) uniform η×10^p: mean rate 10^(−p)/τ == 1-D advance at Δt·10^(−p) (same eigensolve)
   call response_init_modal(refB, e, sht, n_modes=-1, mode_rank=1, dt_be=5.0_wp*kyr)
   call response_set_dt(refB, dt*10.0_wp**(-p))
   call response_init_modal(latB, e, sht, n_modes=-1, mode_rank=1, dt_be=5.0_wp*kyr)
   call response_set_dt(latB, dt)
   pert = p
   latB%lat_method = LAT_LIE_CHAR
   latB%visc3d_tol = -1.0_wp
   call response_enable_lateral_visc_modal(latB, sht, pert)
   call drive_compare(refB, latB, dscaled)

   ! (C) non-uniform field: must run and change the response vs 1-D
   call response_init_modal(refC, e, sht, n_modes=-1, mode_rank=1, dt_be=5.0_wp*kyr)
   call response_set_dt(refC, dt)
   call response_init_modal(latC, e, sht, n_modes=-1, mode_rank=1, dt_be=5.0_wp*kyr)
   call response_set_dt(latC, dt)
   call make_pert_nonuniform(sht, latC%ne, pert)
   latC%lat_method = LAT_LIE_CHAR
   call response_enable_lateral_visc_modal(latC, sht, pert)     ! default tol -> selects active ranks
   call drive_compare(refC, latC, dnon)

   ! (D) COUPLED path, uniform η×10^p forced active: the per-rank coupled rate operator
   !     exp(−Δt·L_i) (Arnoldi) must reduce to the 1-D advance at Δt·10^(−p), to ~SHT
   !     precision — exercises modal_lateral_coupled / expm_small under a known answer.
   call response_init_modal(refD, e, sht, n_modes=-1, mode_rank=1, dt_be=5.0_wp*kyr)
   call response_set_dt(refD, dt*10.0_wp**(-p))
   call response_init_modal(latD, e, sht, n_modes=-1, mode_rank=1, dt_be=5.0_wp*kyr)
   call response_set_dt(latD, dt)
   latD%lat_method = LAT_COUPLED
   latD%visc3d_tol = -1.0_wp                                    ! force all ranks active (coupled)
   pert = p
   call response_enable_lateral_visc_modal(latD, sht, pert)
   call drive_compare(refD, latD, dcoup)

   write(*,'(a)')        ' RESP_MODAL lateral viscosity (split-operator, M3-L70-V01)'
   write(*,'(a,i0)')     '   active anomaly ranks (case C)   = ', latC%nrank3d
   write(*,'(a,es10.3)') '   (A) uniform η=0,  |lat - 1D|    = ', dzero
   write(*,'(a,es10.3)') '   (B) uniform η×10^p vs Δt-scaled = ', dscaled
   write(*,'(a,es10.3)') '   (C) non-uniform, |lat - 1D|     = ', dnon
   write(*,'(a,es10.3)') '   (D) COUPLED uniform vs Δt-scaled= ', dcoup

   if (dzero   > 1.0e-9_wp) then; write(*,'(a)') '   FAIL: uniform-zero lateral /= 1-D (SHT identity)'; ok=.false.; end if
   if (dscaled > 1.0e-9_wp) then; write(*,'(a)') '   FAIL: uniform-scaled lateral /= Δt-scaled 1-D';   ok=.false.; end if
   if (latC%nrank3d < 1)    then; write(*,'(a)') '   FAIL: non-uniform field activated no anomaly ranks'; ok=.false.; end if
   if (dnon    < 1.0e-3_wp) then; write(*,'(a)') '   FAIL: non-uniform lateral has no effect on the response'; ok=.false.; end if
   if (dcoup   > 1.0e-7_wp) then; write(*,'(a)') '   FAIL: COUPLED uniform /= Δt-scaled 1-D (operator/expm)'; ok=.false.; end if

   call response_destroy(refA);  call response_destroy(latA)
   call response_destroy(refB);  call response_destroy(latB)
   call response_destroy(refC);  call response_destroy(latC)
   call response_destroy(refD);  call response_destroy(latD)
   call sht_grid_destroy(sht);    call radial_fe_finalize()

   write(*,'(a)') ''
   if (ok) then
      write(*,'(a)') ' PASS: modal lateral viscosity reduces to 1-D for uniform η and is active otherwise'
   else
      write(*,'(a)') ' FAIL: modal lateral-viscosity checks failed'
      error stop 1
   end if

contains

   subroutine make_load(sht, slm)
      !! A held multi-(l,m) load (incl. m>0) to exercise the general-order path.
      type(sht_grid), intent(in)  :: sht
      complex(wp),    intent(out) :: slm(:)
      slm = (0.0_wp, 0.0_wp)
      slm(sht_grid_lmidx(sht, 2, 0)) = (1.0_wp,  0.0_wp)
      slm(sht_grid_lmidx(sht, 2, 1)) = (0.5_wp,  0.3_wp)
      slm(sht_grid_lmidx(sht, 4, 2)) = (0.4_wp, -0.2_wp)
      slm(sht_grid_lmidx(sht, 6, 0)) = (0.3_wp,  0.0_wp)
   end subroutine make_load

   subroutine make_pert_nonuniform(sht, ne, pert)
      !! Smooth ±0.5-dex lateral log10(η) field with a mild depth gradient, so different
      !! depth bands (hence different modes) see different lateral structure.
      type(sht_grid), intent(in)  :: sht
      integer,        intent(in)  :: ne
      real(wp),       intent(out) :: pert(:,:,:)
      integer  :: i, j, ee
      real(wp) :: lat_pat, depth_fac
      do ee = 1, ne
         depth_fac = 0.5_wp + 0.5_wp*real(ee-1, wp)/real(max(ne-1,1), wp)   ! 0.5 -> 1.0 with depth
         do j = 1, sht%nlat
            do i = 1, sht%nphi
               lat_pat = cos(2.0_wp*sht%lon(i))*sin(sht%colat(j))
               pert(i,j,ee) = 0.5_wp*depth_fac*lat_pat
            end do
         end do
      end do
   end subroutine make_pert_nonuniform

   subroutine drive_compare(ra, rb, drel)
      !! Drive both responses NSTEP held-load steps and return the worst uplift-field
      !! difference relative to the reference swing.
      type(response), intent(inout) :: ra, rb
      real(wp),       intent(out)   :: drel
      complex(wp), allocatable :: ua(:), na(:), ub(:), nb(:)
      real(wp) :: dmax, umax
      integer  :: i
      allocate(ua(sht%nlm), na(sht%nlm), ub(sht%nlm), nb(sht%nlm))
      dmax = 0.0_wp;  umax = 0.0_wp
      do i = 1, NSTEP
         call response_begin_step(ra, sht);  call response_apply(ra, sht, slm, ua, na)
         call response_commit_step(ra, sht, slm)
         call response_begin_step(rb, sht);  call response_apply(rb, sht, slm, ub, nb)
         call response_commit_step(rb, sht, slm)
         dmax = max(dmax, maxval(abs(ua - ub)))
         umax = max(umax, maxval(abs(ua)))
      end do
      drel = dmax/max(umax, tiny(1.0_wp))
   end subroutine drive_compare

end program test_modal_visc3d
