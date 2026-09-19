program diag_modal_ramp
   !! DIAGNOSTIC (not pass/fail): is the modal-vs-VE discrepancy TEMPORAL (the
   !! one-exact-step-per-couple ZOH load treatment) rather than spectral (mode
   !! count)? The companion diag_modal_pblock showed the per-degree modal operator
   !! matches RESP_VE to ~1e-7 m for a HELD load at any n_krylov — so the ensemble
   !! gap must come from how a TIME-VARYING load is integrated.
   !!
   !! In the ensemble, modal takes 1 exact-exponential step per 100-yr couple with
   !! the load held at the step endpoint (ZOH), while RESP_VE sub-steps (n_solve
   !! 1.6 radial → 8.6 deglac3d). This probe drives a deglaciation-like UNLOADING
   !! ramp (hold loaded, ramp to zero over 26 ka, then rebound) through four
   !! integrations of the SAME load history and compares them to a temporally
   !! converged truth (VE, many sub-steps):
   !!
   !!   (1) modal,  1 step/couple   (ZOH)        — what the model does now
   !!   (2) VE,     1 step/couple   (ZOH)        — single-step reference
   !!   (3) modal,  N substeps/couple            — proposed-fix proxy (sub-stepping)
   !!   (4) VE,     N substeps/couple   == TRUTH
   !!
   !! Hypothesis confirmed if  err(modal-1) ≈ err(VE-1) ≫ err(modal-N) ≈ 0:
   !! i.e. modal's error is the SAME temporal ZOH error a single-step VE makes,
   !! and sub-stepping (or a linear-σ modal advance) closes it. If instead
   !! err(modal-1) ≫ err(VE-1), the error is extra modal/spectral, not temporal.
   !!
   !!   usage:  diag_modal_ramp.x [lmax] [nsub] [dt_couple_yr]
   !!   default:                   24    8      100
   use fe_precision,       only: wp
   use fe_constants,       only: kyr
   use fe_earth_structure, only: earth_model, build_M3L70V01
   use fe_radial_fe,       only: radial_fe_finalize
   use fe_response,        only: response, response_init_modal, response_init_ve, &
                                 response_begin_step, response_apply, &
                                 response_commit_step, response_set_dt, response_destroy
   use fe_sht,             only: sht_grid, sht_grid_init, sht_grid_destroy, sht_grid_lmidx
   implicit none

   integer            :: lmax, nsub, ncpl, n, l, ireport, ip
   real(wp)           :: dtc, t_hold, t_ramp, t_reb, t_total, tnp1, signp1
   character(len=64)  :: arg
   type(sht_grid)     :: sht
   type(earth_model)  :: e
   type(response)     :: mz, vz, ms, vf       ! modal-ZOH, VE-ZOH, modal-sub, VE-fine(truth)
   complex(wp), allocatable :: slm(:), ulm(:), nlm(:)
   real(wp),    allocatable :: u_mz(:), u_vz(:), u_ms(:), u_vf(:)   ! (0:lmax) at report time
   ! running max over time of |traj - truth|, per trajectory (aggregate over degrees)
   real(wp) :: mx_mz, mx_vz, mx_ms
   ! per-degree error snapshots at end-of-ramp (PD) — the regime the ensemble scores
   real(wp), allocatable :: epd_mz(:), epd_vz(:), epd_ms(:), upd_vf(:)
   integer  :: ntab

   ! ---- args ---------------------------------------------------------------
   lmax = 24;  nsub = 8;  dtc = 100.0_wp
   if (command_argument_count() >= 1) then; call get_command_argument(1, arg); read(arg,*) lmax; end if
   if (command_argument_count() >= 2) then; call get_command_argument(2, arg); read(arg,*) nsub; end if
   if (command_argument_count() >= 3) then; call get_command_argument(3, arg); read(arg,*) dtc;  end if
   dtc    = dtc*1.0e-3_wp*kyr          ! yr -> s
   t_hold = 40.0_wp*kyr               ! pre-ramp loaded hold (build the relaxed state)
   t_ramp = 26.0_wp*kyr               ! deglaciation window (load 1 -> 0)
   t_reb  = 40.0_wp*kyr               ! post-glacial rebound (load 0)
   t_total = t_hold + t_ramp + t_reb
   ncpl    = nint(t_total/dtc)
   ireport = nint((t_hold + t_ramp)/dtc)   ! couple index at end-of-ramp (PD)

   ! ---- setup --------------------------------------------------------------
   call sht_grid_init(sht, lmax, nlat=2*lmax, nphi=4*lmax)
   e = build_M3L70V01()
   allocate(slm(sht%nlm), ulm(sht%nlm), nlm(sht%nlm))
   allocate(u_mz(0:lmax), u_vz(0:lmax), u_ms(0:lmax), u_vf(0:lmax))
   allocate(epd_mz(0:lmax), epd_vz(0:lmax), epd_ms(0:lmax), upd_vf(0:lmax))

   call response_init_modal(mz, e, sht, n_modes=-1, mode_rank=1, dt_be=5.0_wp*kyr)
   call response_init_modal(ms, e, sht, n_modes=-1, mode_rank=1, dt_be=5.0_wp*kyr)
   call response_init_ve(vz, e, sht, dtc)
   call response_init_ve(vf, e, sht, dtc)

   mx_mz = 0.0_wp;  mx_vz = 0.0_wp;  mx_ms = 0.0_wp

   write(*,'(a)') ' === modal vs VE: TEMPORAL (sub-stepping) probe — M3-L70-V01, deglaciation ramp ==='
   write(*,'(a,i0,a,i0,a,f6.1,a)') '   lmax=', lmax, '   nsub=', nsub, '   dt_couple=', dtc/kyr*1.0e3_wp, ' yr'
   write(*,'(a,i0,a)') '   load: hold 40 ka -> ramp 1->0 over 26 ka -> rebound 40 ka   (', ncpl, ' couples)'
   write(*,'(a)') ''

   ! ---- couple loop --------------------------------------------------------
   do n = 1, ncpl
      tnp1   = real(n,wp)*dtc
      signp1 = ramp(tnp1)

      call advance_zoh(mz, dtc, signp1)
      call advance_zoh(vz, dtc, signp1)
      call advance_sub(ms, real(n-1,wp)*dtc, dtc, nsub)
      call advance_sub(vf, real(n-1,wp)*dtc, dtc, nsub)

      call read_u(mz, signp1, u_mz)
      call read_u(vz, signp1, u_vz)
      call read_u(ms, signp1, u_ms)
      call read_u(vf, signp1, u_vf)      ! truth

      do l = 1, lmax
         mx_mz = max(mx_mz, abs(u_mz(l) - u_vf(l)))
         mx_vz = max(mx_vz, abs(u_vz(l) - u_vf(l)))
         mx_ms = max(mx_ms, abs(u_ms(l) - u_vf(l)))
      end do
      if (n == ireport) then
         do l = 0, lmax
            epd_mz(l) = u_mz(l) - u_vf(l);  epd_vz(l) = u_vz(l) - u_vf(l)
            epd_ms(l) = u_ms(l) - u_vf(l);  upd_vf(l) = u_vf(l)
         end do
      end if
   end do

   ! ---- report -------------------------------------------------------------
   write(*,'(a)') ' Max over time of |trajectory - truth(VE,nsub)|, aggregate over degrees [m]:'
   write(*,'(a,es11.3,a)') '   modal, 1 step/couple  (ZOH, current model)   = ', mx_mz, '   <- the modal "error"'
   write(*,'(a,es11.3,a)') '   VE,    1 step/couple  (ZOH)                  = ', mx_vz, '   <- same temporal error?'
   write(*,'(a,es11.3,a)') '   modal, N substeps/couple (proposed fix)      = ', mx_ms, '   <- closes the gap?'
   write(*,'(a)') ''
   ntab = min(12, lmax)
   write(*,'(a)') ' Per-degree error at end-of-ramp (present day) [m]:'
   write(*,'(a)') '   l    u_truth     modal-1     VE-1       modal-N'
   do l = 1, ntab
      write(*,'(i4,es12.3,es12.3,es11.3,es12.3)') l, upd_vf(l), epd_mz(l), epd_vz(l), epd_ms(l)
   end do
   write(*,'(a)') ''
   write(*,'(a)') ' Read: if modal-1 ≈ VE-1 (same column-2/3 magnitudes) and modal-N ≈ 0, the modal'
   write(*,'(a)') ' error is the single-step TEMPORAL (ZOH) error, removed by sub-stepping / linear-σ.'

   call response_destroy(mz);  call response_destroy(vz)
   call response_destroy(ms);  call response_destroy(vf)
   call sht_grid_destroy(sht);  call radial_fe_finalize()

contains

   real(wp) function ramp(t) result(f)
      !! Deglaciation-like load factor: 1 while loaded, linear 1->0 over the ramp
      !! window, 0 afterward.
      real(wp), intent(in) :: t
      if (t <= t_hold) then
         f = 1.0_wp
      else if (t <= t_hold + t_ramp) then
         f = 1.0_wp - (t - t_hold)/t_ramp
      else
         f = 0.0_wp
      end if
   end function ramp

   subroutine set_white(sigval)
      !! White m=0 load: every degree l=1..lmax carries coefficient sigval.
      real(wp), intent(in) :: sigval
      slm = (0.0_wp,0.0_wp)
      do l = 1, lmax
         slm(sht_grid_lmidx(sht,l,0)) = cmplx(sigval, 0.0_wp, wp)
      end do
   end subroutine set_white

   subroutine advance_zoh(r, dt, sigval)
      !! One step of size dt with the load held at sigval (endpoint ZOH).
      type(response), intent(inout) :: r
      real(wp),       intent(in)    :: dt, sigval
      call response_set_dt(r, dt)
      call set_white(sigval)
      call response_begin_step(r, sht)
      call response_apply(r, sht, slm, ulm, nlm)
      call response_commit_step(r, sht, slm)
   end subroutine advance_zoh

   subroutine advance_sub(r, t_n, dt, ns)
      !! Advance the couple [t_n, t_n+dt] in ns equal substeps, sampling the ramp
      !! at each substep endpoint (a fine-ZOH ≈ exact-ramp integration).
      type(response), intent(inout) :: r
      real(wp),       intent(in)    :: t_n, dt
      integer,        intent(in)    :: ns
      integer  :: s
      real(wp) :: dts, ts
      dts = dt/real(ns,wp)
      do s = 1, ns
         ts = t_n + real(s,wp)*dts
         call advance_zoh(r, dts, ramp(ts))
      end do
   end subroutine advance_sub

   subroutine read_u(r, sigval, u)
      !! Physical displacement at the current memory time: elastic(sigval) +
      !! drift(memory). begin+apply only (no state change).
      type(response), intent(inout) :: r
      real(wp),       intent(in)    :: sigval
      real(wp),       intent(out)   :: u(0:)
      call set_white(sigval)
      call response_begin_step(r, sht)
      call response_apply(r, sht, slm, ulm, nlm)
      do l = 0, lmax
         u(l) = real(ulm(sht_grid_lmidx(sht,l,0)), wp)
      end do
   end subroutine read_u

end program diag_modal_ramp
