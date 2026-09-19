program test_modal_resp
   !! RESP_MODAL wired into the response interface (fe_response).
   !!
   !! (1) Plumbing: driving a held single-degree load through
   !!     begin_step/apply/commit_step must reproduce the analytic step response
   !!     of the response's OWN stored modal spectrum, gu + Σ_i Cu_i(1−e^{−t/τ_i}),
   !!     to machine precision — confirms the φ recurrence + ragged storage.
   !! (2) End-to-end: the same load driven through RESP_MODAL vs RESP_VE (the full
   !!     FE model) must agree to the modal-approximation tolerance, confirming the
   !!     reduced response is a faithful drop-in for the SLE/coupling interface.
   use fe_precision,       only: wp
   use fe_constants,       only: kyr
   use fe_earth_structure, only: earth_model, build_M3L70V01
   use fe_radial_fe,       only: radial_fe_finalize
   use fe_response,        only: response, response_init_modal, response_init_ve, &
                                 response_begin_step, response_apply, &
                                 response_commit_step, response_set_dt, response_destroy
   use fe_sht,             only: sht_grid, sht_grid_init, sht_grid_destroy, sht_grid_lmidx
   implicit none

   integer, parameter :: LMAX = 8, J = 2, NSTEP = 400
   type(sht_grid)    :: sht
   type(earth_model) :: e
   type(response)    :: md, ve
   complex(wp), allocatable :: slm(:), ulm(:), nlm(:)
   real(wp) :: dt, t, umod, uana, uve, swing, plumb_err, e2e_err, u0, uinf
   integer  :: lm, i, k, base
   logical  :: ok

   ok = .true.
   dt = 0.025_wp*kyr                          ! 25 yr (stable for RESP_VE FE)
   call sht_grid_init(sht, LMAX, nlat=2*LMAX, nphi=4*LMAX)
   e = build_M3L70V01()
   lm = sht_grid_lmidx(sht, J, 0)
   allocate(slm(sht%nlm), ulm(sht%nlm), nlm(sht%nlm))
   slm = (0.0_wp,0.0_wp);  slm(lm) = (1.0_wp,0.0_wp)

   call response_init_modal(md, e, sht, n_modes=-1, mode_rank=1, dt_be=5.0_wp*kyr)
   call response_set_dt(md, dt)
   call response_init_ve(ve, e, sht, dt)      ! FE scheme (default) — full reference

   base = md%spec_off(J)
   plumb_err = 0.0_wp;  e2e_err = 0.0_wp;  swing = 0.0_wp
   u0 = md%gu(J);  uinf = md%gu(J) + sum(md%mCu(base+1:base+md%nmode_deg(J)))
   swing = abs(uinf - u0)

   do i = 1, NSTEP
      t = real(i-1, wp)*dt
      ! analytic modal step response from the stored spectrum (held σ=1 from t=0)
      uana = md%gu(J)
      do k = 1, md%nmode_deg(J)
         uana = uana + md%mCu(base+k)*(1.0_wp - exp(-t/md%mtau(base+k)))
      end do
      ! RESP_MODAL through the interface
      call response_begin_step(md, sht)
      call response_apply(md, sht, slm, ulm, nlm)
      umod = real(ulm(lm), wp)
      call response_commit_step(md, sht, slm)
      ! RESP_VE reference through the interface
      call response_begin_step(ve, sht)
      call response_apply(ve, sht, slm, ulm, nlm)
      uve = real(ulm(lm), wp)
      call response_commit_step(ve, sht, slm)

      plumb_err = max(plumb_err, abs(umod - uana))
      e2e_err   = max(e2e_err,   abs(umod - uve))
   end do

   write(*,'(a)')        ' RESP_MODAL via the response interface (M3-L70-V01, degree 2)'
   write(*,'(a,i0)')     '   modes kept at degree 2          = ', md%nmode_deg(J)
   write(*,'(a,es10.3)') '   elastic U0, fluid Uinf swing    = ', swing
   write(*,'(a,es10.3)') '   plumbing err |modal - analytic| = ', plumb_err/swing
   write(*,'(a,es10.3)') '   end-to-end |modal - RESP_VE|    = ', e2e_err/swing

   if (plumb_err/swing > 1.0e-10_wp) then
      write(*,'(a)') '   FAIL: interface stepping /= analytic modal step response';  ok = .false.
   end if
   if (e2e_err/swing > 3.0e-2_wp) then
      write(*,'(a)') '   FAIL: RESP_MODAL disagrees with RESP_VE beyond approx tol';  ok = .false.
   end if

   call response_destroy(md);  call response_destroy(ve)
   call sht_grid_destroy(sht);  call radial_fe_finalize()
   write(*,'(a)') ''
   if (ok) then
      write(*,'(a)') ' PASS: RESP_MODAL interface = analytic spectrum, tracks RESP_VE'
   else
      write(*,'(a)') ' FAIL: RESP_MODAL interface checks failed'
      error stop 1
   end if

end program test_modal_resp
