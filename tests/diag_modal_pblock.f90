program diag_modal_pblock
   !! DIAGNOSTIC (not a pass/fail test): does RESP_MODAL with n_modes=all converge
   !! to RESP_VE as the Arnoldi/Krylov block size n_krylov (p_block) grows?
   !!
   !! Motivation. In the modal-vs-VE ensemble the radial set — where modal(all)
   !! should reproduce VE exactly — plateaus at ~1.2 m RMSE, and n_modes=4/8/all
   !! all hit the SAME floor. The hypothesis: "all" is not all physical modes, it
   !! is "all modes recoverable from a p_block-dimensional load-Krylov subspace",
   !! hardcoded at 20. Truncated modes remove relaxation amplitude → the response
   !! is biased toward the (smaller) elastic limit → displacement is ATTENUATED,
   !! worst at low degree (continental scale, e.g. North America).
   !!
   !! This sweeps p_block over the FULL degree spectrum of the M3-L70-V01 radial
   !! benchmark and reports, per candidate:
   !!   (1) the fully-relaxed (fluid) sum-rule deficit  U_inf^VE(l) − U_inf^modal(l),
   !!       where U_inf^modal(l) = gu(l) + Σ_k Cu_k(l) and U_inf^VE(l) is RESP_VE
   !!       stepped to near-steady-state — the attenuation, quantified;
   !!   (2) the transient error vs RESP_VE driven by a unit white (m=0) load over
   !!       all degrees: the worst-case over time, and the value at the 26 ka
   !!       (present-day) mark of the deglaciation window.
   !! If the deficit collapses as p_block grows, the Krylov ceiling is confirmed
   !! and we know how high to set n_krylov.
   !!
   !!   usage:  diag_modal_pblock.x [lmax] [t_total_kyr] [dt_yr]
   !!   default:                     32     200           25
   use fe_precision,       only: wp
   use fe_constants,       only: kyr
   use fe_earth_structure, only: earth_model, build_M3L70V01
   use fe_radial_fe,       only: radial_fe_finalize
   use fe_response,        only: response, response_init_modal, response_init_ve, &
                                 response_begin_step, response_apply, &
                                 response_commit_step, response_set_dt, response_destroy
   use fe_sht,             only: sht_grid, sht_grid_init, sht_grid_destroy, sht_grid_lmidx
   use fe_viscoelastic,    only: NLAM
   implicit none

   integer, parameter :: NP = 5
   integer            :: pblocks(NP) = [10, 20, 40, 60, 80]
   integer            :: npk
   integer            :: lmax, nstep, ip, l, m, lm, k, base, istep, ntab
   real(wp)           :: dt, t_total, t_pd, t, swing
   character(len=64)  :: arg
   type(sht_grid)     :: sht
   type(earth_model)  :: e
   type(response)     :: ve
   type(response)     :: md(NP)
   complex(wp), allocatable :: slm(:), ulm(:), nlm(:)
   real(wp),    allocatable :: uve(:)              ! (0:lmax) VE u(l,0) this step
   real(wp),    allocatable :: uel(:)              ! (0:lmax) elastic gain gu(l)
   real(wp),    allocatable :: uinf_ve(:)          ! (0:lmax) VE stepped fluid gain (final)
   real(wp),    allocatable :: uinf_md(:,:)        ! (0:lmax,NP) modal analytic fluid sum gu+ΣCu
   real(wp),    allocatable :: def_step(:,:)       ! (0:lmax,NP) stepped deficit uve-umd at T
   ! geoid-N channel — the quantity rsl=N-u depends on but the u-only checks never tested
   real(wp),    allocatable :: nve(:), nel(:), ninf_ve(:), ndef_step(:,:)
   real(wp),    allocatable :: maxtr(:), pdrms(:), pdmax(:)
   integer,     allocatable :: modes(:)
   real(wp)                 :: umd, nmd, d, sse, cnt
   integer                  :: ipd, ip20

   ! ---- args ---------------------------------------------------------------
   lmax = 32;  t_total = 200.0_wp;  dt = 25.0_wp
   if (command_argument_count() >= 1) then; call get_command_argument(1, arg); read(arg,*) lmax;     end if
   if (command_argument_count() >= 2) then; call get_command_argument(2, arg); read(arg,*) t_total;  end if
   if (command_argument_count() >= 3) then; call get_command_argument(3, arg); read(arg,*) dt;       end if
   dt      = dt*1.0e-3_wp*kyr          ! yr -> s  (kyr is seconds-per-kyr)
   t_total = t_total*kyr
   t_pd    = 26.0_wp*kyr               ! deglaciation window length (LGM->PD)
   nstep   = max(1, nint(t_total/dt))
   ipd     = max(1, nint(t_pd/dt))     ! step index nearest the PD mark

   ! ---- setup --------------------------------------------------------------
   call sht_grid_init(sht, lmax, nlat=2*lmax, nphi=4*lmax)
   e = build_M3L70V01()
   allocate(slm(sht%nlm), ulm(sht%nlm), nlm(sht%nlm))
   allocate(uve(0:lmax), uel(0:lmax), uinf_ve(0:lmax))
   allocate(uinf_md(0:lmax,NP), def_step(0:lmax,NP))
   allocate(nve(0:lmax), nel(0:lmax), ninf_ve(0:lmax), ndef_step(0:lmax,NP))
   allocate(maxtr(NP), pdrms(NP), pdmax(NP), modes(NP))
   maxtr = 0.0_wp;  pdrms = 0.0_wp;  pdmax = 0.0_wp

   ! unit WHITE load: one unit coefficient at every degree (m=0), l=1..lmax
   slm = (0.0_wp,0.0_wp)
   do l = 1, lmax
      slm(sht_grid_lmidx(sht, l, 0)) = (1.0_wp,0.0_wp)
   end do

   ! RESP_VE reference (FE scheme) + the modal candidates, all on the same dt
   call response_init_ve(ve, e, sht, dt)
   call response_set_dt(ve, dt)
   do ip = 1, NP
      call response_init_modal(md(ip), e, sht, n_modes=-1, mode_rank=1, &
                               dt_be=5.0_wp*kyr, p_block=pblocks(ip))
      call response_set_dt(md(ip), dt)
      modes(ip) = sum(md(ip)%nmode_deg(1:lmax))
      do l = 1, lmax
         base = md(ip)%spec_off(l)
         uinf_md(l,ip) = md(ip)%gu(l) + sum(md(ip)%mCu(base+1:base+md(ip)%nmode_deg(l)))
      end do
   end do
   uel(0:lmax) = md(1)%gu(0:lmax)     ! elastic gain (same for all)
   nel(0:lmax) = md(1)%gn(0:lmax)     ! elastic geoid gain
   npk = 3*NLAM*count(md(1)%lat_mw)   ! per-degree memory dimension = hard mode ceiling
   ip20 = 2                           ! the default-block (p=20) candidate, for the u-vs-N table

   write(*,'(a)')        ' === RESP_MODAL p_block (Krylov) sweep vs RESP_VE — M3-L70-V01, radial ==='
   write(*,'(a,i0,a,f6.1,a,f6.1,a,i0,a)') '   lmax=', lmax, '   T=', t_total/kyr, ' ka   dt=', &
        dt/kyr*1.0e3_wp, ' yr   (', nstep, ' steps, white m=0 load over all degrees)'
   write(*,'(a,i0,a)') '   per-degree memory dimension npk = ', npk, &
        ' (the hard mode ceiling; p_block>~npk reintroduces spurious modes)'
   write(*,'(a)') ''

   ! ---- lockstep time integration -----------------------------------------
   do istep = 1, nstep
      t = real(istep-1, wp)*dt
      ! RESP_VE
      call response_begin_step(ve, sht)
      call response_apply(ve, sht, slm, ulm, nlm)
      do l = 1, lmax
         uve(l) = real(ulm(sht_grid_lmidx(sht,l,0)), wp)
         nve(l) = real(nlm(sht_grid_lmidx(sht,l,0)), wp)
      end do
      call response_commit_step(ve, sht, slm)
      ! each modal candidate vs VE
      do ip = 1, NP
         call response_begin_step(md(ip), sht)
         call response_apply(md(ip), sht, slm, ulm, nlm)
         sse = 0.0_wp;  cnt = 0.0_wp
         do l = 1, lmax
            umd = real(ulm(sht_grid_lmidx(sht,l,0)), wp)
            nmd = real(nlm(sht_grid_lmidx(sht,l,0)), wp)
            d   = umd - uve(l)
            maxtr(ip) = max(maxtr(ip), abs(d))
            if (istep == ipd) then
               sse = sse + d*d;  cnt = cnt + 1.0_wp
               pdmax(ip) = max(pdmax(ip), abs(d))
            end if
            if (istep == nstep) then
               def_step(l,ip)  = uve(l) - umd
               ndef_step(l,ip) = nve(l) - nmd
            end if
         end do
         if (istep == ipd) pdrms(ip) = sqrt(sse/max(cnt,1.0_wp))
         call response_commit_step(md(ip), sht, slm)
      end do
      if (istep == nstep) then
         do l = 1, lmax;  uinf_ve(l) = uve(l);  ninf_ve(l) = nve(l);  end do
      end if
   end do

   ! ---- report -------------------------------------------------------------
   write(*,'(a)') ' Summary (per p_block):'
   write(*,'(a)') '   p_block  Σmodes   transient max|Δu|   Δu@26ka rms    Δu@26ka max   fluid deficit max'
   write(*,'(a)') '                        [m]                [m]            [m]           [m] (@ degree)'
   do ip = 1, NP
      block
         real(wp) :: dmax;  integer :: lworst
         dmax = 0.0_wp;  lworst = 0
         do l = 1, lmax
            if (abs(def_step(l,ip)) > dmax) then;  dmax = abs(def_step(l,ip));  lworst = l;  end if
         end do
         write(*,'(i8,i9,es18.3,es15.3,es15.3,es14.3,a,i0,a)') &
              pblocks(ip), modes(ip), maxtr(ip), pdrms(ip), pdmax(ip), dmax, '  (l=', lworst, ')'
      end block
   end do
   write(*,'(a)') ''

   ! per-degree fluid deficit (VE stepped to T  −  modal stepped to T), low degrees
   ntab = min(12, lmax)
   write(*,'(a,f5.0,a)') ' Per-degree fluid-limit deficit  U_inf^VE − U_inf^modal  at T=', t_total/kyr, &
        ' ka  [m]  (load-relevant, low degrees):'
   write(*,'(a)', advance='no') '   l   U_inf^VE   swing'
   do ip = 1, NP;  write(*,'(a,i0)', advance='no') '    p=', pblocks(ip);  end do
   write(*,'(a)') ''
   do l = 1, ntab
      swing = uinf_ve(l) - uel(l)
      write(*,'(i4,es11.3,es9.2)', advance='no') l, uinf_ve(l), swing
      do ip = 1, NP;  write(*,'(es10.2)', advance='no') def_step(l,ip);  end do
      write(*,'(a)') ''
   end do
   write(*,'(a)') ''
   write(*,'(a)') ' Note: positive deficit = modal under-relaxes (attenuated vs VE). If the deficit'
   write(*,'(a)') '       shrinks left-to-right, the Krylov block size (n_krylov), not mode ranking,'
   write(*,'(a)') '       is the accuracy ceiling for n_modes=all.'
   write(*,'(a)') ''

   ! ---- the decisive check: geoid N deficit alongside displacement u ----------
   ! rsl = N - u. The u channel matches VE to ~1e-8 m; if the N channel does NOT,
   ! the SLE-coupled rsl gap is a GEOID-response error, not a coupling error.
   write(*,'(a,i0,a)') ' Geoid N vs displacement u — fluid-limit deficit at p_block=', pblocks(ip20), &
        ' (VE − modal) [m], and as a fraction of each channel''s own swing:'
   write(*,'(a)') '   l      u_swing   u_def    u_def/sw      N_swing   N_def    N_def/sw'
   do l = 1, ntab
      block
         real(wp) :: usw, nsw
         usw = uinf_ve(l) - uel(l);  nsw = ninf_ve(l) - nel(l)
         write(*,'(i4,es12.3,es10.2,es12.2,a,es12.3,es10.2,es12.2)') l, &
              usw, def_step(l,ip20),  def_step(l,ip20)/sign(max(abs(usw),tiny(1.0_wp)),usw), '   ', &
              nsw, ndef_step(l,ip20), ndef_step(l,ip20)/sign(max(abs(nsw),tiny(1.0_wp)),nsw)
      end block
   end do
   write(*,'(a)') ''
   write(*,'(a)') ' Read: if N_def/sw ≫ u_def/sw, the modal GEOID response is the weak link'
   write(*,'(a)') ' (fix in mode_residue rn / elastic_gains gn); if both are ~1e-4, the rsl gap'
   write(*,'(a)') ' is the SLE memory<->load coupling, not the per-degree N response.'

   do ip = 1, NP;  call response_destroy(md(ip));  end do
   call response_destroy(ve)
   call sht_grid_destroy(sht);  call radial_fe_finalize()

end program diag_modal_pblock
