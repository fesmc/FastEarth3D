program diag_modal_sle
   !! DIAGNOSTIC (not pass/fail): localise the modal-vs-VE rsl gap. Earlier probes
   !! showed the per-degree response (both u AND geoid N) matches VE to ~1e-4 of
   !! swing for a HELD load, so the gap is not the response operator. This drives
   !! the FULL SLE (fe_sle) with a deglaciation history and compares FOUR
   !! integrations of the same load against a temporally converged truth:
   !!
   !!   (1) modal, 1 step/couple           (exact-exponential, what the model does)
   !!   (2) VE,    1 step/couple           (FE, single step)
   !!   (3) modal, N substeps/couple       (does sub-stepping move modal?)
   !!   (4) VE,    N substeps/couple   ==  TRUTH (temporally converged FE)
   !!
   !! Disambiguation:
   !!   * modal-1 ≈ truth  and  VE-1 is the outlier  ⇒ modal's exact-exponential is
   !!     ALREADY accurate; the ensemble "gap" is VE's own under-resolved stepping
   !!     (VE sub-steps 1.6–8.6×, modal once), i.e. VE-1's FE truncation, not modal.
   !!   * modal-1 lags truth and modal-N closes it ⇒ modal has a genuine temporal
   !!     error; the fix is sub-stepping or the linear-σ (FOH) modal advance.
   !!
   !! Rotation OFF (s_rot absent). Earth = M3-L70-V01 (= the ensemble radial earth).
   !!
   !!   usage:  diag_modal_sle.x [lmax] [nsub] [n_spin] [n_degla] [dt_yr]
   !!   default:                  16    8      60       200       100
   use fe_precision,       only: wp
   use fe_constants,       only: kyr, pi
   use fe_earth_structure, only: earth_model, build_M3L70V01
   use fe_radial_fe,       only: radial_fe_finalize
   use fe_response,        only: response, response_init_ve, response_init_modal, &
                                 response_set_dt, response_destroy
   use fe_sht,             only: sht_grid, sht_grid_init, sht_grid_destroy, &
                                 sht_grid_surface_integral
   use fe_sle,             only: sle_solve, sle_solver, sle_result
   use fe_timestep,        only: stepper_advance, adaptive_stepper
   implicit none

   integer            :: lmax, nsub, n_spin, n_degla, k
   real(wp)           :: dt, f_a, f_b, rtol
   character(len=64)  :: arg
   type(sht_grid)     :: sht
   type(earth_model)  :: e
   type(sle_result)   :: res
   real(wp), allocatable :: topo0(:,:), ice_lgm(:,:), d_ice(:,:), ice(:,:)
   ! five trajectories: VE-1, modal-1, VE-sub(truth), modal-sub, modal-ADAPTIVE (A3)
   type(response)   :: ve1, md1, veN, mdN, ma
   type(sle_solver) :: slv_ve1, slv_md1, slv_veN, slv_mdN, slv_ma
   type(adaptive_stepper) :: stp_ma
   real(wp), allocatable :: Sve1(:,:), Smd1(:,:), SveN(:,:), SmdN(:,:), Sma(:,:)
   real(wp), allocatable :: Cve1(:,:), Cmd1(:,:), CveN(:,:), CmdN(:,:), Cma(:,:)
   real(wp), allocatable :: Struth0(:,:), ice_a(:,:), ice_b(:,:), ice_ref0(:,:)
   real(wp) :: tcur

   lmax = 16;  nsub = 8;  n_spin = 60;  n_degla = 200;  dt = 100.0_wp;  rtol = 1.0e-4_wp
   if (command_argument_count() >= 1) then; call get_command_argument(1,arg); read(arg,*) lmax;    end if
   if (command_argument_count() >= 2) then; call get_command_argument(2,arg); read(arg,*) nsub;    end if
   if (command_argument_count() >= 3) then; call get_command_argument(3,arg); read(arg,*) n_spin;  end if
   if (command_argument_count() >= 4) then; call get_command_argument(4,arg); read(arg,*) n_degla; end if
   if (command_argument_count() >= 5) then; call get_command_argument(5,arg); read(arg,*) dt;      end if
   if (command_argument_count() >= 6) then; call get_command_argument(6,arg); read(arg,*) rtol;    end if
   dt = dt*1.0e-3_wp*kyr

   call sht_grid_init(sht, lmax, nlat=2*lmax, nphi=4*lmax)
   e = build_M3L70V01()
   allocate(topo0(sht%nphi,sht%nlat), ice_lgm(sht%nphi,sht%nlat))
   allocate(d_ice(sht%nphi,sht%nlat), ice(sht%nphi,sht%nlat))
   allocate(Sve1(sht%nphi,sht%nlat), Smd1(sht%nphi,sht%nlat), &
            SveN(sht%nphi,sht%nlat), SmdN(sht%nphi,sht%nlat), Struth0(sht%nphi,sht%nlat))
   allocate(Cve1(sht%nphi,sht%nlat), Cmd1(sht%nphi,sht%nlat), &
            CveN(sht%nphi,sht%nlat), CmdN(sht%nphi,sht%nlat))
   allocate(Sma(sht%nphi,sht%nlat), Cma(sht%nphi,sht%nlat), &
            ice_a(sht%nphi,sht%nlat), ice_b(sht%nphi,sht%nlat), ice_ref0(sht%nphi,sht%nlat))
   ice_ref0 = 0.0_wp                       ! ice-free reference
   call make_fields(topo0, ice_lgm)

   call response_init_ve(ve1, e, sht, dt)
   call response_init_ve(veN, e, sht, dt)
   call response_init_modal(md1, e, sht, n_modes=-1, mode_rank=1, dt_be=5.0_wp*kyr)
   call response_init_modal(mdN, e, sht, n_modes=-1, mode_rank=1, dt_be=5.0_wp*kyr)
   call response_init_modal(ma,  e, sht, n_modes=-1, mode_rank=1, dt_be=5.0_wp*kyr)
   ma%modal_adaptive = .true.              ! opt this trajectory into A3 (off by default)
   tcur = 0.0_wp                           ! running clock for the adaptive stepper
   stp_ma%rtol = rtol                      ! A3 accuracy tolerance (the cost knob)

   write(*,'(a)') ' === modal vs VE through the SLE: temporal attribution — M3-L70-V01, deglaciation ==='
   write(*,'(a,i0,a,i0,a,i0,a,i0,a,f6.1,a,es8.1,a)') '   lmax=', lmax, '  nsub=', nsub, &
        '  n_spin=', n_spin, '  n_degla=', n_degla, '  dt=', dt/kyr*1.0e3_wp, ' yr  rtol=', rtol, '  (rot off)'
   write(*,'(a)') ''

   ! ---- LGM spin-up: hold the cap so each trajectory relaxes toward equilibrium --
   do k = 1, n_spin
      call step_couple(ve1, slv_ve1, Sve1, Cve1, 1,    1.0_wp, 1.0_wp)
      call step_couple(md1, slv_md1, Smd1, Cmd1, 1,    1.0_wp, 1.0_wp)
      call step_couple(veN, slv_veN, SveN, CveN, nsub, 1.0_wp, 1.0_wp)
      call step_couple(mdN, slv_mdN, SmdN, CmdN, nsub, 1.0_wp, 1.0_wp)
      call step_adaptive(1.0_wp, 1.0_wp)
   end do
   Struth0 = SveN                         ! truth rsl entering the deglaciation

   ! ---- deglaciation: unload the cap linearly to zero -----------------------
   do k = 1, n_degla
      f_a = 1.0_wp - real(k-1,wp)/real(n_degla,wp)
      f_b = 1.0_wp - real(k,  wp)/real(n_degla,wp)
      call step_couple(ve1, slv_ve1, Sve1, Cve1, 1,    f_a, f_b)
      call step_couple(md1, slv_md1, Smd1, Cmd1, 1,    f_a, f_b)
      call step_couple(veN, slv_veN, SveN, CveN, nsub, f_a, f_b)
      call step_couple(mdN, slv_mdN, SmdN, CmdN, nsub, f_a, f_b)
      call step_adaptive(f_a, f_b)
   end do

   ! ---- report: each trajectory vs truth (VE, nsub substeps) ----------------
   write(*,'(a)') ' Area-weighted rsl gap vs TRUTH (VE, N substeps/couple) at present day [m]:'
   write(*,'(a)') '   trajectory                       gap rms        gap max'
   write(*,'(a,es15.3,es15.3)') '   modal, 1 step/couple   (model) ', warea_rms(Smd1-SveN), maxval(abs(Smd1-SveN))
   write(*,'(a,es15.3,es15.3)') '   VE,    1 step/couple           ', warea_rms(Sve1-SveN), maxval(abs(Sve1-SveN))
   write(*,'(a,es15.3,es15.3)') '   modal, N substeps/couple       ', warea_rms(SmdN-SveN), maxval(abs(SmdN-SveN))
   write(*,'(a,es15.3,es15.3,a,f5.2,a)') '   modal, ADAPTIVE (A3)           ', &
        warea_rms(Sma-SveN), maxval(abs(Sma-SveN)), &
        '   [', real(stp_ma%n_solve,wp)/real(n_spin+n_degla,wp), ' SLE/couple avg]'
   write(*,'(a,es15.3,a)')      '   (truth relax span over run     ', warea_rms(SveN-Struth0), ' m)'
   write(*,'(a)') ''
   write(*,'(a)') ' Read: A3 (adaptive) should match modal-N accuracy by sub-stepping only where σ'
   write(*,'(a)') ' varies (≈1 SLE/couple in spin-up, more during deglaciation), driven by rtol.'

   call response_destroy(ve1);  call response_destroy(md1)
   call response_destroy(veN);  call response_destroy(mdN);  call response_destroy(ma)
   call sht_grid_destroy(sht);  call radial_fe_finalize()

contains

   subroutine step_adaptive(f_a, f_b)
      !! Advance one coupling interval through the adaptive controller (A3): the
      !! stepper ramps the ice f_a·cap → f_b·cap over [tcur, tcur+dt] and sub-steps
      !! modal to rtol via its own step-doubling. Advances the running clock tcur.
      real(wp), intent(in) :: f_a, f_b
      ice_a = f_a*ice_lgm;  ice_b = f_b*ice_lgm
      call stepper_advance(stp_ma, sht, ma, slv_ma, topo0, ice_a, ice_b, ice_ref0, &
                           tcur, tcur + dt, Sma, Cma)
      tcur = tcur + dt
   end subroutine step_adaptive

   subroutine step_couple(resp, sle, S, C, nss, f_a, f_b)
      !! Advance one coupling interval in nss substeps, the ice load ramping
      !! linearly from f_a·cap to f_b·cap (sampled at each substep endpoint).
      type(response),   intent(inout) :: resp
      type(sle_solver), intent(inout) :: sle
      real(wp),         intent(inout) :: S(:,:), C(:,:)
      integer,          intent(in)    :: nss
      real(wp),         intent(in)    :: f_a, f_b
      integer  :: iss
      real(wp) :: f
      call response_set_dt(resp, dt/real(nss,wp))
      do iss = 1, nss
         f = f_a + (f_b - f_a)*real(iss,wp)/real(nss,wp)
         d_ice = f*ice_lgm;  ice = d_ice
         call sle_solve(sle, sht, resp, d_ice, ice, topo0, S, C, res)
      end do
   end subroutine step_couple

   real(wp) function warea_rms(d) result(r)
      !! Area-weighted (cos lat) RMS of a grid field.
      real(wp), intent(in) :: d(:,:)
      r = sqrt(sht_grid_surface_integral(sht, d*d) / (4.0_wp*pi))
   end function warea_rms

   subroutine make_fields(topo0, ice_lgm)
      !! Land cap (colat<60°, +500 m) over ocean (−4000 m); a 2 km grounded LGM ice
      !! cap (colat<40°), unloaded to zero over the deglaciation.
      real(wp), intent(out) :: topo0(:,:), ice_lgm(:,:)
      integer  :: i, j
      real(wp) :: th
      do j = 1, sht%nlat
         th = sht%colat(j)
         do i = 1, sht%nphi
            topo0(i,j)   = merge(500.0_wp, -4000.0_wp, th < 60.0_wp*pi/180.0_wp)
            ice_lgm(i,j) = merge(2000.0_wp, 0.0_wp,    th < 40.0_wp*pi/180.0_wp)
         end do
      end do
   end subroutine make_fields

end program diag_modal_sle
