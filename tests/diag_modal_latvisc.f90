program diag_modal_latvisc
   !! DIAGNOSTIC (not pass/fail): isolate the modal LATERAL-VISCOSITY error against the
   !! VE tensor-SH ground truth, with lateral viscosity ACTIVE on BOTH paths, AND
   !! attribute it across the three modal lateral methods:
   !!
   !!   LIE     — 1st-order Lie split, rank-characteristic τ̂ (the original production path)
   !!   STRANG  — symmetric 2nd-order split, same τ̂  (isolates the split-ORDER error)
   !!   COUPLED — per-rank coupled rate operator exp(−Δt·L_i), no split, per-degree-exact τ
   !!             (the accurate path; remaining gap = modal-basis + characteristic-weight)
   !!
   !! Earlier modal-vs-VE probes were all RADIAL and missed this entirely. A controlled
   !! polar viscosity CAP (depth-uniform Δ=log10(η/η_ref), 5° taper) + co-located disc
   !! load drives modal+LV vs VE+LV through the IDENTICAL pert_elem — only the algorithm
   !! differs. Free-rebound protocol: glaciate (hold load) then unload (σ=0) and watch
   !! free relaxation, whose decay rate IS the lateral operator (fluid limits are
   !! contrast-independent, so any gap is a transient RATE error).
   !!
   !! Per (theta_cap, Δ) the table reports, for each method: end-glaciation modal/VE
   !! cap-peak ratio (1.0 = perfect), field error grms/pk (global area-rms gap over the VE
   !! peak), and the rebound peak ratio at 5 and 20 kyr. Rotation OFF; M3-L70-V01.
   !!
   !!   usage:  diag_modal_latvisc.x [lmax] [dt_yr] [n_glac] [n_rebound] [m_krylov]
   !!   default:                      48    200     200      100        12
   use fe_precision,       only: wp
   use fe_constants,       only: kyr, pi
   use fe_earth_structure, only: earth_model, build_M3L70V01
   use fe_radial_fe,       only: radial_fe_finalize
   use fe_response,        only: response, response_init_ve, response_init_modal, &
                                 response_set_dt, response_destroy, &
                                 response_begin_step, response_apply, response_commit_step, &
                                 response_enable_lateral_visc, response_enable_lateral_visc_modal, &
                                 LAT_LIE_CHAR, LAT_STRANG_CHAR, LAT_COUPLED
   use fe_sht,             only: sht_grid, sht_grid_init, sht_grid_destroy, &
                                 sht_grid_synthesis, sht_grid_analysis, &
                                 sht_grid_surface_integral
   implicit none

   integer            :: lmax, n_glac, n_rebound, ic, it, m_krylov
   real(wp)           :: dt
   character(len=64)  :: arg
   type(sht_grid)     :: sht
   type(earth_model)  :: e
   complex(wp), allocatable :: load_lm(:)
   real(wp),    allocatable :: load_g(:,:)
   ! sweep grids
   integer,  parameter :: NC = 4, NT = 3
   real(wp), parameter :: contrast(NC)  = [0.5_wp, 1.0_wp, 1.5_wp, 2.0_wp]
   real(wp), parameter :: thcap_deg(NT) = [10.0_wp, 20.0_wp, 40.0_wp]
   ! rebound checkpoints (kyr after unload)
   integer,  parameter :: NCHK = 2
   real(wp) :: chk_kyr(NCHK) = [5.0_wp, 20.0_wp]
   integer  :: methods(3) = [LAT_LIE_CHAR, LAT_STRANG_CHAR, LAT_COUPLED]
   character(len=7) :: mname(3) = ['lie    ', 'strang ', 'coupled']

   lmax = 48;  dt = 200.0_wp;  n_glac = 200;  n_rebound = 100;  m_krylov = 12
   if (command_argument_count() >= 1) then; call get_command_argument(1,arg); read(arg,*) lmax;      end if
   if (command_argument_count() >= 2) then; call get_command_argument(2,arg); read(arg,*) dt;        end if
   if (command_argument_count() >= 3) then; call get_command_argument(3,arg); read(arg,*) n_glac;    end if
   if (command_argument_count() >= 4) then; call get_command_argument(4,arg); read(arg,*) n_rebound; end if
   if (command_argument_count() >= 5) then; call get_command_argument(5,arg); read(arg,*) m_krylov;  end if
   dt = dt*1.0e-3_wp*kyr

   call sht_grid_init(sht, lmax, nlat=2*lmax, nphi=4*lmax)
   e = build_M3L70V01()
   allocate(load_lm(sht%nlm), load_g(sht%nphi,sht%nlat))

   write(*,'(a)')  ' === modal+LV (lie/strang/coupled) vs VE+LV (tensor-SH) — M3-L70-V01 ==='
   write(*,'(a,i0,a,f6.1,a,i0,a,f6.1,a,i0,a,f6.1,a)') &
        '   lmax=', lmax, '  dt=', dt/kyr*1.0e3_wp, ' yr   glaciation=', n_glac, &
        ' steps (', n_glac*dt/kyr, ' kyr)   rebound=', n_rebound, ' steps (', n_rebound*dt/kyr, ' kyr)'
   write(*,'(a,i0)') '   stiff cap Δ=log10(η/η_ref)>0 over colat<theta_cap; co-located disc load.  m_krylov=', m_krylov
   write(*,'(a)')  '   ratio = modal/VE cap-peak (1.0 perfect); grms/pk = field gap over VE peak.'
   write(*,'(a)')  ''

   do ic = 1, NT
      call make_disc_load(sht, thcap_deg(ic)*pi/180.0_wp, load_g)
      call sht_grid_analysis(sht, load_g, load_lm)
      write(*,'(a,f5.1,a)') ' --- theta_cap = theta_load = ', thcap_deg(ic), ' deg ----------------------------'
      write(*,'(a)') '   Δdex  method    endglac-ratio   grms/pk  | rebound 5kyr  20kyr'
      do it = 1, NC
         call run_config(thcap_deg(ic), contrast(it), load_lm)
      end do
      write(*,'(a)') ''
   end do

   call sht_grid_destroy(sht);  call radial_fe_finalize()
   deallocate(load_lm, load_g)

contains

   subroutine run_config(thcap_deg_, delta, load_lm)
      !! For one (theta_cap, Δ): run VE+LV once (truth), then modal+LV under each of the
      !! three lateral methods, and print one row per method comparing to the same VE.
      real(wp),    intent(in) :: thcap_deg_, delta
      complex(wp), intent(in) :: load_lm(:)
      type(response) :: ve, md
      real(wp), allocatable :: pert(:,:,:), veUg_gl(:,:)
      real(wp) :: theta_cap, ve_pk_gl, ve_pk_chk(NCHK)
      real(wp) :: md_pk_gl, md_pk_chk(NCHK), grms
      integer  :: im, ichk

      theta_cap = thcap_deg_*pi/180.0_wp

      ! ---- VE+LV ground truth ------------------------------------------------
      call response_init_ve(ve, e, sht, dt);  call response_set_dt(ve, dt)
      allocate(pert(sht%nphi, sht%nlat, ve%ne))
      call make_cap_pert(sht, theta_cap, delta, pert)
      call response_enable_lateral_visc(ve, sht, pert)
      allocate(veUg_gl(sht%nphi, sht%nlat))
      call run_traj(ve, load_lm, theta_cap, ve_pk_gl, veUg_gl, ve_pk_chk)
      call response_destroy(ve)

      ! ---- modal+LV, three methods, each vs the same VE ----------------------
      call response_init_modal(md, e, sht, n_modes=-1, mode_rank=1, dt_be=5.0_wp*kyr)
      call response_set_dt(md, dt)
      call response_enable_lateral_visc_modal(md, sht, pert)
      do im = 1, 3
         md%phi = (0.0_wp, 0.0_wp);  md%phi_n = (0.0_wp, 0.0_wp)   ! reset to relaxed reference
         md%lat_method = methods(im);  md%m_krylov_lat = m_krylov
         call run_traj_md(md, load_lm, theta_cap, veUg_gl, md_pk_gl, grms, md_pk_chk)
         write(*,'(f6.2,2x,a,2x,f10.3,4x,f8.4,a,2(2x,f7.3))') &
              delta, mname(im), md_pk_gl/max(abs(ve_pk_gl),tiny(1.0_wp)), &
              grms/max(abs(ve_pk_gl),tiny(1.0_wp)), '  |', &
              (md_pk_chk(ichk)/max(abs(ve_pk_chk(ichk)),tiny(1.0_wp)), ichk=1,NCHK)
      end do
      call response_destroy(md)
      deallocate(pert, veUg_gl)
   end subroutine run_config

   subroutine run_traj(resp, load_lm, theta_cap, pk_gl, Ug_gl, pk_chk)
      !! Drive one response through glaciation (hold load) + free rebound (σ=0). Returns
      !! the cap-peak uplift and full uplift grid at end of glaciation, and the cap-peak
      !! at the rebound checkpoints.
      type(response), intent(inout) :: resp
      complex(wp),    intent(in)    :: load_lm(:)
      real(wp),       intent(in)    :: theta_cap
      real(wp),       intent(out)   :: pk_gl, Ug_gl(:,:), pk_chk(NCHK)
      complex(wp), allocatable :: zero_lm(:)
      integer :: k, ichk, nchk_step(NCHK)
      allocate(zero_lm(sht%nlm));  zero_lm = (0.0_wp, 0.0_wp)
      do k = 1, n_glac
         call step(resp, load_lm)
      end do
      call ugrid(resp, Ug_gl)
      call peak_over_cap(Ug_gl, theta_cap, pk_gl)
      do ichk = 1, NCHK
         nchk_step(ichk) = nint(chk_kyr(ichk)*kyr/dt)
      end do
      pk_chk = -1.0_wp
      do k = 1, n_rebound
         call step(resp, zero_lm)
         do ichk = 1, NCHK
            if (k == nchk_step(ichk)) call peak_over_cap_resp(resp, theta_cap, pk_chk(ichk))
         end do
      end do
      deallocate(zero_lm)
   end subroutine run_traj

   subroutine run_traj_md(resp, load_lm, theta_cap, veUg_gl, pk_gl, grms, pk_chk)
      !! Same as run_traj but also returns the global area-rms gap vs the VE end-glaciation
      !! uplift field (modal only).
      type(response), intent(inout) :: resp
      complex(wp),    intent(in)    :: load_lm(:)
      real(wp),       intent(in)    :: theta_cap, veUg_gl(:,:)
      real(wp),       intent(out)   :: pk_gl, grms, pk_chk(NCHK)
      real(wp), allocatable :: Ug_gl(:,:)
      allocate(Ug_gl(sht%nphi, sht%nlat))
      call run_traj(resp, load_lm, theta_cap, pk_gl, Ug_gl, pk_chk)
      grms = sqrt(sht_grid_surface_integral(sht, (Ug_gl-veUg_gl)**2) / (4.0_wp*pi))
      deallocate(Ug_gl)
   end subroutine run_traj_md

   subroutine step(resp, sigma_lm)
      !! One held-load step: apply (read state) then commit (advance memory).
      type(response), intent(inout) :: resp
      complex(wp),    intent(in)    :: sigma_lm(:)
      complex(wp), allocatable :: u_lm(:), n_lm(:)
      allocate(u_lm(sht%nlm), n_lm(sht%nlm))
      call response_begin_step(resp, sht)
      call response_apply(resp, sht, sigma_lm, u_lm, n_lm)
      call response_commit_step(resp, sht, sigma_lm)
      deallocate(u_lm, n_lm)
   end subroutine step

   subroutine ugrid(resp, ug)
      !! Current uplift field on the grid (reads frozen state; load irrelevant to apply).
      type(response), intent(inout) :: resp
      real(wp),       intent(out)   :: ug(:,:)
      complex(wp), allocatable :: u_lm(:), n_lm(:), z(:)
      allocate(u_lm(sht%nlm), n_lm(sht%nlm), z(sht%nlm));  z = (0.0_wp, 0.0_wp)
      call response_begin_step(resp, sht)
      call response_apply(resp, sht, z, u_lm, n_lm)
      call sht_grid_synthesis(sht, u_lm, ug)
      deallocate(u_lm, n_lm, z)
   end subroutine ugrid

   subroutine peak_over_cap_resp(resp, theta_cap, pk)
      type(response), intent(inout) :: resp
      real(wp),       intent(in)    :: theta_cap
      real(wp),       intent(out)   :: pk
      real(wp), allocatable :: ug(:,:)
      allocate(ug(sht%nphi,sht%nlat))
      call ugrid(resp, ug)
      call peak_over_cap(ug, theta_cap, pk)
      deallocate(ug)
   end subroutine peak_over_cap_resp

   subroutine peak_over_cap(ug, theta_cap, pk)
      !! Max |uplift| over the cap region (colat < theta_cap).
      real(wp), intent(in)  :: ug(:,:)
      real(wp), intent(in)  :: theta_cap
      real(wp), intent(out) :: pk
      integer :: i, j
      pk = 0.0_wp
      do j = 1, sht%nlat
         if (sht%colat(j) >= theta_cap) cycle
         do i = 1, sht%nphi
            pk = max(pk, abs(ug(i,j)))
         end do
      end do
   end subroutine peak_over_cap

   subroutine make_disc_load(sht, theta_load, g)
      !! Axisymmetric disc load: 1 inside colat<theta_load, 5° cosine taper to 0.
      type(sht_grid), intent(in)  :: sht
      real(wp),       intent(in)  :: theta_load
      real(wp),       intent(out) :: g(:,:)
      real(wp) :: taper, th, x
      integer  :: i, j
      taper = 5.0_wp*pi/180.0_wp
      do j = 1, sht%nlat
         th = sht%colat(j)
         if (th < theta_load) then;            x = 1.0_wp
         else if (th < theta_load + taper) then; x = 0.5_wp*(1.0_wp + cos(pi*(th-theta_load)/taper))
         else;                                 x = 0.0_wp
         end if
         do i = 1, sht%nphi
            g(i,j) = x
         end do
      end do
   end subroutine make_disc_load

   subroutine make_cap_pert(sht, theta_cap, delta, pert)
      !! Depth-uniform polar viscosity cap: pert = Δ inside colat<theta_cap (5° taper), 0
      !! outside, identical across elements (non-Maxwell elements ignored by both enablers).
      type(sht_grid), intent(in)  :: sht
      real(wp),       intent(in)  :: theta_cap, delta
      real(wp),       intent(out) :: pert(:,:,:)
      real(wp) :: taper, th, x
      integer  :: i, j
      taper = 5.0_wp*pi/180.0_wp
      do j = 1, sht%nlat
         th = sht%colat(j)
         if (th < theta_cap) then;            x = 1.0_wp
         else if (th < theta_cap + taper) then; x = 0.5_wp*(1.0_wp + cos(pi*(th-theta_cap)/taper))
         else;                                 x = 0.0_wp
         end if
         do i = 1, sht%nphi
            pert(i,j,:) = delta*x
         end do
      end do
   end subroutine make_cap_pert

end program diag_modal_latvisc
