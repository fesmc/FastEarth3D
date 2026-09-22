program diag_visc3d_paths
   !! Diagnostic (session 32): the 1-D and 3-D memory-advance paths side by side on
   !! the PRODUCTION configuration, with the memory readable at every step.
   !!
   !! test_response_3d drives both paths on M3-L70-V01 at lmax 8 with ONE uniform
   !! perturbation p applied to every element, and they agree to 1e-13. The block D
   !! runs use the 11-layer PREM table on its 220-element mesh at lmax 170, and a
   !! perturbation that VARIES from about 0 to +9.4 dex across the mesh because the
   !! viscosity file and the namelist eta column are different Earths. On that
   !! configuration the two paths disagree by 4e-5 relative in rsl after ONE coupling
   !! step (LOG.md session 32i). Nothing in the unit harness reproduces it.
   !!
   !! Both responses read the SAME visc_node, through the same driver call
   !! (load_visc_3d + response_enable_lateral_visc_from_nodes); they differ only in
   !! visc3d_tol -- 99 demotes every element to the scalar path, -1 forces every
   !! Maxwell element through the pseudo-spectral kernel. With a laterally uniform
   !! file the two must agree.
   !!
   !! The degree-1 frame is an argument because it is NOT neutral here: Z6 has no
   !! degree-1 harmonic (norm6 = 2*Jr*(Jr-2) = 0 at l=1), so the scalar path keeps a
   !! phantom lambda=6 memory at l=1 that the tensor path correctly zeros. In the
   !! "cf" frame rsl carries no degree 1 at all and the phantom is invisible -- which
   !! is the frame every benchmark and test_response_3d runs in. Block D production
   !! runs "cm", where degree 1 is retained.
   !!
   !!   diag_visc3d_paths.x <visc_file.nc> [lmax] [nsteps] [cf|cm] [skip_tol]
   !!
   !! skip_tol gates solve_drift: a coefficient whose mnorm falls below
   !! skip_tol*maxval(mnorm) has its drift zeroed instead of solved. mnorm is a RAW
   !! max over the four lambda channels, unweighted by the Martinec norms
   !! nrmc = [1, Jr/2, 2Jr^2, 2Jr(Jr-2)] that dissipative_rhs actually applies. Pass
   !! 0 to disable skipping entirely.
   use fe_precision,       only: wp
   use fe_constants,       only: kyr
   use fe_params,          only: fe_param_class
   use fe_earth_structure, only: earth_model, build_earth, load_visc_3d
   use fe_radial_fe,       only: radial_fe_finalize
   use fe_response,        only: response, response_init_ve, response_destroy, &
                                 response_enable_lateral_visc_from_nodes, &
                                 response_begin_step, response_apply, response_commit_step
   use fe_sht,             only: sht_grid, sht_grid_init, sht_grid_lmidx, sht_grid_destroy
   implicit none

   type(fe_param_class) :: p
   type(sht_grid)       :: sht
   type(earth_model)    :: em
   type(response)       :: v1, v3
   real(wp), allocatable :: visc_node(:,:)
   complex(wp), allocatable :: slm(:), u1(:), n1(:), u3(:), n3(:)
   character(len=512) :: arg, viscfile
   character(len=8)   :: frame
   real(wp) :: dt
   logical  :: kill_l1
   integer  :: lmax, nsteps, i, lam, e, k, nlat, nphi
   integer  :: ne3d_1, ne3d_3

   call get_command_argument(1, arg);  viscfile = trim(arg)
   if (len_trim(viscfile) == 0) error stop 'usage: diag_visc3d_paths.x <visc_file.nc> [lmax] [nsteps]'
   lmax = 170;  nsteps = 5
   call get_command_argument(2, arg);  if (len_trim(arg) > 0) read(arg,*) lmax
   call get_command_argument(3, arg);  if (len_trim(arg) > 0) read(arg,*) nsteps
   frame = 'cf';  kill_l1 = .false.
   call get_command_argument(4, arg);  if (len_trim(arg) > 0) frame = trim(arg)
   call get_command_argument(5, arg);  if (len_trim(arg) > 0) kill_l1 = (trim(arg) == 'kill')

   ! --- block D's Earth, verbatim from experiments/timing/blockD.nml.tmpl --------
   p%earth   = "PREM"
   p%n_layer = 11
   p%r_earth = 6371.0e3_wp;  p%r_core = 3480.0e3_wp
   p%r_bot(1:11) = [6291.0e3_wp, 6151.0e3_wp, 5971.0e3_wp, 5771.0e3_wp, 5701.0e3_wp, &
                    5600.0e3_wp, 4943.0e3_wp, 4287.0e3_wp, 3630.0e3_wp, 3480.0e3_wp, 0.0_wp]
   p%r_top(1:11) = [6371.0e3_wp, 6291.0e3_wp, 6151.0e3_wp, 5971.0e3_wp, 5771.0e3_wp, &
                    5701.0e3_wp, 5600.0e3_wp, 4943.0e3_wp, 4287.0e3_wp, 3630.0e3_wp, 3480.0e3_wp]
   p%eta(1:11)   = [1.0e40_wp, 4.0e20_wp, 4.0e20_wp, 4.0e20_wp, 4.0e20_wp, &
                    1.0e22_wp, 1.0e22_wp, 1.0e22_wp, 1.0e22_wp, 1.0e22_wp, 0.0_wp]
   p%rheology(1:11) = [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2]
   p%visc_3d_file = viscfile

   nlat = 2*(lmax + 1);  nphi = 4*lmax           ! the driver's production grid
   call sht_grid_init(sht, lmax, nlat=nlat, nphi=nphi, mmax=lmax)
   em = build_earth(p)
   dt = 0.05_wp*kyr                              ! 50 yr, the production sub-step
   ! deg1_cm must be set BEFORE init: the elastic and viscoelastic gains are computed
   ! there and the frame is part of them (fe_coupling does the same).
   v1%deg1_cm = (trim(frame) == 'cm');  v3%deg1_cm = v1%deg1_cm
   call response_init_ve(v1, em, sht, dt)
   call response_init_ve(v3, em, sht, dt)

   call load_visc_3d(p, sht, v1%r, visc_node)
   v1%visc3d_tol =  99.0_wp                      ! every element demoted  -> 1-D path
   v3%visc3d_tol =  -1.0_wp                      ! every Maxwell element  -> 3-D path
   call response_enable_lateral_visc_from_nodes(v1, sht, visc_node)
   call response_enable_lateral_visc_from_nodes(v3, sht, visc_node)
   ne3d_1 = v1%ne3d;  ne3d_3 = v3%ne3d

   write(*,'(a,i0,a,i0,a,i0,a,a)') ' lmax=', lmax, '  nlat=', sht%nlat, '  nphi=', sht%nphi, &
        '  deg1_frame=', trim(frame)
   write(*,'(a,es9.2)') ' skip_tol=', v1%skip_tol
   write(*,'(a,i0,a,i0,a)')    ' ne=', v1%ne, '   Maxwell elements promoted: 1-D path ', ne3d_1, ','
   write(*,'(a,i0)')           '                                            3-D path ', ne3d_3
   write(*,'(a,es10.3,a,es10.3)') ' lateral spread of the loaded field [dex]: ', &
        maxval(visc_node) - minval(visc_node), '   max M = ', maxval(v1%Mk)
   write(*,'(a)') ''

   allocate(slm(sht%nlm), u1(sht%nlm), n1(sht%nlm), u3(sht%nlm), n3(sht%nlm))
   call build_load(slm)
   write(*,'(a5,a14,4a13,a14)') 'step', 'rel uplift', 'lam1', 'lam2', 'lam5', 'lam6', 'worst elem'
   do i = 1, nsteps
      call response_begin_step(v1, sht);  call response_apply(v1, sht, slm, u1, n1)
      call response_begin_step(v3, sht);  call response_apply(v3, sht, slm, u3, n3)
      call report(i, u1, u3)
      call response_commit_step(v1, sht, slm)
      call response_commit_step(v3, sht, slm)
      ! Optional: annihilate the lambda=6 memory at l<2 in the SCALAR path, which is
      ! the one place the two paths differ by construction. Z6 has no degree-1
      ! harmonic -- ve_strain_constants gives it norm6 = 2*Jr*(Jr-2) = 0 at l=1 -- so
      ! the tensor path zeros it while the scalar advance keeps a phantom value.
      ! test_response_3d masks exactly this slot and calls it "the null space of the
      ! observable". If that is true the line below changes nothing; if the cm
      ! degree-1 term reads the memory without the norm, it is not true in cm.
      if (kill_l1) call zero_l1_lam6(v1)
   end do

   call response_destroy(v1);  call response_destroy(v3)
   call sht_grid_destroy(sht);  call radial_fe_finalize()

contains

   subroutine report(i, ua, ub)
      integer,     intent(in) :: i
      complex(wp), intent(in) :: ua(:), ub(:)
      real(wp) :: dl(4), sl(4), du, su, dbest
      integer  :: ebest
      do lam = 1, 4
         dl(lam) = max(chan(v1%Bre, v3%Bre, lam), chan(v1%Bim, v3%Bim, lam), &
                       chan(v1%Cre, v3%Cre, lam), chan(v1%Cim, v3%Cim, lam), &
                       chan(v1%Are, v3%Are, lam), chan(v1%Aim, v3%Aim, lam))
         sl(lam) = max(maxval(abs(v1%Bre(lam,:,:))), maxval(abs(v1%Cre(lam,:,:))), &
                       maxval(abs(v1%Are(lam,:,:))), tiny(1.0_wp))
      end do
      du = maxval(abs(ua - ub));  su = max(maxval(abs(ua)), tiny(1.0_wp))
      dbest = 0.0_wp;  ebest = 0
      do e = 1, v1%ne
         if (maxval(abs(v1%Bre(:,e,:) - v3%Bre(:,e,:))) > dbest) then
            dbest = maxval(abs(v1%Bre(:,e,:) - v3%Bre(:,e,:)));  ebest = e
         end if
      end do
      write(*,'(i5,es14.3,4es13.3,i8,a,es9.2)') i, du/su, &
           dl(1)/sl(1), dl(2)/sl(2), dl(3)/sl(3), dl(4)/sl(4), ebest, '  M=', v1%Mk(max(ebest,1))
   end subroutine report

   subroutine zero_l1_lam6(v)
      type(response), intent(inout) :: v
      integer :: kk
      do kk = 1, v%nk
         if (v%kdeg(kk) >= 2) cycle
         v%Are(4,:,kk) = 0.0_wp;  v%Aim(4,:,kk) = 0.0_wp
         v%Bre(4,:,kk) = 0.0_wp;  v%Bim(4,:,kk) = 0.0_wp
         v%Cre(4,:,kk) = 0.0_wp;  v%Cim(4,:,kk) = 0.0_wp
      end do
   end subroutine zero_l1_lam6

   real(wp) function chan(a, b, lam) result(d)
      real(wp), intent(in) :: a(:,:,:), b(:,:,:)
      integer,  intent(in) :: lam
      real(wp), allocatable :: da(:,:)
      da = a(lam,:,:) - b(lam,:,:)
      ! mask the physically inert lambda=6 degree-1 slot (norm6 = 2*Jr*(Jr-2) = 0 at
      ! l=1): the scalar path keeps a phantom value there, the tensor path zeros it,
      ! and the dissipation never sees either. Same masking as test_response_3d.
      if (lam == 4) then
         do k = 1, v1%nk
            if (v1%kdeg(k) < 2) da(:,k) = 0.0_wp
         end do
      end if
      d = maxval(abs(da))
   end function chan

   subroutine build_load(s)
      complex(wp), intent(out) :: s(:)
      integer :: l, m, lm
      s = (0.0_wp, 0.0_wp)
      do l = 1, min(lmax, 64)
         do m = 0, l
            lm = sht_grid_lmidx(sht, l, m)
            if (m == 0) then
               s(lm) = cmplx(1000.0_wp/real(l*l, wp), 0.0_wp, wp)
            else
               s(lm) = cmplx(700.0_wp/real(l*l, wp), 400.0_wp/real(l*(m+1), wp), wp)
            end if
         end do
      end do
   end subroutine build_load

end program diag_visc3d_paths
