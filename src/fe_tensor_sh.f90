module fe_tensor_sh
   !! Tensor spherical-harmonic dyadic transforms (Martinec 2000, Appendix B), the
   !! machinery rung 6 (laterally-varying viscosity) needs — GENERAL order (mmax≥0).
   !!
   !! The Maxwell memory τ and strain ε are second-order symmetric tensors on the
   !! sphere, expanded in the spheroidal tensor spherical harmonics Z^λ, λ∈{1,2,5,6}
   !! (toroidal λ=3,4 drop). For lateral viscosity the update τ⁺=(1−M)τ−2μM·ε is
   !! pointwise in PHYSICAL space, so the tensor is reconstructed on the grid via its
   !! six dyadic components (eqs 90/91, B10/B11) and projected back. The six physical
   !! components rr, rθ, rφ, θθ, θφ, φφ relate to the coefficients by
   !!   rr        = Σ T¹ Y                          (Z¹)
   !!   rθ, rφ    = Σ T² (E, F)                      (Z²,  E=∂_θY, F=(1/sinθ)∂_φY)
   !!   θθ        = Σ [−l(l+1)T⁵ Y + T⁶ G]          (Z⁵ trace + Z⁶)
   !!   φφ        = Σ [−l(l+1)T⁵ Y − T⁶ G]
   !!   θφ        = Σ 4 T⁶ H                          (Z⁶,  G,H = B11 second derivatives)
   !!
   !! Synthesis is EXACT via grid identities — no recurrence, no re-analysis. The
   !! spin-2 G,H (the only pieces SHTns has no routine for) come from scalar + vector
   !! synths plus algebraic grid factors, using ∂_φ = "multiply coeffs by im" (exact
   !! on the known input) and ∂_θθ via ∇₁²Y = −l(l+1)Y:
   !!   Sg ≡ Σ T⁶ G = ∇₁²f − 2cotθ·g_θ − 2(1/sinθ)·∂_φ g_φ
   !!   Sh ≡ Σ T⁶ H = (1/sinθ)·∂_φ g_θ − cotθ·g_φ
   !! with f=synth(T⁶), (g_θ,g_φ)=sph_synth(T⁶), ∂_φ(·)=sph_synth(im·T⁶).
   !!
   !! Analysis: channels 1,2,5 invert through SHTns's own scalar/vector analyses (a
   !! synth/analysis pair) with the −l(l+1) factor for the trace; the spin-2 channel 6
   !! uses the adjoint of its synthesis (`Sg*`,`Sh*`) — the same ops with synth↔analysis
   !! swapped — normalised by a per-degree factor calibrated once at init. Validated by
   !! the round trip and the physical ∫τ:ε double-dot vs the B13 norms (test_tensor_sh).
   use fe_precision, only: wp
   use fe_constants, only: pi
   use fe_sht,       only: sht_grid, sht_free_cfg, sht_grid_lmidx, sht_grid_clone_cfg, sht_grid_synthesis, sht_grid_sph_synthesis, sht_grid_analysis, sht_grid_sph_analysis
   use, intrinsic :: iso_c_binding, only: c_ptr
   !$ use omp_lib
   implicit none
   private

   public :: tensor_sh
   public :: tensor_sh_init, tensor_sh_synth, tensor_sh_analysis, tensor_sh_thread_cfg, tensor_sh_destroy
   integer, parameter, public :: TLAM = 4          ! λ = 1,2,5,6 → local 1..4
   ! Dyadic-field plane indices (the third dimension of the dyad array).
   integer, parameter, public :: DY_RR = 1, DY_RT = 2, DY_RP = 3, &
                                 DY_TT = 4, DY_TP = 5, DY_PP = 6

   type :: tensor_sh
      integer :: lmax = 0, nlm = 0, nphi = 0, nlat = 0
      integer,  allocatable :: ldeg(:)   !! (nlm) degree l of each coefficient
      integer,  allocatable :: mord(:)   !! (nlm) order m of each coefficient
      real(wp), allocatable :: llp1(:)   !! (nlm) l(l+1)
      real(wp), allocatable :: cott(:)   !! (nlat) cotθ at the Gauss latitudes
      real(wp), allocatable :: invsin(:) !! (nlat) 1/sinθ
      real(wp), allocatable :: n6(:)     !! (0:lmax) spin-2 channel norm S₆*S₆ (calibrated)
      ! Per-thread SHTns config pool: one config per OpenMP thread, so the element
      ! loop in the memory advance can run the dyadic transforms concurrently (a
      ! single config is NOT safe for concurrent calls). Built serially at init.
      type(c_ptr), allocatable :: pool(:)   !! (npool) independent configs
      integer :: npool = 0
   end type tensor_sh

contains

   subroutine tensor_sh_init(self, sht)
      type(tensor_sh), intent(out) :: self
      type(sht_grid),   intent(in)  :: sht
      complex(wp), allocatable :: c6(:), craw(:)
      real(wp),    allocatable :: tt(:,:), pp(:,:), tp(:,:)
      integer :: l, m, lm, i

      self%lmax = sht%lmax;  self%nlm = sht%nlm
      self%nphi = sht%nphi;  self%nlat = sht%nlat
      allocate(self%ldeg(self%nlm), self%mord(self%nlm), self%llp1(self%nlm))
      do m = 0, sht%mmax*sht%mres, sht%mres
         do l = m, sht%lmax
            lm = sht_grid_lmidx(sht, l, m)
            self%ldeg(lm) = l;  self%mord(lm) = m
            self%llp1(lm) = real(l,wp)*real(l+1,wp)
         end do
      end do
      allocate(self%cott(self%nlat), self%invsin(self%nlat))
      do i = 1, self%nlat
         self%cott(i)   = cos(sht%colat(i))/sin(sht%colat(i))
         self%invsin(i) = 1.0_wp/sin(sht%colat(i))
      end do

      ! Calibrate the spin-2 channel norm per degree: send a unit Z⁶ coefficient
      ! through synthesis then its adjoint; the diagonal response is n6(l). Use a
      ! representative order m (sectoral, m=l) so both G and H are present (l≥2).
      allocate(self%n6(0:self%lmax));  self%n6 = 0.0_wp
      allocate(c6(self%nlm), craw(self%nlm))
      allocate(tt(self%nphi,self%nlat), pp(self%nphi,self%nlat), tp(self%nphi,self%nlat))
      do l = 2, self%lmax
         m = min(l, sht%mmax*sht%mres)
         if (m < 1) cycle                     ! need m≥1 for H (else stays axisymmetric)
         lm = sht_grid_lmidx(sht, l, m)
         c6 = (0.0_wp, 0.0_wp);  c6(lm) = (1.0_wp, 0.0_wp)
         call spin2_synth(self, sht, c6, tt, pp, tp)      ! Sg→±(tt,pp), 4Sh→tp
         call spin2_adjoint(self, sht, tt, pp, tp, craw)  ! raw S₆* (unnormalised)
         self%n6(l) = real(craw(lm), wp)
      end do
      ! Axisymmetric fallback: for runs with mmax=0 the sectoral calibration above is
      ! skipped (m<1); calibrate at m=0 (H≡0, only G contributes) so 1-D still works.
      if (sht%mmax == 0) then
         do l = 2, self%lmax
            lm = sht_grid_lmidx(sht, l, 0)
            c6 = (0.0_wp, 0.0_wp);  c6(lm) = (1.0_wp, 0.0_wp)
            call spin2_synth(self, sht, c6, tt, pp, tp)
            call spin2_adjoint(self, sht, tt, pp, tp, craw)
            self%n6(l) = real(craw(lm), wp)
         end do
      end if
      deallocate(c6, craw, tt, pp, tp)

      ! Per-thread config pool (serial creation — FFTW planning is not thread-safe).
      self%npool = 1
      !$ self%npool = omp_get_max_threads()
      allocate(self%pool(self%npool))
      do i = 1, self%npool
         self%pool(i) = sht_grid_clone_cfg(sht)
      end do
   end subroutine tensor_sh_init

   function tensor_sh_thread_cfg(self) result(cfg)
      !! The calling OpenMP thread's private SHTns config (1-based pool index =
      !! thread id + 1). Serial / non-OpenMP builds always get pool(1).
      type(tensor_sh), intent(in) :: self
      type(c_ptr) :: cfg
      integer :: tid
      tid = 0
      !$ tid = omp_get_thread_num()
      cfg = self%pool(tid+1)
   end function tensor_sh_thread_cfg

   ! --- synthesis -------------------------------------------------------------

   subroutine tensor_sh_synth(self, sht, c, dyad, cfg)
      !! Tensor-harmonic coefficients c(λ=1..4, nlm) → six dyadic grid fields.
      !! Pass `cfg` (a thread_cfg handle) to transform on a thread-local config.
      type(tensor_sh), intent(in)  :: self
      type(sht_grid),   intent(in)  :: sht
      complex(wp),      intent(in)  :: c(:,:)        !! (TLAM, nlm)
      real(wp),         intent(out) :: dyad(:,:,:)   !! (nphi, nlat, 6)
      type(c_ptr), intent(in), optional :: cfg
      complex(wp) :: scaled(self%nlm)
      real(wp)    :: tr(self%nphi,self%nlat), sg(self%nphi,self%nlat), sh(self%nphi,self%nlat)
      real(wp)    :: tt(self%nphi,self%nlat), pp(self%nphi,self%nlat), tp(self%nphi,self%nlat)
      ! rr (Z¹) and rθ,rφ (Z²)
      call sht_grid_synthesis(sht, c(1,:), dyad(:,:,DY_RR), cfg)
      call sht_grid_sph_synthesis(sht, c(2,:), dyad(:,:,DY_RT), dyad(:,:,DY_RP), cfg)
      ! trace from Z⁵:  −l(l+1) Y
      scaled = -self%llp1*c(3,:)
      call sht_grid_synthesis(sht, scaled, tr, cfg)
      ! spin-2 from Z⁶
      call spin2_synth(self, sht, c(4,:), tt, pp, tp, cfg)   ! tt=Sg, pp=−Sg, tp=4Sh
      sg = tt
      dyad(:,:,DY_TT) = tr + sg
      dyad(:,:,DY_PP) = tr - sg
      dyad(:,:,DY_TP) = tp
   end subroutine tensor_sh_synth

   subroutine spin2_synth(self, sht, c6, tt, pp, tp, cfg)
      !! Z⁶ contribution: tt=Sg=Σc6·G, pp=−Sg, tp=4Sh=4Σc6·H, via the exact grid
      !! identities (no recurrence, no re-analysis).
      type(tensor_sh), intent(in)  :: self
      type(sht_grid),   intent(in)  :: sht
      complex(wp),      intent(in)  :: c6(:)
      real(wp),         intent(out) :: tt(:,:), pp(:,:), tp(:,:)
      type(c_ptr), intent(in), optional :: cfg
      complex(wp) :: imc(self%nlm)
      real(wp) :: f(self%nphi,self%nlat), gt(self%nphi,self%nlat), gp(self%nphi,self%nlat)
      real(wp) :: gtf(self%nphi,self%nlat), gpf(self%nphi,self%nlat)
      real(wp) :: lap(self%nphi,self%nlat), sg(self%nphi,self%nlat), sh(self%nphi,self%nlat)
      imc = cmplx(0.0_wp, real(self%mord,wp), wp)*c6            ! im·c6  (= ∂_φ on coeffs)
      call sht_grid_synthesis(sht, c6, f, cfg)
      call sht_grid_sph_synthesis(sht, c6, gt, gp, cfg)                   ! g_θ, g_φ
      call sht_grid_sph_synthesis(sht, imc, gtf, gpf, cfg)                ! ∂_φ g_θ, ∂_φ g_φ
      call sht_grid_synthesis(sht, cmplx(-self%llp1,0.0_wp,wp)*c6, lap, cfg)   ! ∇₁²f = −l(l+1)f
      ! Sg = ∇₁²f − 2cotθ g_θ − 2(1/sinθ)∂_φ g_φ ;  Sh = (1/sinθ)∂_φ g_θ − cotθ g_φ
      sg = lap - 2.0_wp*byprof(gt, self%cott) - 2.0_wp*byprof(gpf, self%invsin)
      sh =        byprof(gtf, self%invsin)     -        byprof(gp,  self%cott)
      tt = sg;  pp = -sg;  tp = 4.0_wp*sh
   end subroutine spin2_synth

   ! --- analysis --------------------------------------------------------------

   subroutine tensor_sh_analysis(self, sht, dyad, c, cfg)
      !! Six dyadic grid fields → tensor-harmonic coefficients.
      !! Pass `cfg` (a thread_cfg handle) to transform on a thread-local config.
      type(tensor_sh), intent(in)    :: self
      type(sht_grid),   intent(in)    :: sht
      real(wp),         intent(inout) :: dyad(:,:,:)   !! (nphi,nlat,6); SHTns overwrites
      complex(wp),      intent(out)   :: c(:,:)        !! (TLAM, nlm)
      type(c_ptr), intent(in), optional :: cfg
      complex(wp) :: craw(self%nlm)
      real(wp)    :: vt(self%nphi,self%nlat), vp(self%nphi,self%nlat)
      integer     :: lm
      ! rr (Z¹): scalar analysis (inverse of synth)
      call sht_grid_analysis(sht, dyad(:,:,DY_RR), c(1,:), cfg)
      ! rθ,rφ (Z²): spheroidal vector analysis (inverse of sph_synth)
      vt = dyad(:,:,DY_RT);  vp = dyad(:,:,DY_RP)
      call sht_grid_sph_analysis(sht, vt, vp, c(2,:), cfg)
      ! trace (Z⁵): analysis(θθ+φφ) = −2 l(l+1) T⁵
      vt = dyad(:,:,DY_TT) + dyad(:,:,DY_PP)
      call sht_grid_analysis(sht, vt, craw, cfg)
      do lm = 1, self%nlm
         if (self%llp1(lm) > 0.0_wp) then
            c(3,lm) = craw(lm)/(-2.0_wp*self%llp1(lm))
         else
            c(3,lm) = (0.0_wp, 0.0_wp)
         end if
      end do
      ! spin-2 (Z⁶): adjoint of spin2_synth, normalised by the calibrated per-degree n6
      call spin2_adjoint(self, sht, dyad(:,:,DY_TT), dyad(:,:,DY_PP), dyad(:,:,DY_TP), craw, cfg)
      do lm = 1, self%nlm
         if (self%n6(self%ldeg(lm)) /= 0.0_wp) then
            c(4,lm) = craw(lm)/self%n6(self%ldeg(lm))
         else
            c(4,lm) = (0.0_wp, 0.0_wp)
         end if
      end do
   end subroutine tensor_sh_analysis

   subroutine spin2_adjoint(self, sht, dtt, dpp, dtp, craw, cfg)
      !! Unnormalised adjoint of spin2_synth: craw = S_g*(θθ−φφ) + S_h*(θφ), where the
      !! *-operators swap each forward op (synth↔analysis, grid-multiply self-adjoint,
      !! im → −im). The ∫dΩ adjoint of synth IS analysis (= Wᵀ·synth), so this returns
      !! ∫(Z⁶-basis):(reconstructed tensor) — the projection numerator. dtt/dpp/dtp are
      !! overwritten by the SHTns analyses.
      type(tensor_sh), intent(in)    :: self
      type(sht_grid),   intent(in)    :: sht
      real(wp),         intent(inout) :: dtt(:,:), dpp(:,:), dtp(:,:)
      complex(wp),      intent(out)   :: craw(:)
      type(c_ptr), intent(in), optional :: cfg
      complex(wp) :: q(self%nlm), s(self%nlm)
      real(wp)    :: D(self%nphi,self%nlat), vt(self%nphi,self%nlat), vp(self%nphi,self%nlat), z(self%nphi,self%nlat)
      ! The ∫dΩ adjoint of the spheroidal vector synth (SHsph_to_spat) is l(l+1)·
      ! spat_to_SHsphtor (its INVERSE differs from its adjoint by the spheroidal norm
      ! l(l+1)); the scalar synth's adjoint is plain analysis (orthonormal). So every
      ! vector-analysis result below is scaled by llp1 to be the true adjoint.
      D = dtt - dpp                                   ! θθ−φφ feeds S_g*
      z = 0.0_wp
      ! S_g*(D) = −l(l+1)·analysis(D) − 2·sphAnal(cotθ·D,0).S + 2 im·sphAnal(0,(1/sinθ)·D).S
      call sht_grid_analysis(sht, D, q, cfg);   craw = -self%llp1*q
      vt = byprof(D, self%cott);  vp = z
      call sht_grid_sph_analysis(sht, vt, vp, s, cfg);   craw = craw - 2.0_wp*self%llp1*s
      vt = z;  vp = byprof(D, self%invsin)
      call sht_grid_sph_analysis(sht, vt, vp, s, cfg)
      craw = craw + 2.0_wp*cmplx(0.0_wp, real(self%mord,wp), wp)*self%llp1*s
      ! θφ contributes ∫T:Z⁶|_θφ = a_θφ·4H·(e_θφ:e_θφ=½) = 2 a_θφ H ⇒ 2·S_h*(θφ), with
      ! S_h*(D) = −im·llp1·sphAnal((1/sinθ)D,0).S − llp1·sphAnal(0,cotθ·D).S.
      D = dtp
      vt = byprof(D, self%invsin);  vp = z
      call sht_grid_sph_analysis(sht, vt, vp, s, cfg)
      craw = craw - 2.0_wp*cmplx(0.0_wp, real(self%mord,wp), wp)*self%llp1*s
      vt = z;  vp = byprof(D, self%cott)
      call sht_grid_sph_analysis(sht, vt, vp, s, cfg)
      craw = craw - 2.0_wp*self%llp1*s
   end subroutine spin2_adjoint

   ! --- helpers ---------------------------------------------------------------

   pure function byprof(field, prof) result(out)
      !! Multiply each Gauss-latitude column of a (nphi,nlat) field by prof(nlat).
      real(wp), intent(in) :: field(:,:), prof(:)
      real(wp) :: out(size(field,1), size(field,2))
      integer  :: i
      do i = 1, size(field,2)
         out(:,i) = field(:,i)*prof(i)
      end do
   end function byprof

   subroutine tensor_sh_destroy(self)
      type(tensor_sh), intent(inout) :: self
      integer :: i
      if (allocated(self%ldeg))   deallocate(self%ldeg, self%mord, self%llp1)
      if (allocated(self%cott))   deallocate(self%cott, self%invsin)
      if (allocated(self%n6))     deallocate(self%n6)
      if (allocated(self%pool)) then
         do i = 1, self%npool;  call sht_free_cfg(self%pool(i));  end do
         deallocate(self%pool)
      end if
      self%npool = 0;  self%lmax = 0;  self%nlm = 0
   end subroutine tensor_sh_destroy

end module fe_tensor_sh
