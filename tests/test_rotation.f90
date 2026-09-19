program test_rotation
   !! Rung-5b validation: rotational feedback / polar motion (fe_rotation) vs the
   !! Spada et al. (2011) Test 3/2 benchmark (data/benchmarks/rotation_spada2011).
   !!
   !! Three checks, increasingly physical:
   !!   (1) NORMALIZATION: the rigid (2,1) inertia of each load, I_rigid/(C−A),
   !!       reproduces the paper's published geometrical factor G (Table 14 caption
   !!       / Spada eq. 31). This isolates the spherical-harmonic → inertia mapping
   !!       from the viscoelastic / Liouville physics.
   !!   (2) SECULAR Love number k_s = k^T_f and elastic k^T_e are physical.
   !!   (3) POLAR MOTION |m(t)| under a Heaviside load matches Table 14 (the
   !!       Chandler-excluded, Cw=0, column — the quasi-static regime fe_rotation
   !!       integrates) at t = 0,1,2,5,10,20 kyr, for the cap AND disc loads.
   use fe_precision,       only: wp
   use fe_earth_structure, only: earth_model, build_M3L70V01
   use fe_radial_fe,       only: radial_fe_finalize
   use fe_rotation,        only: rotation_update, rotation_destroy, rotation_init, rotation_state
   use fe_sht,             only: sht_grid, sht_grid_init, sht_grid_destroy, sht_grid_surface_integral
   implicit none

   real(wp), parameter :: deg = acos(-1.0_wp)/180.0_wp
   real(wp), parameter :: rad2deg = 180.0_wp/acos(-1.0_wp)
   real(wp), parameter :: yr = 3.15576e7_wp                ! Spada Table 2 sec/yr
   real(wp), parameter :: kyr = 1.0e3_wp*yr
   real(wp), parameter :: rho_i = 931.0_wp, alpha = 10.0_wp*deg
   real(wp), parameter :: thetac = 25.0_wp*deg, lambdac = 75.0_wp*deg
   ! benchmark times and Table 14 |m| [deg], Cw=0 (Gs column)
   real(wp), parameter :: t_ref(6) = [0.0_wp,1.0_wp,2.0_wp,5.0_wp,10.0_wp,20.0_wp]
   real(wp), parameter :: m_cap(6) = [0.0132_wp,0.0161_wp,0.0180_wp,0.0211_wp,0.0240_wp,0.0292_wp]
   real(wp), parameter :: m_dsc(6) = [0.0131_wp,0.0160_wp,0.0179_wp,0.0210_wp,0.0239_wp,0.0290_wp]
   ! published geometrical factors (Table 14 caption)
   complex(wp), parameter :: G_cap = -cmplx(0.541e-4_wp, 0.202e-3_wp, wp)
   complex(wp), parameter :: G_dsc = -cmplx(0.539e-4_wp, 0.201e-3_wp, wp)

   type(earth_model)    :: earth
   type(sht_grid)       :: sht
   type(rotation_state) :: rot
   real(wp), allocatable :: load(:,:)
   real(wp) :: dt, k_s, kTe
   logical  :: ok

   ok = .true.
   earth = build_M3L70V01()
   call sht_grid_init(sht, 128, nlat=256, nphi=512)
   allocate(load(sht%nphi, sht%nlat))
   dt = 25.0_wp*yr

   ! Love numbers are load-independent: report once.
   call rotation_init(rot, earth, sht, dt)
   k_s = rot%k_s;  kTe = rot%kTe
   write(*,'(a)')       ' (2) tidal Love numbers (degree 2)'
   write(*,'(a,f9.4)')  '      elastic  k^T_e     = ', kTe
   write(*,'(a,f9.4)')  '      secular  k_s=k^T_f = ', k_s
   if (k_s <= 0.5_wp .or. k_s >= 1.5_wp) then
      write(*,'(a)') '      FAIL: secular k_s outside the physical range'; ok = .false.
   end if
   if (kTe <= 0.0_wp .or. kTe >= k_s) then
      write(*,'(a)') '      FAIL: elastic k^T_e unphysical (expect 0 < k^T_e < k_s)'; ok = .false.
   end if
   call rotation_destroy(rot)

   call run_load('cap',  1.5e3_wp, G_cap, m_cap)
   call run_load('disc', 1.0e3_wp, G_dsc, m_dsc)

   write(*,'(a)') ''
   if (ok) then
      write(*,'(a)') ' PASS: rotational feedback reproduces Spada 2011 Test 3/2 polar'
      write(*,'(a)') '       motion (rigid inertia, secular k_s, |m(t)|) for cap + disc'
   else
      write(*,'(a)') ' FAIL: rotation validation did not all pass'
      call sht_grid_destroy(sht);  call radial_fe_finalize();  error stop 1
   end if
   call sht_grid_destroy(sht);  call radial_fe_finalize()

contains

   subroutine run_load(name, h, Gref, mref)
      character(len=*), intent(in) :: name
      real(wp),         intent(in) :: h, mref(6)
      complex(wp),      intent(in) :: Gref
      complex(wp) :: Gnum
      real(wp) :: t_now, mdeg, relG, relm
      integer  :: istep, nsteps, iref

      call build_load(name, h, load)

      ! (1) normalization vs published G
      Gnum = inertia21(sht, load, earth%r_earth)/2.63e35_wp
      relG = abs(Gnum - Gref)/abs(Gref)
      write(*,'(a)') ''
      write(*,'(3a)') ' (1) ', name, ': rigid (2,1) inertia vs published G (eq. 31)'
      write(*,'(a,2es12.4)') '      G ref   = ', real(Gref), aimag(Gref)
      write(*,'(a,2es12.4)') '      G model = ', real(Gnum), aimag(Gnum)
      write(*,'(a,f8.3,a)')  '      rel.err = ', 100.0_wp*relG, ' %'
      if (relG > 0.02_wp) then
         write(*,'(3a)') '      FAIL: ', name, ' rigid inertia off published G (>2%)'; ok = .false.
      end if

      ! (3) polar motion |m(t)| vs Table 14 (Cw=0)
      call rotation_init(rot, earth, sht, dt)
      rot%enabled = .true.
      write(*,'(3a)') ' (3) ', name, ': polar motion |m(t)| vs Table 14 (Cw=0)'
      write(*,'(a)')  '      t[kyr]   |m| model[deg]   |m| ref[deg]    rel.err'
      nsteps = nint(20.0_wp*kyr/dt)
      do istep = 0, nsteps
         t_now = rot%time
         call rotation_update(rot, sht, load, dt)            ! rot%m is now m(t_now)
         iref = ref_index(t_now/kyr)
         if (iref > 0) then
            mdeg = abs(rot%m)*rad2deg
            relm = abs(mdeg - mref(iref))/max(mref(iref), tiny(1.0_wp))
            write(*,'(f8.1,2f16.5,f12.3)') t_now/kyr, mdeg, mref(iref), relm
            if (relm > 0.03_wp) then
               write(*,'(3a)') '      FAIL: ', name, ' |m| off the benchmark (>3%)'; ok = .false.
            end if
         end if
      end do
      call rotation_destroy(rot)
   end subroutine run_load

   subroutine build_load(name, h, load)
      !! Cap (parabolic, Table 4): σ = ρ_i h √[(cosγ−cosα)/(1−cosα)]; disc (uniform):
      !! σ = ρ_i h, both for angular distance γ ≤ α from the centroid (θc,λc).
      character(len=*), intent(in)  :: name
      real(wp),         intent(in)  :: h
      real(wp),         intent(out) :: load(:,:)
      real(wp) :: ca, cg, frac
      integer  :: il, ip
      ca = cos(alpha);  load = 0.0_wp
      do il = 1, sht%nlat
         do ip = 1, sht%nphi
            cg = cos(thetac)*cos(sht%colat(il)) &
               + sin(thetac)*sin(sht%colat(il))*cos(sht%lon(ip) - lambdac)
            if (cg >= ca) then
               if (name == 'cap') then
                  frac = (cg - ca)/(1.0_wp - ca)
                  load(ip,il) = rho_i*h*sqrt(max(frac, 0.0_wp))
               else
                  load(ip,il) = rho_i*h
               end if
            end if
         end do
      end do
   end subroutine build_load

   complex(wp) function inertia21(sht, load, a) result(I21)
      !! (2,1) inertia I₁₃+iI₂₃ = −a⁴∫σ sinθcosθ e^{iφ}dΩ (the check mirrors the
      !! integral fe_rotation uses internally).
      type(sht_grid), intent(in) :: sht
      real(wp),       intent(in) :: load(:,:), a
      real(wp), allocatable :: w13(:,:), w23(:,:)
      integer  :: il, ip
      allocate(w13(sht%nphi,sht%nlat), w23(sht%nphi,sht%nlat))
      do il = 1, sht%nlat
         do ip = 1, sht%nphi
            w13(ip,il) = load(ip,il)*sin(sht%colat(il))*cos(sht%colat(il))*cos(sht%lon(ip))
            w23(ip,il) = load(ip,il)*sin(sht%colat(il))*cos(sht%colat(il))*sin(sht%lon(ip))
         end do
      end do
      I21 = cmplx(-a**4*sht_grid_surface_integral(sht, w13), -a**4*sht_grid_surface_integral(sht, w23), wp)
   end function inertia21

   integer function ref_index(tk) result(idx)
      real(wp), intent(in) :: tk
      integer :: i
      idx = 0
      do i = 1, size(t_ref)
         if (abs(tk - t_ref(i)) < 0.4_wp*dt/kyr) then;  idx = i;  return;  end if
      end do
   end function ref_index

end program test_rotation
