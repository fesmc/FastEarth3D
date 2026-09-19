program test_tidal
   !! Rung-5a validation: the TIDAL forcing path (response to an EXTERNAL degree-j
   !! potential, fe_radial_fe%tidal_rhs / tidal_love) — the mechanism the rotational
   !! feedback (centrifugal potential, fe_rotation) drives the solid Earth with.
   !!
   !! Pinned by the homogeneous incompressible self-gravitating sphere, which has
   !! closed-form tidal Love numbers (Munk & MacDonald 1960; Kelvin):
   !!
   !!   fluid (μ→0):   k^T → 3/(2(n−1)),   h^T → (2n+1)/(2(n−1))
   !!   rigid (μ→∞):   h^T, l^T, k^T → 0
   !!   degree-2:      k^T = (3/2)/(1+μ̃),  h^T = (5/2)/(1+μ̃),  μ̃ = 19μ/(2ρga)
   !!
   !! and the elastic benchmark model M3-L70-V01 tidal Love numbers are printed
   !! (the secular k^T_f = k_s of the Liouville rotational-feedback term, Spada
   !! et al. 2011 eq 11, comes from the fluid limit of this same path).
   use fe_precision,       only: wp
   use fe_earth_structure, only: earth_gravity_at, earth_model, earth_layer, build_M3L70V01, &
                                 RHEOL_ELASTIC
   use fe_radial_fe,       only: radial_operator_tidal_rhs, radial_operator_destroy, radial_operator_solve_vec, radial_operator_assemble, radial_mesh_build, radial_mesh, radial_operator, tidal_love, &
                                 radial_fe_finalize, idx_u, idx_v, idx_f
   implicit none

   real(wp), parameter :: km = 1.0e3_wp
   real(wp), parameter :: rho_h = 5511.0_wp, a_h = 6371.0_wp*km
   integer  :: j
   real(wp) :: h, l, k, ua, va, fa, rsd, hf, kf, g, mut, ke, he
   logical  :: ok
   type(radial_operator) :: op

   ok = .true.

   ! --- 1. fluid limit: k^T → 3/(2(n−1)), h^T → (2n+1)/(2(n−1)) ----------------
   write(*,'(a)') ' (1) homogeneous near-fluid sphere -> tidal fluid limits'
   write(*,'(a)') '      j      h        h_f=(2j+1)/2(j-1)     k       k_f=3/2(j-1)   resid'
   do j = 2, 6
      call solve_homog_tidal(1.0_wp, j, ua, va, fa, rsd)        ! mu ~ 0
      call love_for(j, ua, va, fa, h, l, k)
      hf = real(2*j+1, wp)/(2.0_wp*real(j-1, wp))
      kf = 3.0_wp/(2.0_wp*real(j-1, wp))
      write(*,'(i7,2f12.5,4x,2f12.5,es11.2)') j, h, hf, k, kf, rsd
      if (abs(h - hf) > 1.0e-3_wp) then
         write(*,'(a)') '      FAIL: h^T off the fluid limit'; ok = .false.
      end if
      if (abs(k - kf) > 1.0e-3_wp) then
         write(*,'(a)') '      FAIL: k^T off the fluid limit'; ok = .false.
      end if
      if (rsd > 1.0e-8_wp) then
         write(*,'(a)') '      FAIL: solver did not converge'; ok = .false.
      end if
   end do

   ! --- 2. rigid limit: h, l, k → 0 -------------------------------------------
   write(*,'(a)') ''
   write(*,'(a)') ' (2) near-rigid sphere (mu=1e18) -> h,l,k -> 0'
   write(*,'(a)') '      j       h           l           k'
   do j = 2, 6
      call solve_homog_tidal(1.0e18_wp, j, ua, va, fa, rsd)
      call love_for(j, ua, va, fa, h, l, k)
      write(*,'(i7,3es13.4)') j, h, l, k
      if (abs(h) > 1.0e-4_wp .or. abs(l) > 1.0e-4_wp .or. abs(k) > 1.0e-4_wp) then
         write(*,'(a)') '      FAIL: not approaching the rigid limit'; ok = .false.
      end if
   end do

   ! --- 3. degree-2 elastic closed form: k=(3/2)/(1+μ̃), h=(5/2)/(1+μ̃) ---------
   write(*,'(a)') ''
   write(*,'(a)') ' (3) homogeneous degree-2 elastic vs Kelvin closed form'
   write(*,'(a)') '      mu[Pa]       h          h_exact       k          k_exact'
   g = grav_homog()
   call check_elastic(1.0e10_wp, g)
   call check_elastic(1.0e11_wp, g)
   call check_elastic(5.0e11_wp, g)

   ! --- 4. elastic M3-L70-V01 tidal Love numbers ------------------------------
   write(*,'(a)') ''
   write(*,'(a)') ' (4) elastic M3-L70-V01 tidal Love numbers'
   write(*,'(a)') '      j       h           l           k'
   call elastic_M3()

   write(*,'(a)') ''
   if (ok) then
      write(*,'(a)') ' PASS: tidal forcing path reproduces the analytic homogeneous-'
      write(*,'(a)') '       sphere tidal Love-number limits (fluid, rigid, Kelvin)'
   else
      write(*,'(a)') ' FAIL: tidal Love-number validation did not all pass'
      call radial_fe_finalize()
      error stop 1
   end if
   call radial_fe_finalize()

contains

   subroutine solve_homog_tidal(mu, j, ua, va, fa, rsd)
      !! Single-layer incompressible sphere of rigidity mu, forced by a unit
      !! external degree-j potential (tidal_rhs), solved for the surface response.
      real(wp), intent(in)  :: mu
      integer,  intent(in)  :: j
      real(wp), intent(out) :: ua, va, fa, rsd
      type(earth_model)     :: e
      type(radial_mesh)     :: m
      real(wp), allocatable :: x(:)
      integer :: it
      e = homog_earth(mu)
      call radial_mesh_build(m, e)
      call radial_operator_assemble(op, e, m, j)
      allocate(x(op%ndof))
      call radial_operator_solve_vec(op, radial_operator_tidal_rhs(op, 1.0_wp), x, iters=it, resid=rsd)
      ua = x(idx_u(m%nr));  va = x(idx_v(m%nr));  fa = x(idx_f(m%nr))
      call radial_operator_destroy(op)
   end subroutine solve_homog_tidal

   function homog_earth(mu) result(e)
      real(wp), intent(in) :: mu
      type(earth_model) :: e
      e%name = "homog";  e%r_earth = a_h;  e%r_core = 0.0_wp
      allocate(e%layers(1))
      e%layers(1) = earth_layer(0.0_wp, a_h, rho_h, mu, huge(1.0_wp), RHEOL_ELASTIC)
   end function homog_earth

   real(wp) function grav_homog() result(g)
      type(earth_model) :: e
      e = homog_earth(1.0_wp)
      g = earth_gravity_at(e, a_h)
   end function grav_homog

   subroutine love_for(j, ua, va, fa, h, l, k)
      integer,  intent(in)  :: j
      real(wp), intent(in)  :: ua, va, fa
      real(wp), intent(out) :: h, l, k
      type(earth_model) :: e
      e = homog_earth(1.0_wp)
      call tidal_love(e, j, 1.0_wp, ua, va, fa, h, l, k)
   end subroutine love_for

   subroutine check_elastic(mu, g)
      !! Degree-2 elastic tidal Love numbers of the homogeneous sphere vs the
      !! Kelvin closed form k=(3/2)/(1+μ̃), h=(5/2)/(1+μ̃), μ̃=19μ/(2ρga).
      real(wp), intent(in) :: mu, g
      real(wp) :: hh, ll, kk, uu, vv, ff, rr
      call solve_homog_tidal(mu, 2, uu, vv, ff, rr)
      call love_for(2, uu, vv, ff, hh, ll, kk)
      mut = 19.0_wp*mu/(2.0_wp*rho_h*g*a_h)
      ke  = 1.5_wp/(1.0_wp + mut)
      he  = 2.5_wp/(1.0_wp + mut)
      write(*,'(es12.2,4f12.6)') mu, hh, he, kk, ke
      if (abs(hh - he) > 2.0e-3_wp .or. abs(kk - ke) > 2.0e-3_wp) then
         write(*,'(a)') '      FAIL: degree-2 elastic off the Kelvin closed form'
         ok = .false.
      end if
   end subroutine check_elastic

   subroutine elastic_M3()
      type(earth_model)     :: e
      type(radial_mesh)     :: m
      real(wp), allocatable :: x(:)
      integer  :: jj, it
      real(wp) :: uu, vv, ff, hh, ll, kk, rr
      e = build_M3L70V01()
      call radial_mesh_build(m, e)
      do jj = 2, 8
         call radial_operator_assemble(op, e, m, jj)
         allocate(x(op%ndof))
         call radial_operator_solve_vec(op, radial_operator_tidal_rhs(op, 1.0_wp), x, iters=it, resid=rr)
         uu = x(idx_u(m%nr));  vv = x(idx_v(m%nr));  ff = x(idx_f(m%nr))
         call tidal_love(e, jj, 1.0_wp, uu, vv, ff, hh, ll, kk)
         write(*,'(i7,3es13.5)') jj, hh, ll, kk
         deallocate(x)
         if (rr > 1.0e-6_wp) then
            write(*,'(a)') '      FAIL: M3 tidal solve did not converge'; ok = .false.
         end if
         if (hh /= hh .or. kk /= kk) then
            write(*,'(a)') '      FAIL: M3 tidal Love numbers not finite'; ok = .false.
         end if
         if (jj == 2 .and. (kk <= 0.0_wp .or. kk >= 0.6_wp)) then
            write(*,'(a)') '      FAIL: M3 elastic k^T_2 outside physical range'; ok = .false.
         end if
         call radial_operator_destroy(op)
      end do
   end subroutine elastic_M3

end program test_tidal
