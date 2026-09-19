program test_benchmark_martinec
   !! Rung-4 validation against the Martinec et al. (2018) SLE benchmark
   !! (data/benchmarks/sle_martinec2018/, see data/benchmarks/PROVENANCE.md).
   !!
   !! CASE A (this test): a spherical ice cap (alpha=10 deg, h0=1500 m, rho=931)
   !! at the pole, applied as a Heaviside at 10 kyr, with NO topography/ocean
   !! (B0=0). It is therefore a pure viscoelastic LOADING response 10 kyr after
   !! application -- no sea-level equation -- and the cleanest Martinec-2018 case.
   !! giapy drops degrees n=0,1 (jmin=2), so the degree-1 frame does not enter.
   !! The reference A_fig10_SBK.dat gives, along colatitude: col1 colat[deg],
   !! col2 uplift u[m], col3 horizontal th-displacement[m], col4 geoid N[m].
   !!
   !! This is the first full-field check of the viscoelastic response synthesized
   !! over many degrees, and the first benchmark of the HORIZONTAL displacement
   !! (hence an independent confirmation of the l Love number's sign/scale).
   !!
   !! The remaining Martinec-2018 cases (C2/D3/E2/F1) exercise the full migrating-
   !! coastline SLE with basin topography and time-evolving caps; they are a
   !! separate, larger effort (full-sphere SLE + arbitrary lon/lat profiles).
   use fe_precision,       only: wp
   use fe_constants,       only: pi, rho_ice, kyr
   use fe_earth_structure, only: earth_gravity_at, earth_model, build_M3L70V01
   use fe_radial_fe,       only: radial_mesh_build, radial_mesh, radial_fe_finalize
   use fe_viscoelastic,    only: ve_destroy, ve_step, ve_init, ve_degree
   implicit none
   character(*), parameter :: REF = 'data/benchmarks/sle_martinec2018/A_fig10_SBK.dat'
   integer,  parameter :: NROW = 721, NMAX = 128, NQ = 8000, NSTEP = 500
   real(wp), parameter :: ALPHA_DEG = 10.0_wp, H0 = 1500.0_wp
   real(wp), parameter :: TOL_U = 1.5e-2_wp, TOL_H = 1.5e-2_wp, TOL_N = 3.0e-2_wp
   real(wp) :: Un(NMAX), Vn(NMAX), Nn(NMAX), sig(NMAX)
   real(wp) :: alpha, ca, g, dt, t1, ua, va, fa
   real(wp) :: colat(NROW), uref(NROW), href(NROW), gref(NROW)
   real(wp) :: eu, eh, en
   logical  :: ok
   integer  :: n, istep
   type(earth_model) :: em
   type(radial_mesh) :: mm
   type(ve_degree)   :: ve

   ok = .true.
   call read_ref()
   em = build_M3L70V01();  call radial_mesh_build(mm, em)
   g  = earth_gravity_at(em, em%r_earth)
   alpha = ALPHA_DEG*pi/180.0_wp;  ca = cos(alpha)
   call cap_coeffs()                       ! spherical-cap load coefficients sig_n

   ! viscoelastic relaxation: held cap load for 10 kyr (dt=20 yr, 500 steps)
   dt = 0.02_wp*kyr
   Un = 0.0_wp;  Vn = 0.0_wp;  Nn = 0.0_wp
   do n = 2, NMAX                          ! n=0,1 dropped (giapy jmin=2)
      call ve_init(ve, em, mm, n, dt)
      do istep = 1, NSTEP
         call ve_step(ve, 1.0_wp, t1, ua, va, fa)
      end do
      Un(n) = ua;  Vn(n) = va;  Nn(n) = -fa/g
      call ve_destroy(ve)
   end do

   call synth_and_compare(eu, eh, en)
   write(*,'(a)') ' (A) cap loading response after 10 kyr vs Martinec-2018 fig.10'
   write(*,'(a,f9.3,a,f6.2,a)') '      uplift     peak ', maxval(abs(uref)), ' m, max err ', 100*eu, ' %'
   write(*,'(a,f9.3,a,f6.2,a)') '      horizontal peak ', maxval(abs(href)), ' m, max err ', 100*eh, ' %'
   write(*,'(a,f9.3,a,f6.2,a)') '      geoid      peak ', maxval(abs(gref)), ' m, max err ', 100*en, ' %'
   if (eu > TOL_U) then; write(*,'(a)') '      FAIL: uplift off the benchmark';     ok = .false.; end if
   if (eh > TOL_H) then; write(*,'(a)') '      FAIL: horizontal off the benchmark'; ok = .false.; end if
   if (en > TOL_N) then; write(*,'(a)') '      FAIL: geoid off the benchmark';      ok = .false.; end if

   write(*,'(a)') ''
   if (ok) then
      write(*,'(a)') ' PASS: Martinec-2018 case A (cap loading u/horizontal/geoid) matches'
   else
      write(*,'(a)') ' FAIL: Martinec-2018 case A validation did not all pass'
      call radial_fe_finalize();  error stop 1
   end if
   call radial_fe_finalize()

contains

   subroutine synth_and_compare(eu, eh, en)
      !! Synthesize u(theta)=Sum sig_n U_n P_n, horizontal=Sum sig_n V_n dP_n/dth,
      !! N=Sum sig_n N_n P_n along the reference colatitudes; peak-normalized error.
      real(wp), intent(out) :: eu, eh, en
      real(wp) :: th, x, st, us, vs, ns, pn, pnm1, pnp1, dpn, du, dh, dg
      integer  :: i
      du = 0.0_wp;  dh = 0.0_wp;  dg = 0.0_wp
      do i = 1, NROW
         th = colat(i)*pi/180.0_wp;  x = cos(th);  st = sin(th)
         pnm1 = 1.0_wp;  pn = x                ! P_0, P_1
         us = 0.0_wp;  vs = 0.0_wp;  ns = 0.0_wp
         do n = 1, NMAX
            ! pn = P_n, pnm1 = P_{n-1}; dP_n/dtheta = n(x P_n - P_{n-1})/sin th
            if (st > 1.0e-12_wp) then;  dpn = real(n,wp)*(x*pn - pnm1)/st;  else;  dpn = 0.0_wp;  end if
            if (n >= 2) then
               us = us + sig(n)*Un(n)*pn
               vs = vs + sig(n)*Vn(n)*dpn
               ns = ns + sig(n)*Nn(n)*pn
            end if
            pnp1 = (real(2*n+1,wp)*x*pn - real(n,wp)*pnm1)/real(n+1,wp)
            pnm1 = pn;  pn = pnp1
         end do
         du = max(du, abs(us-uref(i)));  dh = max(dh, abs(vs-href(i)));  dg = max(dg, abs(ns-gref(i)))
      end do
      eu = du/maxval(abs(uref));  eh = dh/maxval(abs(href));  en = dg/maxval(abs(gref))
   end subroutine synth_and_compare

   subroutine cap_coeffs()
      !! sig_n = (2n+1)/2 \int_0^alpha rho h0 sqrt((cos d-ca)/(1-ca)) P_n(cos d) sin d dd
      !! by midpoint quadrature in colatitude d (the sqrt cap profile is smooth
      !! except for a sqrt cusp at the edge; NQ=8000 resolves it to <1e-4).
      integer  :: iq
      real(wp) :: d, dd, w, s, cd, pnm1, pn, pnp1
      sig = 0.0_wp
      dd  = alpha/real(NQ,wp)
      do iq = 1, NQ
         d  = (real(iq,wp)-0.5_wp)*dd
         cd = cos(d)
         s  = rho_ice*H0*sqrt(max((cd-ca)/(1.0_wp-ca), 0.0_wp))
         w  = s*sin(d)*dd
         pnm1 = 1.0_wp;  pn = cd                 ! P_0, P_1
         do n = 1, NMAX
            sig(n) = sig(n) + w*pn               ! pn = P_n(cos d)
            pnp1 = (real(2*n+1,wp)*cd*pn - real(n,wp)*pnm1)/real(n+1,wp)
            pnm1 = pn;  pn = pnp1
         end do
      end do
      do n = 1, NMAX
         sig(n) = sig(n)*real(2*n+1,wp)*0.5_wp
      end do
   end subroutine cap_coeffs

   subroutine read_ref()
      integer :: u, i, ios
      open(newunit=u, file=REF, status='old', action='read', iostat=ios)
      if (ios /= 0) then;  write(*,'(2a)') ' FAIL: cannot read ', REF;  error stop 1;  end if
      read(u,*)                                  ! header comment line
      do i = 1, NROW
         read(u,*,iostat=ios) colat(i), uref(i), href(i), gref(i)
         if (ios /= 0) then;  write(*,'(a,i0)') ' FAIL: short read at row ', i;  error stop 1;  end if
      end do
      close(u)
   end subroutine read_ref

end program test_benchmark_martinec
