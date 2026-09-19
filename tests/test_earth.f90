program test_earth
   !! Validate the built-in Earth-structure datasets for internal consistency.
   !!
   !! M3-L70-V01: the unperturbed gravity computed from the density profile,
   !! g(r) = G*M(<r)/r^2, must reproduce the gravity values Spada et al. (2011)
   !! Table 3 states at each interface — tying the hard-coded densities and radii
   !! to published numbers. Mass and the moment-of-inertia factor are checked for
   !! physical plausibility.
   !!
   !! PREM: the 11-layer incompressible-PREM model (rho/mu auto-filled by the
   !! r^2-weighted shell average, build_earth(earth="PREM")) must reproduce the
   !! real Earth's total mass and MOI factor to <1% and surface gravity to ~9.8 —
   !! this is the mass-consistency check the volume-average buys over midpoint
   !! sampling. The pointwise prem_rho_mu evaluator is checked at the surface crust.
   use fe_precision,       only: wp
   use fe_constants,       only: m_earth
   use fe_params,          only: fe_param_class
   use fe_earth_structure, only: earth_moi, earth_total_mass, earth_gravity_at, earth_n_layers, &
                                 earth_model, build_M3L70V01, build_earth, prem_rho_mu
   implicit none

   type(earth_model)    :: em
   type(fe_param_class) :: p
   integer  :: i
   logical  :: ok
   real(wp) :: mass, moi_factor, g, gtol, rho_c, mu_c
   ! Benchmark interface radii [m] and stated gravities [m s^-2] (Spada 2011 Tab 3)
   real(wp), parameter :: r_if(5) = [6371.0_wp, 6301.0_wp, 5951.0_wp, &
                                     5701.0_wp, 3480.0_wp]*1.0e3_wp
   real(wp), parameter :: g_if(5) = [9.815_wp, 9.854_wp, 9.978_wp, &
                                     10.024_wp, 10.457_wp]

   ok = .true.
   em = build_M3L70V01()
   print '(a,a,a,i0,a)', ' model: ', em%name, '  (', earth_n_layers(em), ' layers)'

   ! --- Gravity profile vs benchmark interface values -------------------------
   ! 0.02 m/s^2 tolerance absorbs the benchmark's 4-sig-fig rounding and any
   ! difference in the gravitational constant convention.
   gtol = 0.02_wp
   print '(a)', '   interface       r [km]      g model     g bench      err'
   do i = 1, 5
      g = earth_gravity_at(em, r_if(i))
      print '(i10,4f12.4)', i, r_if(i)/1.0e3_wp, g, g_if(i), abs(g - g_if(i))
      if (abs(g - g_if(i)) >= gtol) ok = .false.
   end do

   ! --- Total mass: plausible vs real Earth (simplified model, so ~few %) -----
   mass = earth_total_mass(em)
   print '(a,es12.5,a,f7.4)', ' total mass = ', mass, ' kg   (M_earth ratio = ', &
        mass/m_earth, ')'
   if (mass/m_earth < 0.95_wp .or. mass/m_earth > 1.05_wp) ok = .false.

   ! --- Moment-of-inertia factor I/(M R^2): physically ~0.30-0.40 -------------
   moi_factor = earth_moi(em)/(mass*em%r_earth**2)
   print '(a,f8.5,a)', ' MOI factor I/(M R^2) = ', moi_factor, '  (Earth ~0.3307)'
   if (moi_factor < 0.30_wp .or. moi_factor > 0.40_wp) ok = .false.

   ! === PREM: 11-layer incompressible model (volume-averaged rho/mu) ==========
   p%earth   = "PREM"
   p%n_layer = 11
   p%r_earth = 6371.0e3_wp
   p%r_core  = 3480.0e3_wp
   p%r_bot(1:11) = [6291.0_wp, 6151.0_wp, 5971.0_wp, 5771.0_wp, 5701.0_wp, 5600.0_wp, &
                    4943.0_wp, 4287.0_wp, 3630.0_wp, 3480.0_wp,    0.0_wp]*1.0e3_wp
   p%r_top(1:11) = [6371.0_wp, 6291.0_wp, 6151.0_wp, 5971.0_wp, 5771.0_wp, 5701.0_wp, &
                    5600.0_wp, 4943.0_wp, 4287.0_wp, 3630.0_wp, 3480.0_wp]*1.0e3_wp
   p%eta(1:11)   = [1.0e40_wp, 4.0e20_wp, 4.0e20_wp, 4.0e20_wp, 4.0e20_wp, 1.0e22_wp, &
                    1.0e22_wp, 1.0e22_wp, 1.0e22_wp, 1.0e22_wp,    0.0_wp]
   p%rheology(1:11) = [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2]

   em = build_earth(p)
   print '(/,a,a,a,i0,a)', ' model: ', em%name, '  (', earth_n_layers(em), ' layers)'

   ! Total mass and MOI to <1% of the real Earth — the payoff of volume-averaging.
   mass = earth_total_mass(em)
   print '(a,es12.5,a,f7.4)', ' total mass = ', mass, ' kg   (M_earth ratio = ', &
        mass/m_earth, ')'
   if (abs(mass/m_earth - 1.0_wp) > 0.01_wp) ok = .false.

   moi_factor = earth_moi(em)/(mass*em%r_earth**2)
   print '(a,f8.5,a)', ' MOI factor I/(M R^2) = ', moi_factor, '  (Earth ~0.3307)'
   if (moi_factor < 0.32_wp .or. moi_factor > 0.34_wp) ok = .false.

   ! Surface gravity within 0.05 m/s^2 of the real ~9.81 (vs 10.10 for midpoint).
   g = earth_gravity_at(em, em%r_earth)
   print '(a,f8.4,a)', ' g(surface) = ', g, ' m/s^2  (Earth ~9.81)'
   if (abs(g - 9.81_wp) > 0.05_wp) ok = .false.

   ! Pointwise evaluator at the upper crust: rho=2600, mu=rho*Vs^2 with Vs=3200.
   call prem_rho_mu(6365.0e3_wp, rho_c, mu_c)
   print '(a,f8.1,a,es11.4,a)', ' prem_rho_mu(crust): rho = ', rho_c, '  mu = ', mu_c, ' Pa'
   if (abs(rho_c - 2600.0_wp) > 1.0e-6_wp .or. &
       abs(mu_c - 3200.0_wp**2*2600.0_wp) > 1.0_wp) ok = .false.

   if (ok) then
      print '(/,a)', ' PASS: M3-L70-V01 and PREM Earth-structure consistency checks'
   else
      print '(/,a)', ' FAIL: Earth-structure consistency check exceeded tolerance'
      error stop 1
   end if
end program test_earth
