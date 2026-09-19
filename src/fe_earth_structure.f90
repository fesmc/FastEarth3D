module fe_earth_structure
   !! Reference Earth structure: a radially layered, incompressible Maxwell Earth,
   !! plus an optional 3D (laterally varying) viscosity field.
   !!
   !! The layered description is the *physical* model (piecewise-constant density,
   !! shear modulus, and viscosity between interface radii) — e.g. the benchmark
   !! model M3-L70-V01. The radial finite-element mesh that the solver discretizes
   !! onto is built from this by fe_radial_fe.
   !!
   !! 3D-ready (project goal): `visc_3d`, when allocated, carries absolute
   !! log10-viscosity on the Gauss-Legendre spatial grid per radial node — how
   !! VILMA injects lateral heterogeneity (Albrecht et al. 2024). 1D runs leave it
   !! unallocated; the same solver path reduces to the spherically symmetric case.
   use fe_precision, only: wp
   use fe_constants, only: pi, grav_G
   use fe_params,    only: fe_param_class, MAX_LAYER
   use fe_sht,       only: sht_grid
   use ncio,         only: nc_read, nc_size
   implicit none
   private

   real(wp), parameter :: DEG2RAD_ = acos(-1.0_wp)/180.0_wp
   real(wp), parameter :: R_PREM   = 6371.0e3_wp   !! PREM reference surface radius [m]

   ! Rheology of a layer.
   integer, parameter, public :: RHEOL_ELASTIC = 0   !! elastic (viscosity -> inf)
   integer, parameter, public :: RHEOL_MAXWELL = 1   !! Maxwell viscoelastic
   integer, parameter, public :: RHEOL_FLUID   = 2   !! inviscid fluid (mu = 0)

   public :: earth_layer, earth_model, build_M3L70V01, build_earth, prem_rho_mu
   public :: fe_read_visc_3d, load_visc_3d
   public :: earth_n_layers, earth_rho_at, earth_mu_at, earth_eta_at, earth_mass_below, earth_total_mass, earth_gravity_at, earth_moi, earth_is_3d

   type :: earth_layer
      !! A homogeneous spherical shell, r_bot <= r <= r_top.
      real(wp) :: r_bot = 0.0_wp   !! inner radius [m]
      real(wp) :: r_top = 0.0_wp   !! outer radius [m]
      real(wp) :: rho   = 0.0_wp   !! density [kg m^-3]
      real(wp) :: mu    = 0.0_wp   !! shear modulus / rigidity [Pa]
      real(wp) :: eta   = 0.0_wp   !! viscosity [Pa s] (huge for elastic, 0 fluid)
      integer  :: rheology = RHEOL_MAXWELL
   end type earth_layer

   type :: earth_model
      character(len=:), allocatable :: name
      real(wp) :: r_earth = 0.0_wp        !! surface radius [m]
      real(wp) :: r_core  = 0.0_wp        !! core-mantle boundary radius [m]
      type(earth_layer), allocatable :: layers(:)   !! surface-first (index 1 = top)
      ! Optional lateral viscosity: ABSOLUTE log10(η [Pa·s]) on the Gauss grid ×
      ! FE radial nodes, (nphi*nlat, nr). Populated by fe_read_visc_3d from a real
      ! lon-lat-r field (rung 6c). The node→element bridge (ve_response) takes the
      ! log10-mean of the two bracketing nodes and forms the perturbation against
      ! the element's radial reference η — storing the absolute field here avoids a
      ! per-node reference ambiguity at layer interfaces.
      real(wp), allocatable :: visc_3d(:,:)
   end type earth_model

contains

   ! --- Construction ----------------------------------------------------------

   function build_earth(p) result(em)
      !! Build the earth model selected by p%earth: a named built-in, "custom"
      !! assembled verbatim from the per-layer arrays (surface-first), or "PREM"
      !! — the same layer geometry/viscosity arrays, but with density and shear
      !! modulus auto-filled from the incompressible-PREM profile (see build_layered).
      type(fe_param_class), intent(in) :: p
      type(earth_model) :: em

      select case (trim(p%earth))
      case ("M3-L70-V01")
         em = build_M3L70V01()
      case ("custom")
         em = build_layered(p, use_prem=.false.)
      case ("PREM")
         em = build_layered(p, use_prem=.true.)
      case default
         error stop 'build_earth: unknown earth model "'//trim(p%earth)//'"'
      end select
   end function build_earth

   function build_layered(p, use_prem) result(em)
      !! Assemble a layered model from the surface-first parameter arrays
      !! (n_layer, r_bot, r_top, eta, rheology). With use_prem=.false. ("custom")
      !! every property — including rho and mu — is taken verbatim from the arrays.
      !! With use_prem=.true. ("PREM") the per-layer density and shear modulus are
      !! instead the r^2-weighted (mass-consistent) shell average of the
      !! incompressible-PREM profile over each layer (prem_shell_avg), so the rho/mu
      !! arrays are ignored; viscosity and rheology still come from the eta/rheology
      !! arrays (PREM defines no viscosity — this mirrors VILMA, where PREM supplies
      !! the elastic+density structure and a separate profile/3-D field supplies η).
      !! Fluid layers keep mu = 0.
      type(fe_param_class), intent(in) :: p
      logical,              intent(in) :: use_prem
      type(earth_model) :: em
      integer  :: k, n
      real(wp) :: rho_k, mu_k

      n = p%n_layer
      if (n < 1 .or. n > MAX_LAYER) &
         error stop 'build_layered: n_layer out of range for a layered earth model'
      if (use_prem) then; em%name = "PREM"; else; em%name = "custom"; end if
      em%r_earth = p%r_earth
      em%r_core  = p%r_core
      allocate(em%layers(n))
      do k = 1, n
         rho_k = p%rho(k)
         mu_k  = p%mu(k)
         if (use_prem) then
            call prem_shell_avg(p%r_bot(k), p%r_top(k), rho_k, mu_k)
            if (p%rheology(k) == RHEOL_FLUID) mu_k = 0.0_wp
         end if
         em%layers(k) = earth_layer(p%r_bot(k), p%r_top(k), rho_k, &
                                    mu_k, p%eta(k), p%rheology(k))
      end do
   end function build_layered

   pure subroutine prem_rho_mu(r, rho, mu)
      !! Incompressible-PREM density and shear modulus at radius r [m]
      !! (Dziewonski & Anderson 1981), as realized by VILMA's mod_polyprem.f90 —
      !! coefficients copied verbatim for parity; the canonical D&A-1981 Table I
      !! differs only in a few digits (e.g. the 5600-5701 km region). Returns
      !! rho [kg m^-3] and the shear modulus mu = rho*Vs^2 [Pa]; the incompressible
      !! model needs no bulk modulus, so Vp is not evaluated. The normalized radius
      !! is x = r / R_PREM with R_PREM = 6371 km (PREM's reference surface), keeping
      !! the polynomials true to PREM regardless of the model's r_earth. The values
      !! returned in the fluid core (r <= 3480 km) follow VILMA's fixed mu; callers
      !! set mu = 0 for fluid layers, so only rho is used there.
      real(wp), intent(in)  :: r
      real(wp), intent(out) :: rho, mu
      real(wp) :: x, vs

      x = r / R_PREM
      if (r <= 1221.5e3_wp) then              ! inner core
         rho = 13088.5_wp - 8838.1_wp*x**2
         mu  = 7.0363e10_wp
      else if (r <= 3480.0e3_wp) then         ! outer core (fluid)
         rho = 12581.5_wp - 1263.8_wp*x - 3642.6_wp*x**2 - 5528.1_wp*x**3
         mu  = 7.0363e10_wp
      else if (r <= 3630.0e3_wp) then         ! lowermost mantle (D'')
         vs  = 6925.4_wp + 1467.2_wp*x - 2083.4_wp*x**2 + 978.3_wp*x**3
         rho = 7956.5_wp - 6476.1_wp*x + 5528.3_wp*x**2 - 3080.7_wp*x**3
         mu  = vs*vs*rho
      else if (r <= 5600.0e3_wp) then         ! lower mantle
         vs  = 11167.1_wp - 13781.8_wp*x + 17457.5_wp*x**2 - 9277.7_wp*x**3
         rho = 7956.5_wp - 6476.1_wp*x + 5528.3_wp*x**2 - 3080.7_wp*x**3
         mu  = vs*vs*rho
      else if (r <= 5701.0e3_wp) then         ! lower mantle, just above 5600 km
         vs  = 22345.9_wp - 17247.3_wp*x - 2083.4_wp*x**2 + 978.3_wp*x**3
         rho = 7956.5_wp - 6476.1_wp*x + 5528.3_wp*x**2 - 3080.7_wp*x**3
         mu  = vs*vs*rho
      else if (r <= 5771.0e3_wp) then         ! transition zone (660-600 km)
         vs  = 9983.9_wp - 4932.4_wp*x
         rho = 5319.7_wp - 1483.6_wp*x
         mu  = vs*vs*rho
      else if (r <= 5971.0e3_wp) then         ! transition zone (600-400 km)
         vs  = 22351.2_wp - 18585.6_wp*x
         rho = 11249.4_wp - 8029.8_wp*x
         mu  = vs*vs*rho
      else if (r <= 6151.0e3_wp) then         ! upper mantle (400-220 km)
         vs  = 8949.6_wp - 4459.7_wp*x
         rho = 7108.9_wp - 3804.5_wp*x
         mu  = vs*vs*rho
      else if (r <= 6346.6e3_wp) then         ! upper mantle/lid (220 km-Moho), isotropic Vs
         vs  = 2151.9_wp + 2348.1_wp*x
         rho = 2691.0_wp + 692.4_wp*x
         mu  = vs*vs*rho
      else if (r <= 6356.0e3_wp) then         ! lower crust
         vs  = 3900.0_wp
         rho = 2900.0_wp
         mu  = vs*vs*rho
      else                                    ! upper crust (PREM ocean replaced by crust)
         vs  = 3200.0_wp
         rho = 2600.0_wp
         mu  = vs*vs*rho
      end if
   end subroutine prem_rho_mu

   pure subroutine prem_shell_avg(r_bot, r_top, rho, mu)
      !! r^2-weighted (volume) average of the incompressible-PREM density and shear
      !! modulus over the spherical shell [r_bot, r_top]: rho = layer mass / layer
      !! volume = ∫ rho(r) r^2 dr / ∫ r^2 dr (and likewise mu). This is the
      !! mass-consistent way to collapse the continuous PREM profile onto a finite
      !! layer (the SELEN/TABOO/ALMA3 convention) — at coarse layer counts it avoids
      !! the few-percent mass/gravity bias of midpoint sampling. Evaluated by
      !! composite Simpson quadrature; the profile is piecewise-polynomial, so with
      !! layer boundaries on the PREM discontinuities this is effectively exact.
      real(wp), intent(in)  :: r_bot, r_top
      real(wp), intent(out) :: rho, mu
      integer, parameter :: n = 200      ! Simpson sub-intervals (even)
      integer  :: i
      real(wp) :: h, r, w, ww, rho_i, mu_i, num_rho, num_mu, den

      h = (r_top - r_bot) / real(n, wp)
      num_rho = 0.0_wp;  num_mu = 0.0_wp;  den = 0.0_wp
      do i = 0, n
         if (i == 0 .or. i == n) then
            w = 1.0_wp
         else if (mod(i,2) == 1) then
            w = 4.0_wp
         else
            w = 2.0_wp
         end if
         r  = r_bot + real(i, wp)*h
         ww = w*r*r
         call prem_rho_mu(r, rho_i, mu_i)
         num_rho = num_rho + ww*rho_i
         num_mu  = num_mu  + ww*mu_i
         den     = den     + ww
      end do
      rho = num_rho / den
      mu  = num_mu  / den
   end subroutine prem_shell_avg

   function build_M3L70V01() result(em)
      !! Benchmark Earth model M3-L70-V01 (Spada et al. 2011, GJI 185:106,
      !! Table 3; values cross-checked against ALMA3 MODELS/M3L70V01.dat).
      !! Incompressible, self-gravitating, Maxwell. 70 km elastic lithosphere,
      !! 3 Maxwell mantle layers, inviscid fluid core.
      type(earth_model) :: em
      real(wp), parameter :: km = 1.0e3_wp

      em%name    = "M3-L70-V01"
      em%r_earth = 6371.0_wp*km
      em%r_core  = 3480.0_wp*km
      allocate(em%layers(5))
      !                       r_bot[m]      r_top[m]      rho      mu[Pa]      eta[Pa s]      rheology
      em%layers(1) = earth_layer(6301.0_wp*km, 6371.0_wp*km, 3037.0_wp,  5.0605e10_wp, huge(1.0_wp), RHEOL_ELASTIC) ! lithosphere
      em%layers(2) = earth_layer(5951.0_wp*km, 6301.0_wp*km, 3438.0_wp,  7.0363e10_wp, 1.0e21_wp,    RHEOL_MAXWELL) ! upper mantle
      em%layers(3) = earth_layer(5701.0_wp*km, 5951.0_wp*km, 3871.0_wp,  1.0549e11_wp, 1.0e21_wp,    RHEOL_MAXWELL) ! transition
      em%layers(4) = earth_layer(3480.0_wp*km, 5701.0_wp*km, 4978.0_wp,  2.2834e11_wp, 2.0e21_wp,    RHEOL_MAXWELL) ! lower mantle
      em%layers(5) = earth_layer(   0.0_wp,    3480.0_wp*km, 10750.0_wp, 0.0_wp,       0.0_wp,       RHEOL_FLUID)   ! core
   end function build_M3L70V01

   ! --- Queries ---------------------------------------------------------------

   integer function earth_n_layers(self) result(n)
      type(earth_model), intent(in) :: self
      n = 0
      if (allocated(self%layers)) n = size(self%layers)
   end function earth_n_layers

   integer function layer_at(self, r) result(k)
      !! Index of the layer containing radius r (-1 if outside the model).
      type(earth_model), intent(in) :: self
      real(wp),           intent(in) :: r
      integer :: i
      k = -1
      do i = 1, earth_n_layers(self)
         if (r >= self%layers(i)%r_bot .and. r <= self%layers(i)%r_top) then
            k = i
            return
         end if
      end do
   end function layer_at

   real(wp) function earth_rho_at(self, r) result(val)
      type(earth_model), intent(in) :: self
      real(wp),           intent(in) :: r
      integer :: k
      k = layer_at(self, r)
      val = 0.0_wp
      if (k > 0) val = self%layers(k)%rho
   end function earth_rho_at

   real(wp) function earth_mu_at(self, r) result(val)
      type(earth_model), intent(in) :: self
      real(wp),           intent(in) :: r
      integer :: k
      k = layer_at(self, r)
      val = 0.0_wp
      if (k > 0) val = self%layers(k)%mu
   end function earth_mu_at

   real(wp) function earth_eta_at(self, r) result(val)
      type(earth_model), intent(in) :: self
      real(wp),           intent(in) :: r
      integer :: k
      k = layer_at(self, r)
      val = 0.0_wp
      if (k > 0) val = self%layers(k)%eta
   end function earth_eta_at

   real(wp) function earth_mass_below(self, r) result(m)
      !! Mass enclosed within radius r [kg], integrating the piecewise-constant
      !! density profile (partial shells handled).
      type(earth_model), intent(in) :: self
      real(wp),           intent(in) :: r
      real(wp) :: lo, hi
      integer  :: i
      m = 0.0_wp
      do i = 1, earth_n_layers(self)
         lo = self%layers(i)%r_bot
         hi = min(self%layers(i)%r_top, r)
         if (hi > lo) m = m + (4.0_wp/3.0_wp)*pi*self%layers(i)%rho*(hi**3 - lo**3)
      end do
   end function earth_mass_below

   real(wp) function earth_total_mass(self) result(m)
      type(earth_model), intent(in) :: self
      m = earth_mass_below(self, self%r_earth)
   end function earth_total_mass

   real(wp) function earth_gravity_at(self, r) result(g)
      !! Unperturbed gravity at radius r [m s^-2]: G * M(<r) / r^2.
      type(earth_model), intent(in) :: self
      real(wp),           intent(in) :: r
      g = 0.0_wp
      if (r > 0.0_wp) g = grav_G*earth_mass_below(self, r)/(r*r)
   end function earth_gravity_at

   real(wp) function earth_moi(self) result(inertia)
      !! Moment of inertia about a diameter [kg m^2] for the spherically
      !! symmetric profile: I = (8 pi / 15) * sum rho * (r_top^5 - r_bot^5).
      type(earth_model), intent(in) :: self
      integer :: i
      inertia = 0.0_wp
      do i = 1, earth_n_layers(self)
         inertia = inertia + (8.0_wp*pi/15.0_wp)*self%layers(i)%rho* &
                   (self%layers(i)%r_top**5 - self%layers(i)%r_bot**5)
      end do
   end function earth_moi

   logical function earth_is_3d(self) result(yes)
      type(earth_model), intent(in) :: self
      yes = allocated(self%visc_3d)
   end function earth_is_3d

   ! --- 3D viscosity loading (real lon-lat-r log10(eta) fields) ----------------

   subroutine load_visc_3d(p, sht, r_node, visc_node)
      !! Read p%visc_3d_file onto the Gauss grid x FE nodes (absolute log10 eta),
      !! optionally perturb by f_visc_sd standard deviations, and clamp to
      !! [visc_log10_min, visc_log10_max]. sigma is read from name_visc_sd if that
      !! variable is named, else taken RELATIVE to the field, f_visc_rel*log10(eta)
      !! (so the perturbation tracks the viscosity structure rather than a constant
      !! floor). The clamp also imposes the viscosity floor on the raw field.
      type(fe_param_class), intent(in)  :: p
      type(sht_grid),       intent(in)  :: sht
      real(wp),             intent(in)  :: r_node(:)
      real(wp), allocatable, intent(out) :: visc_node(:,:)
      real(wp), allocatable :: sd_node(:,:)
      call fe_read_visc_3d(p%visc_3d_file, sht, r_node, visc_node, varname=p%name_visc, &
           lonname=p%name_visc_lon, latname=p%name_visc_lat, rname=p%name_visc_r)
      if (p%f_visc_sd /= 0.0_wp) then
         if (len_trim(p%name_visc_sd) > 0) then
            call fe_read_visc_3d(p%visc_3d_file, sht, r_node, sd_node, varname=p%name_visc_sd, &
                 lonname=p%name_visc_lon, latname=p%name_visc_lat, rname=p%name_visc_r)
         else
            allocate(sd_node, source=abs(visc_node)*p%f_visc_rel)   ! relative sigma [log10 dex]
         end if
         visc_node = visc_node + p%f_visc_sd*sd_node
      end if
      visc_node = min(max(visc_node, p%visc_log10_min), p%visc_log10_max)
   end subroutine load_visc_3d

   subroutine fe_read_visc_3d(filename, sht, r_node, visc_node, varname, &
                              lonname, latname, rname)
      !! Read a real 3D viscosity field (lon, lat, radius) from netCDF and
      !! interpolate it onto the Gauss-Legendre grid × the FE radial nodes,
      !! returning ABSOLUTE log10(η [Pa·s]) as visc_node(nphi*nlat, nr). Horizontal
      !! interp is bilinear in log10(η) with periodic longitude; vertical interp is
      !! linear in radius (clamped to the source range), also on log10(η). The
      !! node→element bridge and the perturbation against the radial reference η
      !! happen later in ve_response (which owns the per-element rheology).
      !!
      !! The stored variable is assumed to already be log10(η). Coordinate names
      !! default to lon/lat/r; the radius coordinate must be metres, ascending.
      character(len=*), intent(in)  :: filename
      type(sht_grid),   intent(in)  :: sht
      real(wp),         intent(in)  :: r_node(:)              !! (nr) FE node radii [m], ascending
      real(wp), allocatable, intent(out) :: visc_node(:,:)    !! (nphi*nlat, nr) log10(η)
      character(len=*), intent(in), optional :: varname, lonname, latname, rname
      character(len=64) :: vnm, lnm, tnm, rnm
      real(wp), allocatable :: lon_s(:), lat_s(:), r_s(:), eta_s(:,:,:)
      integer,  allocatable :: il0(:), il1(:), jl0(:), jl1(:)
      real(wp), allocatable :: wl(:), wt(:)
      integer  :: nlon, nlat_s, nr_s, nr, nphi, nlat
      integer  :: i, j, k, k0, k1, sp
      real(wp) :: lon_t, lat_t, span, wr, v0, v1

      vnm = 'eta';  lnm = 'lon';  tnm = 'lat';  rnm = 'r'
      if (present(varname)) vnm = varname
      if (present(lonname)) lnm = lonname
      if (present(latname)) tnm = latname
      if (present(rname))   rnm = rname

      nlon   = nc_size(filename, trim(lnm))
      nlat_s = nc_size(filename, trim(tnm))
      nr_s   = nc_size(filename, trim(rnm))
      allocate(lon_s(nlon), lat_s(nlat_s), r_s(nr_s), eta_s(nlon, nlat_s, nr_s))
      call nc_read(filename, trim(lnm), lon_s)
      call nc_read(filename, trim(tnm), lat_s)
      call nc_read(filename, trim(rnm), r_s)
      call nc_read(filename, trim(vnm), eta_s)        ! eta_s(lon, lat, r) = log10(η)

      nphi = sht%nphi;  nlat = sht%nlat;  nr = size(r_node)
      allocate(visc_node(nphi*nlat, nr))

      ! Horizontal interpolation weights, computed once (reused over all radii).
      span = lon_s(nlon) - lon_s(1) + (lon_s(2) - lon_s(1))   ! ≈ 360°
      allocate(il0(nphi), il1(nphi), wl(nphi))
      do i = 1, nphi
         lon_t = sht%lon(i) / DEG2RAD_                        ! [0,360)
         call locate_periodic(lon_s, lon_t, span, il0(i), il1(i), wl(i))
      end do
      allocate(jl0(nlat), jl1(nlat), wt(nlat))
      do j = 1, nlat
         lat_t = 90.0_wp - sht%colat(j)/DEG2RAD_
         call locate_clamped(lat_s, lat_t, jl0(j), jl1(j), wt(j))
      end do

      ! Per node: vertical bracket (clamped), then bilinear blend of the two levels.
      do k = 1, nr
         call locate_clamped(r_s, r_node(k), k0, k1, wr)
         do j = 1, nlat
            do i = 1, nphi
               sp = i + (j-1)*nphi
               v0 = bilin(eta_s(:,:,k0), il0(i), il1(i), wl(i), jl0(j), jl1(j), wt(j))
               v1 = bilin(eta_s(:,:,k1), il0(i), il1(i), wl(i), jl0(j), jl1(j), wt(j))
               visc_node(sp, k) = (1.0_wp - wr)*v0 + wr*v1
            end do
         end do
      end do
   end subroutine fe_read_visc_3d

   pure real(wp) function bilin(f, i0, i1, wi, j0, j1, wj) result(v)
      !! Bilinear sample of f(lon,lat) given precomputed bracketing indices/weights.
      real(wp), intent(in) :: f(:,:), wi, wj
      integer,  intent(in) :: i0, i1, j0, j1
      v = (1.0_wp-wi)*(1.0_wp-wj)*f(i0,j0) + wi*(1.0_wp-wj)*f(i1,j0) &
        + (1.0_wp-wi)*wj         *f(i0,j1) + wi*wj         *f(i1,j1)
   end function bilin

   pure subroutine locate_clamped(x, xt, i0, i1, w)
      !! Bracket xt in the ascending array x, clamping to the ends (w in [0,1]).
      real(wp), intent(in)  :: x(:), xt
      integer,  intent(out) :: i0, i1
      real(wp), intent(out) :: w
      integer :: n, i
      n = size(x)
      if (xt <= x(1))  then; i0 = 1; i1 = 1; w = 0.0_wp; return; end if
      if (xt >= x(n))  then; i0 = n; i1 = n; w = 0.0_wp; return; end if
      i0 = 1
      do i = 1, n-1
         if (xt >= x(i) .and. xt <= x(i+1)) then; i0 = i; exit; end if
      end do
      i1 = i0 + 1
      w  = (xt - x(i0)) / (x(i1) - x(i0))
   end subroutine locate_clamped

   pure subroutine locate_periodic(x, xt0, span, i0, i1, w)
      !! Bracket xt0 in the ascending periodic array x (period `span`), wrapping
      !! between x(n) and x(1)+span. xt0 is normalised into [x(1), x(1)+span).
      real(wp), intent(in)  :: x(:), xt0, span
      integer,  intent(out) :: i0, i1
      real(wp), intent(out) :: w
      integer  :: n, i
      real(wp) :: xt
      n = size(x)
      xt = xt0
      do while (xt <  x(1));        xt = xt + span; end do
      do while (xt >= x(1) + span); xt = xt - span; end do
      if (xt >= x(n)) then          ! wrap interval [x(n), x(1)+span)
         i0 = n; i1 = 1
         w  = (xt - x(n)) / (x(1) + span - x(n))
         return
      end if
      i0 = 1
      do i = 1, n-1
         if (xt >= x(i) .and. xt < x(i+1)) then; i0 = i; exit; end if
      end do
      i1 = i0 + 1
      w  = (xt - x(i0)) / (x(i1) - x(i0))
   end subroutine locate_periodic

end module fe_earth_structure
