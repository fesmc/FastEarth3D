module fe_params
   !! Physics + numerics configuration record for the solid-Earth model, loaded from
   !! one namelist group `&fe3d` (yelmo convention: a flat parameter type filled by
   !! nml_read, see fesm-utils/utils/src/nml.f90). This is the host API contract:
   !! the high-level system init (solid_earth_init) consumes the whole record and
   !! distributes the values to the sub-solvers, while the specific component inits
   !! keep their direct-argument signatures.
   !!
   !! A host model (CLIMBER-X) fills this record in memory and drives the model
   !! through the API (z_bed_eq/h_ice_eq passed to solid_earth_init, each interval to
   !! solid_earth_update). The standalone driver reads it from a file. Run-management
   !! settings the *executables* need — forcing/reference/output paths, the time
   !! window, online-remap and i_eq selectors, restart-in path — live in the separate
   !! `&ctl` group (fe_control), which a host never reads.
   use fe_precision, only: wp
   use fe_constants, only: kyr, sec_per_year
   use nml
   implicit none
   private

   public :: fe_param_class, fe_par_load, fe_par_print
   public :: MAX_LAYER

   integer, parameter :: MAX_LAYER = 16   !! cap on custom earth-structure layers

   type :: fe_param_class
      ! --- spectral grid (fe_sht) ------------------------------------------------
      integer  :: lmax  = 0       !! maximum spherical-harmonic degree (required)
      integer  :: nlat  = 0       !! Gauss latitudes  (0 => SHTns default = lmax+2)
      integer  :: nphi  = 0       !! longitudes       (0 => SHTns default)
      integer  :: mmax  = -1      !! maximum order    (<0 => = lmax)
      integer  :: mres  = 1       !! order stride
      real(wp) :: eps_polar = -1.0_wp   !! polar-optimization threshold (<0 => library default)

      ! --- earth structure (fe_earth_structure) ---------------------------------
      character(len=64) :: earth = "M3-L70-V01"
         !! named built-in model, or "custom" to build from the layer arrays below
      integer  :: n_layer = 0                  !! # custom layers (surface-first)
      real(wp) :: r_earth = 6371.0e3_wp        !! custom: surface radius [m]
      real(wp) :: r_core  = 3480.0e3_wp        !! custom: core-mantle boundary radius [m]
      real(wp) :: r_bot(MAX_LAYER) = 0.0_wp    !! custom: layer inner radii [m]
      real(wp) :: r_top(MAX_LAYER) = 0.0_wp    !! custom: layer outer radii [m]
      real(wp) :: rho(MAX_LAYER)   = 0.0_wp    !! custom: layer densities [kg m^-3]
      real(wp) :: mu(MAX_LAYER)    = 0.0_wp    !! custom: layer shear moduli [Pa]
      real(wp) :: eta(MAX_LAYER)   = 0.0_wp    !! custom: layer viscosities [Pa s]
      integer  :: rheology(MAX_LAYER) = 1      !! custom: 0=elastic 1=Maxwell 2=fluid

      ! --- viscoelastic memory scheme (fe_response / fe_viscoelastic) ------------
      character(len=8) :: scheme = "fe"        !! fe | etd1 | trap | be
      integer  :: max_couple_iter = 20         !! SLE<->memory co-convergence cap (implicit schemes)

      ! --- response kind selector (fe_response / fe_modal) -----------------------
      character(len=8)  :: earth_response = "ve"  !! ve | modal | elastic | null
      integer  :: n_modes   = -1                  !! modal: # modes/degree (<=0 => all above tol)
      character(len=12) :: mode_rank = "isostatic" !! modal rank metric: isostatic | rate | residue
      character(len=8)  :: lat_method = "coupled" !! modal lateral-η method: coupled | lie | strang
      real(wp) :: dt_be     = kyr                 !! modal eigensolve backward-Euler Δt [s] (nml in YEARS)
      integer  :: n_krylov  = 20                  !! modal: Arnoldi/Krylov block size (caps modes/degree found)
      logical  :: modal_adaptive = .false.        !! modal: sub-step each coupling interval to rtol (A3).
                                                  !! .false. (default) = 1 exact step/interval (fast). Only
                                                  !! worth it for 1-D temporal accuracy; does NOT help 3-D.

      ! --- sea-level equation (fe_sle) ------------------------------------------
      integer  :: sle_n_outer      = 3
      integer  :: sle_n_inner      = 20
      real(wp) :: sle_tol          = 1.0e-7_wp
      integer  :: sle_max_mem_iter = 20
      logical  :: sle_fixed_ocean  = .false.
      logical  :: sle_subgrid      = .true.

      ! --- adaptive time stepping (fe_timestep) ----------------------------------
      ! The Δt fields are SI [s] in the record; the nml supplies them in YEARS and
      ! fe_par_load converts on read (so the in-memory record is uniformly SI). The
      ! coupling cadence is NOT a parameter: the host (or the standalone driver's
      ! forcing axis) passes each interval to solid_earth_update.
      real(wp) :: dt_init   = 0.0_wp         !! first trial Δt [s] (0 => try the whole interval)
      real(wp) :: dt_min    = 0.0_wp         !! Δt floor [s] (0 => none)
      real(wp) :: dt_max    = huge(1.0_wp)   !! Δt ceiling [s]
      real(wp) :: rtol      = 1.0e-4_wp      !! relative local-error tolerance (memory ∞-norm)
      real(wp) :: atol      = 1.0e-3_wp      !! absolute local-error floor
      real(wp) :: safety    = 0.9_wp         !! step-size safety factor
      real(wp) :: grow_max  = 5.0_wp         !! max Δt growth per accepted step
      real(wp) :: shrink_min = 0.2_wp        !! min Δt shrink per step
      real(wp) :: cfl       = 1.0_wp         !! explicit (fe) sub-step Maxwell-number ceiling M=μΔt/η

      ! --- rotational feedback (fe_rotation) ------------------------------------
      logical  :: rotation = .true.          !! TPW feedback (on for real runs; off for non-rotating benchmarks)

      ! --- LGM-memory spin-up (fe_coupling solid_earth_spinup) -------------------
      ! A model capability the host opts into. Relax under the start-slice ice while
      ! HOLDING the reference (z_bed_eq, h_ice_eq) as the datum, so the transient
      ! enters with viscous memory. The standalone driver triggers it; a host
      ! (CLIMBER-X) leaves equil_time_max=0 / pre_spinup_1d=.false. to skip it and
      ! takes whatever state it passes in as the (zero-memory) equilibrium.
      real(wp) :: equil_time_max = 0.0_wp
         !! >0: spin up the LGM memory state by holding the start-slice ice (relaxing to
         !! isostatic equilibrium) in the FULL model BEFORE the transient, with the
         !! reference held as the datum. A cap [years]: the relaxation exits early once
         !! the bed stops moving, or at this time with a warning if it has not
         !! converged. =0 skips the full-model equilibration phase. Non-default.
      real(wp) :: equil_rate_tol = 1.0e-3_wp
         !! spin-up convergence criterion: the mean bed velocity over a relaxation pass
         !! [m/yr]. Relaxation exits once the rate drops below this.
      logical :: pre_spinup_1d = .false.
         !! .true.: before the full-model phase, run a cheap 1-D pre-equilibration to
         !! bed-stationary convergence (the 1-D radial viscosity is the lateral
         !! geometric mean of the 3-D field), then switch to the full model carrying the
         !! spun-up memory. Independent of equil_time_max.

      ! --- 3D viscosity field + uncertainty sampling (fe_earth_structure) --------
      ! Mirrors the CLIMBER-X VILMA scheme (src/geo/vilma.F90) but with a RELATIVE
      ! 1-sigma instead of a constant floor: perturb log10(eta) by f_visc_sd*sigma,
      ! sigma read from the file if name_visc_sd is set, else f_visc_rel*log10(eta).
      logical  :: l_visc_3d   = .false.   !! load a lateral log10(eta) field
      character(len=512) :: visc_3d_file  = ""       !! lon-lat-r log10(eta) field
      character(len=64)  :: name_visc     = "eta"    !! viscosity var (log10 Pa s)
      character(len=64)  :: name_visc_lon = "lon"
      character(len=64)  :: name_visc_lat = "lat"
      character(len=64)  :: name_visc_r   = "r"
      character(len=64)  :: name_visc_sd  = ""       !! optional sigma(log10 eta) var; "" => relative
      real(wp) :: f_visc_sd      = 0.0_wp   !! perturbation in units of sigma (0 = mean field)
      real(wp) :: f_visc_rel     = 0.1_wp   !! relative sigma = f_visc_rel*log10(eta) when no sd var
      real(wp) :: visc_log10_min = 19.5_wp  !! floor on log10(eta) after read + perturbation [dex]
      real(wp) :: visc_log10_max = 30.0_wp  !! ceiling on log10(eta) after read + perturbation [dex]
      real(wp) :: visc3d_tol     = 1.0e-3_wp !! lateral log10(eta) spread [dex] above which a radial
         !! element is treated as genuinely 3-D (pays the dyadic SHT round-trip); below it the
         !! element collapses to its lateral-mean scalar rate (cheap degree-diagonal path). Raising
         !! it demotes weakly-3-D elements to 1-D and cuts the memory-advance cost (the dominant cost).
   end type fe_param_class

contains

   subroutine fe_par_load(p, filename, defaults_file, group)
      !! Fill the whole parameter record from the `&fe3d` group of `filename`,
      !! overlaid on a complete `defaults_file` (yelmo convention): every parameter
      !! must exist in the defaults file, but the user `filename` may set only the
      !! subset it wants to override. If `defaults_file` is omitted, `filename` IS
      !! its own defaults — so it must then be complete. Override `group` to read a
      !! differently-named namelist.
      type(fe_param_class), intent(inout) :: p
      character(len=*),     intent(in)    :: filename
      character(len=*),     intent(in), optional :: defaults_file
      character(len=*),     intent(in), optional :: group
      character(len=64)  :: g
      character(len=512) :: df
      real(wp) :: dt_init_yr, dt_min_yr, dt_max_yr
      real(wp) :: equil_time_max_yr, dt_be_yr

      g  = "fe3d";      if (present(group))         g  = group
      df = filename;    if (present(defaults_file)) df = defaults_file
      call nml_set_verbose(.false.)             ! fe_par_print echoes a concise summary instead

      ! grid
      call nml_read(filename, g, "lmax",      p%lmax,      defaults_file=df)
      call nml_read(filename, g, "nlat",      p%nlat,      defaults_file=df)
      call nml_read(filename, g, "nphi",      p%nphi,      defaults_file=df)
      call nml_read(filename, g, "mmax",      p%mmax,      defaults_file=df)
      call nml_read(filename, g, "mres",      p%mres,      defaults_file=df)
      call nml_read(filename, g, "eps_polar", p%eps_polar, defaults_file=df)

      ! earth structure
      call nml_read(filename, g, "earth",     p%earth,     defaults_file=df)
      call nml_read(filename, g, "n_layer",   p%n_layer,   defaults_file=df)
      call nml_read(filename, g, "r_earth",   p%r_earth,   defaults_file=df)
      call nml_read(filename, g, "r_core",    p%r_core,    defaults_file=df)
      call nml_read(filename, g, "r_bot",     p%r_bot,     defaults_file=df)
      call nml_read(filename, g, "r_top",     p%r_top,     defaults_file=df)
      call nml_read(filename, g, "rho",       p%rho,       defaults_file=df)
      call nml_read(filename, g, "mu",        p%mu,        defaults_file=df)
      call nml_read(filename, g, "eta",       p%eta,       defaults_file=df)
      call nml_read(filename, g, "rheology",  p%rheology,  defaults_file=df)

      ! memory scheme
      call nml_read(filename, g, "scheme",          p%scheme,          defaults_file=df)
      call nml_read(filename, g, "max_couple_iter", p%max_couple_iter, defaults_file=df)

      ! response kind selector (dt_be given in YEARS, converted to SI below)
      call nml_read(filename, g, "earth_response",  p%earth_response,  defaults_file=df)
      call nml_read(filename, g, "n_modes",         p%n_modes,         defaults_file=df)
      call nml_read(filename, g, "mode_rank",       p%mode_rank,       defaults_file=df)
      call nml_read(filename, g, "lat_method",      p%lat_method,      defaults_file=df)
      call nml_read(filename, g, "n_krylov",        p%n_krylov,        defaults_file=df)
      call nml_read(filename, g, "modal_adaptive",  p%modal_adaptive,  defaults_file=df)
      dt_be_yr = p%dt_be/sec_per_year
      call nml_read(filename, g, "dt_be",           dt_be_yr,          defaults_file=df)
      p%dt_be = dt_be_yr*sec_per_year

      ! sea-level equation
      call nml_read(filename, g, "sle_n_outer",      p%sle_n_outer,      defaults_file=df)
      call nml_read(filename, g, "sle_n_inner",      p%sle_n_inner,      defaults_file=df)
      call nml_read(filename, g, "sle_tol",          p%sle_tol,          defaults_file=df)
      call nml_read(filename, g, "sle_max_mem_iter", p%sle_max_mem_iter, defaults_file=df)
      call nml_read(filename, g, "sle_fixed_ocean",  p%sle_fixed_ocean,  defaults_file=df)
      call nml_read(filename, g, "sle_subgrid",      p%sle_subgrid,      defaults_file=df)

      ! adaptive time stepping. The Δt fields are given in YEARS in the nml and
      ! converted to SI seconds here (the record is uniformly SI internally).
      dt_init_yr = p%dt_init/sec_per_year
      dt_min_yr  = p%dt_min /sec_per_year
      dt_max_yr  = p%dt_max /sec_per_year
      call nml_read(filename, g, "dt_init",    dt_init_yr,   defaults_file=df)
      call nml_read(filename, g, "dt_min",     dt_min_yr,    defaults_file=df)
      call nml_read(filename, g, "dt_max",     dt_max_yr,    defaults_file=df)
      p%dt_init = dt_init_yr*sec_per_year
      p%dt_min  = dt_min_yr *sec_per_year
      p%dt_max  = dt_max_yr *sec_per_year
      call nml_read(filename, g, "rtol",       p%rtol,       defaults_file=df)
      call nml_read(filename, g, "atol",       p%atol,       defaults_file=df)
      call nml_read(filename, g, "safety",     p%safety,     defaults_file=df)
      call nml_read(filename, g, "grow_max",   p%grow_max,   defaults_file=df)
      call nml_read(filename, g, "shrink_min", p%shrink_min, defaults_file=df)
      call nml_read(filename, g, "cfl",        p%cfl,        defaults_file=df)

      ! rotation
      call nml_read(filename, g, "rotation",   p%rotation,   defaults_file=df)

      ! LGM-memory spin-up (equil_time_max given in YEARS, converted to SI below)
      equil_time_max_yr = p%equil_time_max/sec_per_year
      call nml_read(filename, g, "equil_time_max", equil_time_max_yr, defaults_file=df)
      p%equil_time_max = equil_time_max_yr*sec_per_year
      call nml_read(filename, g, "equil_rate_tol", p%equil_rate_tol, defaults_file=df)
      call nml_read(filename, g, "pre_spinup_1d",  p%pre_spinup_1d,  defaults_file=df)

      ! 3D viscosity + uncertainty
      call nml_read(filename, g, "l_visc_3d",      p%l_visc_3d,      defaults_file=df)
      call nml_read(filename, g, "visc_3d_file",   p%visc_3d_file,   defaults_file=df)
      call nml_read(filename, g, "visc3d_tol",     p%visc3d_tol,     defaults_file=df)
      call nml_read(filename, g, "name_visc",      p%name_visc,      defaults_file=df)
      call nml_read(filename, g, "name_visc_lon",  p%name_visc_lon,  defaults_file=df)
      call nml_read(filename, g, "name_visc_lat",  p%name_visc_lat,  defaults_file=df)
      call nml_read(filename, g, "name_visc_r",    p%name_visc_r,    defaults_file=df)
      call nml_read(filename, g, "name_visc_sd",   p%name_visc_sd,   defaults_file=df)
      call nml_read(filename, g, "f_visc_sd",      p%f_visc_sd,      defaults_file=df)
      call nml_read(filename, g, "f_visc_rel",     p%f_visc_rel,     defaults_file=df)
      call nml_read(filename, g, "visc_log10_min", p%visc_log10_min, defaults_file=df)
      call nml_read(filename, g, "visc_log10_max", p%visc_log10_max, defaults_file=df)
   end subroutine fe_par_load

   subroutine fe_par_print(p, unit)
      !! Echo the active configuration (to stdout, or `unit` if given).
      type(fe_param_class), intent(in) :: p
      integer, intent(in), optional :: unit
      integer :: u, k
      u = 6;  if (present(unit)) u = unit

      write(u,'(a)')          ' [fe3d] configuration'
      write(u,'(a,i0,a,i0,a,i0)') '   grid:   lmax=', p%lmax, '  nlat=', p%nlat, '  nphi=', p%nphi
      write(u,'(a,a)')        '   earth:  ', trim(p%earth)
      if (trim(p%earth) == "custom") then
         do k = 1, p%n_layer
            write(u,'(a,i0,a,es9.2,a,es9.2,a,f8.1,a,es9.2,a,es9.2,a,i0)') &
                 '     layer ', k, ': r=[', p%r_bot(k), ',', p%r_top(k), &
                 ']  rho=', p%rho(k), '  mu=', p%mu(k), '  eta=', p%eta(k), &
                 '  rheol=', p%rheology(k)
         end do
      end if
      write(u,'(a,a)')        '   response: ', trim(p%earth_response)
      if (trim(p%earth_response) == "modal") then
         write(u,'(a,i0,a,a,a,es9.2,a,i0,a,l1)') '     modal: n_modes=', p%n_modes, &
              '  mode_rank=', trim(p%mode_rank), '  dt_be=', p%dt_be, &
              '  n_krylov=', p%n_krylov, '  adaptive=', p%modal_adaptive
         write(u,'(a,a)')     '     modal: lat_method=', trim(p%lat_method)
      end if
      write(u,'(a,a,a,i0)')   '   scheme: ', trim(p%scheme), '   max_couple_iter=', p%max_couple_iter
      write(u,'(a,i0,a,i0,a,es8.1,a,l1,a,l1)') &
           '   sle:    n_outer=', p%sle_n_outer, '  n_inner=', p%sle_n_inner, &
           '  tol=', p%sle_tol, '  fixed_ocean=', p%sle_fixed_ocean, '  subgrid=', p%sle_subgrid
      write(u,'(a,es9.2,a,es8.1,a,es8.1,a,f5.2)') &
           '   dt:     init=', p%dt_init, &
           '  rtol=', p%rtol, '  atol=', p%atol, '  cfl=', p%cfl
      write(u,'(a,l1)')       '   rotation: ', p%rotation
      write(u,'(a,es9.2,a,es9.2,a,l1)') '   spinup: equil_time_max=', p%equil_time_max, &
           '  equil_rate_tol=', p%equil_rate_tol, '  pre_spinup_1d=', p%pre_spinup_1d
      if (p%l_visc_3d) then
         write(u,'(a,a)')   '   visc_3d: ', trim(p%visc_3d_file)
         write(u,'(a,f6.2,a,f6.2,a,f5.2,a,f5.2,a)') &
              '            f_visc_sd=', p%f_visc_sd, '  f_visc_rel=', p%f_visc_rel, &
              '  clamp=[', p%visc_log10_min, ',', p%visc_log10_max, ']'
         write(u,'(a,es9.2,a)') '            visc3d_tol=', p%visc3d_tol, ' dex (3-D split)'
      end if
   end subroutine fe_par_print

end module fe_params
