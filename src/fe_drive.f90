module fe_drive
   !! Standalone forced-run driver. Given a &fe3d physics configuration (fe_params)
   !! and a &ctl run-control configuration (fe_control), build the transform grid and
   !! the solid-Earth model, read an ice-thickness forcing time series (and a bedrock
   !! reference) from netCDF, optionally remap it onto the Gauss grid on the fly, set
   !! up the reference / equilibration state, and march the model across the forcing,
   !! writing the diagnostic surface fields each step. All file/var/window/i_eq/remap
   !! knobs come from the &ctl record (c); the physics record (p) supplies the grid and
   !! the model. A host model (CLIMBER-X) replaces this driver with its own loop.
   !!
   !! General forced-transient driver: fields are read per (file, var, time-index) so
   !! experiments can mix and match. Two input modes:
   !!   remap_input=.true.  (default): the forcing is on a regular lon-lat grid and is
   !!                       conservatively remapped onto the model Gauss grid per slice
   !!                       (fe_remap); raw datasets stay on disk, no preprocessing.
   !!   remap_input=.false.: the forcing is already on the Gauss grid (host remapped,
   !!                       or produced by the offline fastearth_remap tool).
   !!
   !! Reference state (i_eq), mirroring CLIMBER-X i_equilibrium. RSL is measured
   !! against the relaxed reference z_bed_eq, so i_eq=1 yields rsl ~0 at present day.
   !!   i_eq=0  start slice is the equilibrium (z_bed_eq=bed[k0], h_ice_eq=ice[k0]).
   !!   i_eq=1  (default) present-day reference: z_bed_eq / h_ice_eq from
   !!           z_bed_ref_file / h_ice_ref_file (e.g. RTopo), remapped online.
   !!   i_eq=2  equilibrium read from z_bed_eq_file / h_ice_eq_file.
   !!   i_eq=3  z_bed_eq = present-day reference bed + rsl from rsl_restart_file.
   !! Optional (non-default) LGM-memory spin-up when pre_spinup_1d or equil_time_max>0:
   !! relax under the start-slice ice while HOLDING the i_eq reference as the datum, so
   !! the model enters the transient spun up (viscous memory) — see solid_earth_spinup.
   use fe_precision, only: wp
   use fe_constants, only: sec_per_year
   use fe_params,    only: fe_param_class, fe_par_load, fe_par_print
   use fe_control,   only: fe_ctl_class, fe_ctl_load, fe_ctl_print
   use fe_sht,       only: sht_grid, sht_grid_destroy, sht_grid_init
   use fe_coupling,  only: solid_earth_finalize, solid_earth_update, solid_earth_init, solid_earth, &
                           solid_earth_spinup
   use fe_response,  only: RESP_MODAL
   use fe_remap,     only: remap_ll_gauss, remap_init, remap_to_gauss
   use fe_io,        only: fe_write_step, fe_restart_write, fe_restart_read
   use ncio,         only: nc_read, nc_size, nc_exists_var
   use iso_fortran_env, only: error_unit
   implicit none
   private

   public :: fastearth_run

   ! diagnostic surface fields written each output step (the prognostic memory is
   ! written separately via fe_restart_write when a restart is wanted)
   character(len=8), parameter :: OUT_VARS(5) = &
        [character(len=8) :: "h_ice", "rsl", "z_bed", "C_ocean", "bsl"]

contains

   subroutine fastearth_run(cfg_file, defaults_file)
      !! Run a forced simulation defined by cfg_file. Its &fe3d (physics) group is
      !! overlaid on the complete physics defaults in defaults_file (the program
      !! passes fe_control's DEFAULTS_FILE = input/fastearth3d_defaults.nml); its &ctl
      !! (run control) group is read from cfg_file alone, with the fe_ctl_class in-code
      !! defaults filling any gaps.
      character(len=*), intent(in) :: cfg_file
      character(len=*), intent(in) :: defaults_file

      type(fe_param_class)   :: p
      type(fe_ctl_class)     :: c
      type(sht_grid), target :: sht
      type(solid_earth)      :: se
      type(remap_ll_gauss)     :: rmap
      real(wp), allocatable  :: z_bed_eq(:,:), h_ice_eq(:,:), h_ice(:,:), ice_lgm(:,:)
      real(wp), allocatable  :: src(:,:), tyr(:)
      real(wp) :: t0, dt
      integer  :: nt, k, k0, k1, np, nl, nlon, nls
      logical  :: remap
      character(len=:), allocatable :: rundir
      integer(kind=8) :: pc0, pc1, prate          ! PROFILE: per-step phase timers
      real(wp) :: t_read = 0.0_wp, t_upd = 0.0_wp, t_wrt = 0.0_wp
      real(wp) :: t_dr, t_mm                      ! PROFILE: solid_earth_update sub-phases
      integer  :: nstep = 0

      ! --- configuration --------------------------------------------------------
      ! &fe3d: cfg_file overlaid on the physics defaults. &ctl: from cfg_file alone
      ! (the defaults file carries no &ctl group — it is the host API contract).
      call fe_par_load(p, cfg_file, defaults_file=defaults_file)
      call fe_ctl_load(c, cfg_file)
      call fe_par_print(p)
      call fe_ctl_print(c)

      ! --- validate inputs up front (fail loudly on a config typo) --------------
      ! ncio's nc_read stop's with exit status 0 on a missing variable, so a typo
      ! in a name_* / file_* config would otherwise abort silently mid-run; check
      ! every (file, var) the run will read before doing any work.
      call validate_inputs(c)

      ! --- transform grid + work arrays -----------------------------------------
      call build_grid(p, sht)
      np = sht%nphi;  nl = sht%nlat
      allocate(z_bed_eq(np,nl), h_ice_eq(np,nl), h_ice(np,nl), ice_lgm(np,nl))

      ! --- build the conservative lon-lat -> Gauss map (online remap mode) -------
      remap = c%remap_input
      call system_clock(pc0, prate)
      if (remap) then
         call build_remap(c, sht, rmap, nlon, nls)
         allocate(src(nlon, nls))
         write(*,'(a,i0,a,i0,a,i0,a,i0,a)') ' remap: lon-lat (', nlon, 'x', nls, &
              ') -> Gauss (', np, 'x', nl, ')'
      end if
      call system_clock(pc1)
      write(*,'(a,f8.2,a)') ' [PROFILE setup] build_remap =', real(pc1-pc0,wp)/prate, ' s'

      ! --- forcing time axis (years) + window -----------------------------------
      nt = nc_size(c%file_forcing, trim(c%name_time))
      allocate(tyr(nt));  call nc_read(c%file_forcing, trim(c%name_time), tyr)
      call select_window(tyr*sec_per_year, c%time_init, c%time_end, k0, k1)
      if (k1 <= k0) error stop 'fastearth_run: forcing window contains < 2 time slices'
      write(*,'(a,i0,a,f0.1,a,f0.1,a)') ' fastearth: ', k1-k0+1, ' slices, t = ', &
           tyr(k0), ' -> ', tyr(k1), ' yr'

      ! --- start-slice ice + reference state (z_bed_eq, h_ice_eq) per i_eq ------
      call read_ice(c, rmap, sht, remap, k0, src, ice_lgm)      ! start-slice (LGM) ice
      call setup_reference(c, rmap, sht, remap, k0, src, ice_lgm, z_bed_eq, h_ice_eq)

      ! --- initialise; restart resume and/or optional LGM-memory spin-up --------
      ! Restarts are written alongside the run output: <rundir>/spinup after the spin-up
      ! and <rundir>/final at the end, where rundir is the directory of file_out.
      rundir = dir_of(c%file_out)
      call system_clock(pc0, prate)
      se%par = p; call solid_earth_init(se, z_bed_eq, h_ice_eq)              ! reference, memory 0
      if (len_trim(c%restart_in_file) > 0) &                                 ! resume saved memory + clock
         call fe_restart_read(se, trim(c%restart_in_file))
      if (p%pre_spinup_1d .or. p%equil_time_max > 0.0_wp) then
         ! LGM-memory spin-up: relax under the start-slice ice while HOLDING the
         ! reference (z_bed_eq, h_ice_eq) as the datum, so the transient measures
         ! rsl/bsl against the reference (-> 0 as ice -> h_ice_eq).
         call solid_earth_spinup(se, ice_lgm)
         call fe_restart_write(se, se%time, folder=trim(rundir)//"/spinup")
      else
         ! seed the entering (LGM) ice and re-solve the diagnostics (also fills
         ! rsl/z_bed/C after a cross-resolution restart); no integration.
         call solid_earth_update(se, ice_lgm, 0.0_wp)
      end if
      call system_clock(pc1)
      write(*,'(a,f8.2,a)') ' [PROFILE setup] solid_earth_init+visc3d+seed =', real(pc1-pc0,wp)/prate, ' s'

      ! Lateral-viscosity diagnostic. RESP_VE: how many radial elements are genuinely
      ! 3-D (pay the tensor-SH advance) vs collapse to the cheap 1-D path. RESP_MODAL:
      ! how many within-degree mode ranks carry a lateral anomaly (the K scalar-SHT path).
      if (p%l_visc_3d) then
         if (se%resp%kind == RESP_MODAL) then
            write(*,'(a,i0,a,i0,a)') ' modal lateral viscosity: ', se%resp%nrank3d, &
                 ' of ', se%resp%maxmode, ' mode ranks carry a lateral anomaly'
         else
            write(*,'(a,i0,a,i0,a)') ' visc3d split: ', se%resp%ne3d, &
                 ' of ', se%resp%ne, ' radial elements laterally 3-D (rest advance as 1-D)'
         end if
      end if

      ! --- march the transient --------------------------------------------------
      t0 = tyr(k0)*sec_per_year
      se%time = tyr(k0);  se%resp%time = t0          ! coupling clock in years; response clock in SI
      call fe_write_step(se, c%file_out, se%time, nms=OUT_VARS, init=.true.)

      se%resp%t_drift = 0.0_wp;  se%resp%t_mem = 0.0_wp   ! PROFILE: time the transient only
      se%resp%n_drift = 0;       se%resp%n_mem = 0
      do k = k0, k1-1
         dt = tyr(k+1) - tyr(k)                    ! coupling interval [years]
         call system_clock(pc0, prate)
         call read_ice(c, rmap, sht, remap, k+1, src, h_ice)
         call system_clock(pc1);  t_read = t_read + real(pc1-pc0,wp)/prate
         call system_clock(pc0)
         call solid_earth_update(se, h_ice, dt)
         call system_clock(pc1);  t_upd = t_upd + real(pc1-pc0,wp)/prate
         call system_clock(pc0)
         call fe_write_step(se, c%file_out, se%time, nms=OUT_VARS, init=.false.)
         call system_clock(pc1);  t_wrt = t_wrt + real(pc1-pc0,wp)/prate
         nstep = nstep + 1
         write(*,'(a,f12.2,a,es10.2,a,es9.2)') '   t=', se%time, &
              ' yr   max|rsl|=', maxval(abs(se%rsl)), '   mass_resid=', se%worst_mass_resid
      end do

      if (nstep > 0) write(*,'(a,i0,a,/,3(a,f8.1,a,f5.1,a,/))') &
         ' [PROFILE] per coupling step (mean over ', nstep, ' steps):', &
         '   read_ice (remap+IO) =', 1.0e3_wp*t_read/nstep, ' ms (', &
            100.0_wp*t_read/(t_read+t_upd+t_wrt), ' %)', &
         '   solid_earth_update  =', 1.0e3_wp*t_upd /nstep, ' ms (', &
            100.0_wp*t_upd /(t_read+t_upd+t_wrt), ' %)', &
         '   fe_write_step (out) =', 1.0e3_wp*t_wrt /nstep, ' ms (', &
            100.0_wp*t_wrt /(t_read+t_upd+t_wrt), ' %)'
      ! solid_earth_update internal split (wall-clock, accumulated in the response).
      ! drift = per-degree band LU; memory advance = the Maxwell update (3-D dyadic SHT
      ! round-trip when laterally 3-D); rest = SLE iteration + load/geoid SHTs.
      if (nstep > 0) then
         t_dr = se%resp%t_drift;  t_mm = se%resp%t_mem
         write(*,'(a,/,3(a,f8.1,a,f5.1,a,/))') &
            ' [PROFILE] solid_earth_update breakdown (per step, wall-clock):', &
            '   drift solve (band LU) =', 1.0e3_wp*t_dr/nstep, ' ms (', &
               100.0_wp*t_dr/max(t_upd,tiny(1.0_wp)), ' % of update)', &
            '   memory advance        =', 1.0e3_wp*t_mm/nstep, ' ms (', &
               100.0_wp*t_mm/max(t_upd,tiny(1.0_wp)), ' % of update)', &
            '   SLE + coupling (rest) =', 1.0e3_wp*(t_upd-t_dr-t_mm)/nstep, ' ms (', &
               100.0_wp*(t_upd-t_dr-t_mm)/max(t_upd,tiny(1.0_wp)), ' % of update)'
      end if
      if (nstep > 0) write(*,'(a,f7.1,a,f7.1,a)') &
         '   sub-steps/interval: n_accept=', real(se%stepper%n_accept,wp)/nstep, &
         '  n_solve=', real(se%stepper%n_solve,wp)/nstep, '  (per coupling step)'
      write(*,'(a,a)') ' fastearth: wrote ', trim(c%file_out)
      call fe_restart_write(se, se%time, folder=trim(rundir)//"/final")
      write(*,'(a,a)') ' fastearth: wrote restart ', trim(rundir)//'/final/fe_restart.nc'
      call solid_earth_finalize(se);  call sht_grid_destroy(sht)
   end subroutine fastearth_run

   pure function dir_of(path) result(d)
      !! Directory part of a path (everything before the last "/"); "." if none.
      character(len=*), intent(in)  :: path
      character(len=:), allocatable :: d
      integer :: i
      i = index(trim(path), "/", back=.true.)
      if (i > 0) then;  d = path(1:i-1);  else;  d = "."; end if
   end function dir_of

   ! --- input helpers ----------------------------------------------------------

   subroutine validate_inputs(c)
      !! Check every (file, variable) the run will actually read exists, before any
      !! work begins. Mirrors the conditional read paths (remap_input, i_eq) so a
      !! mistyped name_* / file_* config fails loudly with nonzero status instead of
      !! aborting silently mid-run (ncio's nc_read stop's with exit status 0 on a
      !! missing variable). 'what' names the offending config key in the message.
      type(fe_ctl_class), intent(in) :: c
      logical :: remap
      remap = c%remap_input

      ! forcing file: time axis + ice always; lon/lat only when remapping online.
      call require_var(c%file_forcing, c%name_time, 'file_forcing/name_time')
      call require_var(c%file_forcing, c%name_ice,  'file_forcing/name_ice')
      if (remap) then
         call require_var(c%file_forcing, c%name_lon, 'file_forcing/name_lon')
         call require_var(c%file_forcing, c%name_lat, 'file_forcing/name_lat')
      end if

      ! relaxed reference state (z_bed_eq, h_ice_eq), per i_eq -- see read_bed /
      ! setup_reference for which (file, var) each case reads.
      select case (c%i_eq)
      case (0)
         if (remap) then
            call require_var(c%file_forcing, c%name_zbed_eq, 'file_forcing/name_zbed_eq')
         else
            call require_var(c%file_ref, c%name_zbed_eq, 'file_ref/name_zbed_eq')
         end if
      case (1)
         call require_var(c%z_bed_ref_file, c%name_z_bed_ref, 'z_bed_ref_file/name_z_bed_ref')
         call require_var(c%h_ice_ref_file, c%name_h_ice_ref, 'h_ice_ref_file/name_h_ice_ref')
      case (2)
         call require_var(c%z_bed_eq_file, c%name_z_bed_ref, 'z_bed_eq_file/name_z_bed_ref')
         call require_var(c%h_ice_eq_file, c%name_h_ice_ref, 'h_ice_eq_file/name_h_ice_ref')
      case (3)
         call require_var(c%z_bed_ref_file,   c%name_z_bed_ref, 'z_bed_ref_file/name_z_bed_ref')
         call require_var(c%rsl_restart_file, c%name_rsl,       'rsl_restart_file/name_rsl')
         call require_var(c%h_ice_ref_file,   c%name_h_ice_ref, 'h_ice_ref_file/name_h_ice_ref')
      case default
         error stop 'fastearth_run: i_eq must be 0, 1, 2, or 3'
      end select
   end subroutine validate_inputs

   subroutine require_var(file, var, what)
      !! error stop (nonzero, flushed) if 'file' is unset/missing or lacks netCDF
      !! variable 'var'. 'what' is the config key(s) involved, named in the message.
      character(len=*), intent(in) :: file, var, what
      logical :: ok
      if (len_trim(file) == 0) then
         write(error_unit,'(a)') ' fastearth: input validation FAILED'
         write(error_unit,'(a)') '   config "'//trim(what)//'": file not set for the chosen i_eq'
         flush(error_unit)
         error stop 'fastearth_run: required input file not set (see message above)'
      end if
      inquire(file=trim(file), exist=ok)
      if (.not. ok) then
         write(error_unit,'(a)') ' fastearth: input validation FAILED'
         write(error_unit,'(a)') '   config "'//trim(what)//'": file not found'
         write(error_unit,'(a)') '   file: '//trim(file)
         flush(error_unit)
         error stop 'fastearth_run: input file not found (see message above)'
      end if
      if (.not. nc_exists_var(trim(file), trim(var))) then
         write(error_unit,'(a)') ' fastearth: input validation FAILED'
         write(error_unit,'(a)') '   config "'//trim(what)//'": variable "'//trim(var)//'" not found'
         write(error_unit,'(a)') '   file: '//trim(file)
         flush(error_unit)
         error stop 'fastearth_run: required netCDF variable not found (see message above)'
      end if
   end subroutine require_var

   subroutine build_remap(c, sht, rmap, nlon, nls)
      !! Read the source lon/lat axes from the forcing file and build the conservative
      !! lon-lat -> Gauss map. Returns the source dimensions for the work buffer.
      type(fe_ctl_class),   intent(in)    :: c
      type(sht_grid),       intent(in)    :: sht
      type(remap_ll_gauss),   intent(out)   :: rmap
      integer,              intent(out)   :: nlon, nls
      real(wp), allocatable :: lon_s(:), lat_s(:)
      nlon = nc_size(c%file_forcing, trim(c%name_lon))
      nls  = nc_size(c%file_forcing, trim(c%name_lat))
      allocate(lon_s(nlon), lat_s(nls))
      call nc_read(c%file_forcing, trim(c%name_lon), lon_s)
      call nc_read(c%file_forcing, trim(c%name_lat), lat_s)
      call remap_init(rmap, sht, lon_s, lat_s)
   end subroutine build_remap

   subroutine read_ice(c, rmap, sht, remap, k, src, h_ice)
      !! Read ice slice k onto the Gauss grid. remap: read lon-lat then conservatively
      !! remap (mass-conserving); else read the Gauss-grid slice directly.
      type(fe_ctl_class),   intent(in)    :: c
      type(remap_ll_gauss),   intent(in)    :: rmap
      type(sht_grid),       intent(in)    :: sht
      logical,              intent(in)    :: remap
      integer,              intent(in)    :: k
      real(wp),             intent(inout) :: src(:,:)
      real(wp),             intent(out)   :: h_ice(:,:)
      if (remap) then
         call nc_read(c%file_forcing, trim(c%name_ice), src, &
                      start=[1,1,k], count=[size(src,1), size(src,2), 1])
         call remap_to_gauss(rmap, sht, src, h_ice, conserve_mass=.true.)
      else
         call nc_read(c%file_forcing, trim(c%name_ice), h_ice, &
                      start=[1,1,k], count=[size(h_ice,1), size(h_ice,2), 1])
      end if
   end subroutine read_ice

   subroutine read_bed(c, rmap, sht, remap, k, src, z_bed)
      !! Read the reference bedrock. remap: slice k of name_zbed_eq from the forcing
      !! file (lon-lat-time), remapped (no mass rescale -- bed is geometry, not mass).
      !! else: the 2D name_zbed_eq from file_ref (legacy Gauss-grid reference file).
      type(fe_ctl_class),   intent(in)    :: c
      type(remap_ll_gauss),   intent(in)    :: rmap
      type(sht_grid),       intent(in)    :: sht
      logical,              intent(in)    :: remap
      integer,              intent(in)    :: k
      real(wp),             intent(inout) :: src(:,:)
      real(wp),             intent(out)   :: z_bed(:,:)
      if (remap) then
         call nc_read(c%file_forcing, trim(c%name_zbed_eq), src, &
                      start=[1,1,k], count=[size(src,1), size(src,2), 1])
         call remap_to_gauss(rmap, sht, src, z_bed, conserve_mass=.false.)
      else
         if (len_trim(c%file_ref) == 0) error stop 'fastearth_run: file_ref not set (remap_input=.false.)'
         call nc_read(c%file_ref, trim(c%name_zbed_eq), z_bed)
      end if
   end subroutine read_bed

   subroutine setup_reference(c, rmap, sht, remap, k0, src, ice_lgm, z_bed_eq, h_ice_eq)
      !! Fill the relaxed reference (z_bed_eq = SLE topo0, h_ice_eq) per c%i_eq,
      !! mirroring CLIMBER-X i_equilibrium. Reference files are lon-lat and remapped
      !! online with a per-file conservative map (their grid differs from the forcing).
      type(fe_ctl_class),   intent(in)    :: c
      type(remap_ll_gauss),   intent(in)    :: rmap
      type(sht_grid),       intent(in)    :: sht
      logical,              intent(in)    :: remap
      integer,              intent(in)    :: k0
      real(wp),             intent(inout) :: src(:,:)
      real(wp),             intent(in)    :: ice_lgm(:,:)
      real(wp),             intent(out)   :: z_bed_eq(:,:), h_ice_eq(:,:)
      real(wp), allocatable :: rsl_r(:,:)
      select case (c%i_eq)
      case (0)
         call read_bed(c, rmap, sht, remap, k0, src, z_bed_eq)   ! data bed at the start slice
         h_ice_eq = ice_lgm
      case (1)
         call read_ref2d(c, sht, remap, c%z_bed_ref_file, c%name_z_bed_ref, .false., z_bed_eq)
         call read_ref2d(c, sht, remap, c%h_ice_ref_file, c%name_h_ice_ref, .true.,  h_ice_eq)
      case (2)
         call read_ref2d(c, sht, remap, c%z_bed_eq_file, c%name_z_bed_ref, .false., z_bed_eq)
         call read_ref2d(c, sht, remap, c%h_ice_eq_file, c%name_h_ice_ref, .true.,  h_ice_eq)
      case (3)
         allocate(rsl_r(size(z_bed_eq,1), size(z_bed_eq,2)))
         call read_ref2d(c, sht, remap, c%z_bed_ref_file,   c%name_z_bed_ref, .false., z_bed_eq)
         call read_ref2d(c, sht, remap, c%rsl_restart_file, c%name_rsl,       .false., rsl_r)
         call read_ref2d(c, sht, remap, c%h_ice_ref_file,   c%name_h_ice_ref, .true.,  h_ice_eq)
         z_bed_eq = z_bed_eq + rsl_r
      case default
         error stop 'fastearth_run: i_eq must be 0, 1, 2, or 3'
      end select
   end subroutine setup_reference

   subroutine read_ref2d(c, sht, remap, file, varname, conserve, field)
      !! Read a 2D lon-lat reference field onto the Gauss grid. If the file is
      !! already on the Gauss grid (its lon/lat dims match nphi/nlat — e.g. the
      !! canonical reference at this run's resolution) it is read directly. Otherwise
      !! a per-file conservative map is built from the file's own axes and applied
      !! (conserve=.true. for ice) — the weights are cached (fe_remap), so the build
      !! cost is paid once. In legacy (remap=.false.) mode the field is read directly.
      type(fe_ctl_class),   intent(in)  :: c
      type(sht_grid),       intent(in)  :: sht
      logical,              intent(in)  :: remap, conserve
      character(len=*),     intent(in)  :: file, varname
      real(wp),             intent(out) :: field(:,:)
      type(remap_ll_gauss)    :: m
      real(wp), allocatable :: lon_s(:), lat_s(:), buf(:,:)
      integer :: nlon, nls
      if (len_trim(file) == 0) &
         error stop 'fastearth_run: reference file not set for the chosen i_eq'
      if (remap) then
         nlon = nc_size(file, trim(c%name_lon));  nls = nc_size(file, trim(c%name_lat))
         if (nlon == sht%nphi .and. nls == sht%nlat) then
            call nc_read(file, trim(varname), field)   ! already on the Gauss grid
            return
         end if
         allocate(lon_s(nlon), lat_s(nls), buf(nlon, nls))
         call nc_read(file, trim(c%name_lon), lon_s)
         call nc_read(file, trim(c%name_lat), lat_s)
         call nc_read(file, trim(varname),    buf)
         call remap_init(m, sht, lon_s, lat_s)
         call remap_to_gauss(m, sht, buf, field, conserve_mass=conserve)
      else
         call nc_read(file, trim(varname), field)
      end if
   end subroutine read_ref2d

   ! --- grid + window ----------------------------------------------------------

   subroutine build_grid(p, sht)
      !! Build the Gauss-Legendre transform grid from p. When nlat/nphi are unset
      !! (<=0) default to a de-aliased grid (nlat=2 lmax+2, nphi=4 lmax) sized for
      !! the SLE's quadratic ocean-function product.
      type(fe_param_class), intent(in)            :: p
      type(sht_grid),       intent(inout), target :: sht
      integer :: nlat, nphi
      nlat = p%nlat;  if (nlat <= 0) nlat = 2*p%lmax + 2
      nphi = p%nphi;  if (nphi <= 0) nphi = 4*p%lmax
      if (p%mmax >= 0 .and. p%eps_polar > 0.0_wp) then
         call sht_grid_init(sht, p%lmax, nlat=nlat, nphi=nphi, mmax=p%mmax, mres=p%mres, eps_polar=p%eps_polar)
      else if (p%mmax >= 0) then
         call sht_grid_init(sht, p%lmax, nlat=nlat, nphi=nphi, mmax=p%mmax, mres=p%mres)
      else if (p%eps_polar > 0.0_wp) then
         call sht_grid_init(sht, p%lmax, nlat=nlat, nphi=nphi, mres=p%mres, eps_polar=p%eps_polar)
      else
         call sht_grid_init(sht, p%lmax, nlat=nlat, nphi=nphi, mres=p%mres)
      end if
   end subroutine build_grid

   subroutine select_window(tsec, t_init, t_end, k0, k1)
      !! First/last forcing indices inside [t_init, t_end] (all in seconds).
      real(wp), intent(in)  :: tsec(:), t_init, t_end
      integer,  intent(out) :: k0, k1
      integer :: k
      k0 = 0;  k1 = 0
      do k = 1, size(tsec)
         if (tsec(k) >= t_init .and. tsec(k) <= t_end) then
            if (k0 == 0) k0 = k
            k1 = k
         end if
      end do
      if (k0 == 0) error stop 'fastearth_run: no forcing times within [time_init, time_end]'
   end subroutine select_window

end module fe_drive
