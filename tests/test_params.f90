program test_params
   !! fe_params: load a &fe3d namelist into fe_param_class and build the earth
   !! model from it. A sparse user file is overlaid on the complete, shipped
   !! defaults file (input/fastearth3d_defaults.nml) — the yelmo defaults_file
   !! convention. Checks
   !! (1) scalar / string / logical overrides, (2) the YEARS->seconds conversion
   !! of the time fields, (3) custom per-layer earth assembly, (4) the named
   !! built-in path (defaults), and (5) an un-overridden value falling through to
   !! the defaults file.
   use fe_precision,       only: wp
   use fe_constants,       only: sec_per_year
   use fe_params,          only: fe_param_class, fe_par_load
   use fe_earth_structure, only: earth_n_layers, earth_model, build_earth, RHEOL_MAXWELL, RHEOL_FLUID
   implicit none

   character(len=*), parameter :: NML  = "obj/test_params.nml"
   character(len=*), parameter :: DEFS = "input/fastearth3d_defaults.nml"   ! shipped complete &fe3d defaults
   type(fe_param_class) :: p, pdef
   type(earth_model)    :: em, emdef
   integer :: u
   logical :: ok
   real(wp) :: tol

   ok = .true.;  tol = 1.0e-6_wp

   ! --- write a SPARSE user file (custom 2-layer earth) overriding the defaults -
   open(newunit=u, file=NML, status="replace", action="write")
   write(u,'(a)') "&fe3d"
   write(u,'(a)') "    lmax    = 32"
   write(u,'(a)') '    earth   = "custom"'
   write(u,'(a)') "    n_layer = 2"
   write(u,'(a)') "    r_bot    = 5701.0e3, 0.0"
   write(u,'(a)') "    r_top    = 6371.0e3, 3480.0e3"
   write(u,'(a)') "    rho      = 4000.0, 10000.0"
   write(u,'(a)') "    mu       = 1.0e11, 0.0"
   write(u,'(a)') "    eta      = 1.0e21, 0.0"
   write(u,'(a)') "    rheology = 1, 2"
   write(u,'(a)') '    scheme       = "be"'
   write(u,'(a)') "    sle_n_outer  = 5"
   write(u,'(a)') "    sle_subgrid  = .false."
   write(u,'(a)') "    equil_time_max = 2000.0"
   write(u,'(a)') "    dt_init      = 500.0"
   write(u,'(a)') "    rotation     = .true."
   write(u,'(a)') "/"
   close(u)

   call fe_par_load(p, NML, defaults_file=DEFS)

   ! --- (1) scalar / string / logical overrides -------------------------------
   call check_int("lmax",        p%lmax,        32)
   call check_str("earth",       trim(p%earth), "custom")
   call check_int("n_layer",     p%n_layer,     2)
   call check_str("scheme",      trim(p%scheme), "be")
   call check_int("sle_n_outer", p%sle_n_outer, 5)
   call check_log("sle_subgrid", p%sle_subgrid, .false.)
   call check_log("rotation",    p%rotation,    .true.)

   ! --- (2) YEARS -> seconds conversion of the time fields ---------------------
   call check_real("equil_time_max [s]", p%equil_time_max, 2000.0_wp*sec_per_year)
   call check_real("dt_init [s]",        p%dt_init,         500.0_wp*sec_per_year)

   ! --- (3) custom per-layer earth assembly ------------------------------------
   em = build_earth(p)
   call check_int("custom n_layers", earth_n_layers(em), 2)
   call check_str("custom name",     em%name,       "custom")
   call check_int("layer1 rheology", em%layers(1)%rheology, RHEOL_MAXWELL)
   call check_int("layer2 rheology", em%layers(2)%rheology, RHEOL_FLUID)
   call check_real("layer1 mu [Pa]", em%layers(1)%mu, 1.0e11_wp)

   ! --- (4) named built-in path (in-code defaults) -----------------------------
   emdef = build_earth(pdef)                     ! pdef left at defaults => "M3-L70-V01"
   call check_str("default earth",    trim(pdef%earth), "M3-L70-V01")
   call check_int("M3-L70-V01 layers", earth_n_layers(emdef), 5)

   ! --- (5) un-overridden values fall through to the defaults file -------------
   call check_int("fallthrough max_couple_iter", p%max_couple_iter, 20)
   call check_int("fallthrough sle_n_inner",     p%sle_n_inner,     20)

   write(*,'(a)') ''
   if (ok) then
      write(*,'(a)') ' PASS: fe_par_load fills the record (years->s) and build_earth'
      write(*,'(a)') '       assembles both custom and named earth models'
   else
      write(*,'(a)') ' FAIL: fe_params did not all pass'
      error stop 1
   end if

contains

   subroutine check_int(name, got, want)
      character(len=*), intent(in) :: name
      integer,          intent(in) :: got, want
      if (got /= want) then
         write(*,'(a,a,a,i0,a,i0)') '   FAIL ', name, ': got ', got, ' want ', want
         ok = .false.
      end if
   end subroutine check_int

   subroutine check_real(name, got, want)
      character(len=*), intent(in) :: name
      real(wp),         intent(in) :: got, want
      if (abs(got - want) > tol*max(abs(want), 1.0_wp)) then
         write(*,'(a,a,a,es12.4,a,es12.4)') '   FAIL ', name, ': got ', got, ' want ', want
         ok = .false.
      end if
   end subroutine check_real

   subroutine check_str(name, got, want)
      character(len=*), intent(in) :: name, got, want
      if (got /= want) then
         write(*,'(a,a,a,a,a,a)') '   FAIL ', name, ': got "', got, '" want "', want//'"'
         ok = .false.
      end if
   end subroutine check_str

   subroutine check_log(name, got, want)
      character(len=*), intent(in) :: name
      logical,          intent(in) :: got, want
      if (got .neqv. want) then
         write(*,'(a,a,a,l1,a,l1)') '   FAIL ', name, ': got ', got, ' want ', want
         ok = .false.
      end if
   end subroutine check_log

end program test_params
