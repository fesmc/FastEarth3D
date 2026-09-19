program test_spinup
   !! LGM-memory spin-up (solid_earth_spinup): hold a start-slice (LGM) ice load and
   !! relax toward isostatic equilibrium BEFORE the transient, HOLDING the reference as
   !! the datum. Uses the default 1-D earth (no 3-D viscosity file needed). Checks:
   !!   (1) full-model phase (equil_time_max): relaxes the bed far more than a single
   !!       cold coupling step, and leaves it near-stationary (a follow-up step moves it
   !!       less than a cold step would);
   !!   (2) pre_spinup_1d phase also runs and drives appreciable subsidence.
   use fe_precision, only: wp
   use fe_constants, only: pi, sec_per_year
   use fe_params,    only: fe_param_class
   use fe_radial_fe, only: radial_fe_finalize
   use fe_sht,       only: sht_grid, sht_grid_init, sht_grid_destroy
   use fe_coupling,  only: solid_earth, solid_earth_init, solid_earth_update, &
                           solid_earth_spinup, solid_earth_finalize
   implicit none

   integer, parameter :: LMAX = 16
   type(sht_grid), target :: sht
   type(fe_param_class)   :: p
   type(solid_earth)      :: se
   real(wp), allocatable  :: z_bed_eq(:,:), h_ice_eq(:,:), h_ice(:,:)
   real(wp) :: cold_step, spin_sub, follow_move, bed0
   integer  :: i, j, jice
   logical  :: ok

   ok = .true.
   call sht_grid_init(sht, LMAX, nlat=2*LMAX, nphi=4*LMAX)
   allocate(z_bed_eq(sht%nphi,sht%nlat), h_ice_eq(sht%nphi,sht%nlat), h_ice(sht%nphi,sht%nlat))
   do j = 1, sht%nlat
      do i = 1, sht%nphi
         z_bed_eq(i,j) = merge(500.0_wp, -4000.0_wp, sht%colat(j)*180.0_wp/pi < 50.0_wp)
         h_ice(i,j)    = merge(2000.0_wp, 0.0_wp,     sht%colat(j)*180.0_wp/pi < 30.0_wp)
      end do
   end do
   h_ice_eq = 0.0_wp
   jice = nearest_row(15.0_wp)

   p%lmax = LMAX;  p%nlat = 2*LMAX;  p%nphi = 4*LMAX

   ! reference: a single cold 2 kyr coupling step (elastic + a little viscous)
   se%par = p;  call solid_earth_init(se, z_bed_eq, h_ice_eq)
   call solid_earth_update(se, h_ice, 2.0e3_wp)
   cold_step = z_bed_eq(1,jice) - se%z_bed(1,jice)
   call solid_earth_finalize(se)

   ! (1) full-model spin-up: relax to (near) equilibrium under the held load
   p%pre_spinup_1d  = .false.
   p%equil_time_max = 5.0e5_wp*sec_per_year       ! 500 kyr cap (par stores SI; >> Maxwell time)
   se%par = p;  call solid_earth_init(se, z_bed_eq, h_ice_eq)
   call solid_earth_spinup(se, h_ice)
   spin_sub = z_bed_eq(1,jice) - se%z_bed(1,jice)
   bed0 = se%z_bed(1,jice)
   call solid_earth_update(se, h_ice, 2.0e3_wp)   ! a further step under the same load
   follow_move = abs(se%z_bed(1,jice) - bed0)
   call solid_earth_finalize(se)

   write(*,'(a)') ' spin-up: 1-D earth, 2 km polar cap held to equilibrium'
   write(*,'(a,f12.3,a)') '   cold single-step subsidence   =', cold_step,   ' m'
   write(*,'(a,f12.3,a)') '   spin-up subsidence            =', spin_sub,    ' m'
   write(*,'(a,f12.3,a)') '   follow-up step move (post-spin)=', follow_move, ' m'

   if (spin_sub < 50.0_wp) then
      write(*,'(a)') '   FAIL: spin-up produced little subsidence'; ok = .false.
   end if
   if (spin_sub < 1.5_wp*cold_step) then
      write(*,'(a)') '   FAIL: spin-up did not relax appreciably beyond a single cold step'; ok = .false.
   end if
   ! the strong equilibration signal: a further step barely moves the bed (<< one cold step)
   if (follow_move > 0.1_wp*cold_step) then
      write(*,'(a)') '   FAIL: post-spin-up state is not near-stationary'; ok = .false.
   end if

   ! (2) pre_spinup_1d phase (1-D pre-equilibration only, no full-model phase)
   p%pre_spinup_1d  = .true.
   p%equil_time_max = 0.0_wp
   se%par = p;  call solid_earth_init(se, z_bed_eq, h_ice_eq)
   call solid_earth_spinup(se, h_ice)
   write(*,'(a,f12.3,a)') '   pre_spinup_1d subsidence       =', &
        z_bed_eq(1,jice) - se%z_bed(1,jice), ' m'
   if (z_bed_eq(1,jice) - se%z_bed(1,jice) < 50.0_wp) then
      write(*,'(a)') '   FAIL: pre_spinup_1d produced little subsidence'; ok = .false.
   end if
   call solid_earth_finalize(se)

   call sht_grid_destroy(sht)
   write(*,'(a)') ''
   if (ok) then
      write(*,'(a)') ' PASS: solid_earth_spinup relaxes the held load to equilibrium'
   else
      write(*,'(a)') ' FAIL: spin-up'
      call radial_fe_finalize()
      error stop 1
   end if
   call radial_fe_finalize()

contains

   integer function nearest_row(colat_deg) result(jbest)
      real(wp), intent(in) :: colat_deg
      integer  :: j
      real(wp) :: d, dbest
      jbest = 1;  dbest = huge(1.0_wp)
      do j = 1, sht%nlat
         d = abs(sht%colat(j)*180.0_wp/pi - colat_deg)
         if (d < dbest) then;  dbest = d;  jbest = j;  end if
      end do
   end function nearest_row

end program test_spinup
