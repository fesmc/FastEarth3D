program test_flotation
   !! Grounded-ice flotation in the SLE ocean function (design.md item 2).
   !!
   !! A cell below sea level is ocean only where the ice floats; ice thick
   !! enough to ground on the bed (ρ_i·I ≥ −ρ_w·topo) is land. We drive the SLE
   !! with NO ice-change (d_ice = 0, so S ≡ 0 and the classification comes purely
   !! from the reference bathymetry + the absolute ice field) and check the
   !! returned ocean function C against the criterion recomputed independently.
   !!
   !! Zones (by colatitude) exercise the headline cases AND the depth dependence
   !! — the SAME ice thickness must ground over shallow water yet float over
   !! deep water, because the test compares the ice draft to the water column,
   !! not to a fixed thickness:
   !!
   !!   A  colat<20   topo=-4000  ice=5000  ->  grounds (C=0): thick ice / deep
   !!   B  20..40     topo=-4000  ice= 100  ->  floats  (C=1): thin ice  / deep
   !!   C  40..60     topo=-4000  ice=   0  ->  ocean   (C=1): ice-free   / deep
   !!   D  60..80     topo=-2000  ice=3000  ->  grounds (C=0): 3 km ice, shallow
   !!   E  colat>80   topo=-4000  ice=3000  ->  floats  (C=1): 3 km ice, deep
   use fe_precision, only: wp
   use fe_constants, only: rho_ice, rho_water, pi
   use fe_response,  only: response, response_init_elastic, response_init_ve, response_init_null
   use fe_sht,       only: sht_grid, sht_grid_init, sht_grid_destroy
   use fe_sle,       only: sle_solve, sle_solver, sle_result
   implicit none

   integer, parameter :: LMAX = 16
   type(sht_grid)      :: sht
   type(sle_solver)    :: sle
   type(sle_result)    :: res
   type(response) :: resp
   real(wp), allocatable :: topo0(:,:), d_ice(:,:), ice(:,:), S(:,:), C(:,:), &
                            Cexpect(:,:)
   integer :: i, j, nmismatch, n_grounded, n_ocean
   real(wp) :: smax
   logical  :: ok

   ok = .true.
   call sht_grid_init(sht, LMAX, nlat=2*LMAX, nphi=4*LMAX)
   allocate(topo0(sht%nphi,sht%nlat), d_ice(sht%nphi,sht%nlat), &
            ice(sht%nphi,sht%nlat), S(sht%nphi,sht%nlat), C(sht%nphi,sht%nlat), &
            Cexpect(sht%nphi,sht%nlat))

   call make_fields(topo0, ice)
   d_ice = 0.0_wp                       ! no load change: isolate classification

   call sle_solve(sle, sht, resp, d_ice, ice, topo0, S, C, res)

   ! S must be identically zero with no ice-change forcing
   smax = maxval(abs(S))

   ! Independently recompute the expected ocean function and compare cell-by-cell
   where (topo0 < 0.0_wp .and. rho_ice*ice < -rho_water*topo0)
      Cexpect = 1.0_wp
   elsewhere
      Cexpect = 0.0_wp
   end where

   nmismatch = 0;  n_grounded = 0;  n_ocean = 0
   do j = 1, sht%nlat
      do i = 1, sht%nphi
         if (abs(C(i,j) - Cexpect(i,j)) > 1.0e-12_wp) nmismatch = nmismatch + 1
         if (C(i,j) < 0.5_wp) then
            n_grounded = n_grounded + 1
         else
            n_ocean = n_ocean + 1
         end if
      end do
   end do

   write(*,'(a)')          ' Grounded-ice flotation in the ocean function'
   write(*,'(a,es10.2)')   '      max|S| with d_ice=0           =', smax
   write(*,'(a,i6)')       '      cells classified grounded (C=0)=', n_grounded
   write(*,'(a,i6)')       '      cells classified ocean    (C=1)=', n_ocean
   write(*,'(a,i6)')       '      C vs. recomputed criterion mismatches =', nmismatch

   ! --- headline assertions ---------------------------------------------------
   call zone_check('A thick ice / deep  -> grounded', 10.0_wp, 0.0_wp, ok)
   call zone_check('B thin  ice / deep  -> ocean   ', 30.0_wp, 1.0_wp, ok)
   call zone_check('C ice-free  / deep  -> ocean   ', 50.0_wp, 1.0_wp, ok)
   call zone_check('D 3km ice / shallow -> grounded', 70.0_wp, 0.0_wp, ok)
   call zone_check('E 3km ice / deep    -> ocean   ', 85.0_wp, 1.0_wp, ok)

   if (smax > 1.0e-12_wp) then
      write(*,'(a)') '      FAIL: S not identically zero with no ice-change'; ok = .false.
   end if
   if (nmismatch /= 0) then
      write(*,'(a)') '      FAIL: ocean function disagrees with flotation criterion'; ok = .false.
   end if
   if (n_grounded == 0 .or. n_ocean == 0) then
      write(*,'(a)') '      FAIL: degenerate classification (all one class)'; ok = .false.
   end if

   write(*,'(a)') ''
   if (ok) then
      write(*,'(a)') ' PASS: grounded ice is excluded from the ocean and the'
      write(*,'(a)') '       criterion is depth-dependent (D grounds, E floats)'
   else
      write(*,'(a)') ' FAIL: flotation classification incorrect'
      call sht_grid_destroy(sht)
      error stop 1
   end if
   call sht_grid_destroy(sht)

contains

   subroutine make_fields(topo0, ice)
      real(wp), intent(out) :: topo0(:,:), ice(:,:)
      integer  :: i, j
      real(wp) :: thd
      do j = 1, sht%nlat
         thd = sht%colat(j)*180.0_wp/pi          ! colatitude [deg]
         do i = 1, sht%nphi
            if (thd < 20.0_wp) then
               topo0(i,j) = -4000.0_wp;  ice(i,j) = 5000.0_wp   ! A
            else if (thd < 40.0_wp) then
               topo0(i,j) = -4000.0_wp;  ice(i,j) =  100.0_wp   ! B
            else if (thd < 60.0_wp) then
               topo0(i,j) = -4000.0_wp;  ice(i,j) =    0.0_wp   ! C
            else if (thd < 80.0_wp) then
               topo0(i,j) = -2000.0_wp;  ice(i,j) = 3000.0_wp   ! D
            else
               topo0(i,j) = -4000.0_wp;  ice(i,j) = 3000.0_wp   ! E
            end if
         end do
      end do
   end subroutine make_fields

   subroutine zone_check(label, colat_deg, cexp, ok)
      !! Check the ocean function at the latitude row nearest a colatitude.
      character(len=*), intent(in)    :: label
      real(wp),         intent(in)    :: colat_deg, cexp
      logical,          intent(inout) :: ok
      integer  :: j, jbest
      real(wp) :: d, dbest
      jbest = 1;  dbest = huge(1.0_wp)
      do j = 1, sht%nlat
         d = abs(sht%colat(j)*180.0_wp/pi - colat_deg)
         if (d < dbest) then;  dbest = d;  jbest = j;  end if
      end do
      write(*,'(a,a,a,f4.1,a,f4.1)') '      ', label, '  got C=', C(1,jbest), &
                                     '  expect ', cexp
      if (abs(C(1,jbest) - cexp) > 1.0e-12_wp) then
         write(*,'(a)') '      FAIL: zone misclassified'; ok = .false.
      end if
   end subroutine zone_check

end program test_flotation
