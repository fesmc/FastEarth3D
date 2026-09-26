program test_visc_lid
   !! The lid rule of response_enable_lateral_visc_from_nodes: a Maxwell element
   !! lying wholly above lid_depth whose log10 η exceeds lid_log10max is set to
   !! log10_cap, and nothing else changes.
   !!
   !! M3-L70-V01 with its 70 km lithosphere made Maxwell (η = 1e30; an elastic
   !! element carries no lateral viscosity, so the rule could not act on it). The
   !! field is log10 η = 25 on half the grid and 21 on the other, at every depth.
   !! With lid_depth = 70 km and lid_log10max = 22: lid elements read 30 where the
   !! field is 25 and keep 21; mantle elements keep 25 and 21. Without the lid
   !! arguments, every element keeps the field.
   use vilma_precision,       only: wp
   use vilma_earth_structure, only: earth_model, build_M3L70V01, RHEOL_MAXWELL
   use vilma_radial_fe,       only: radial_fe_finalize
   use vilma_response,        only: response, response_init_ve, response_destroy, &
                                    response_enable_lateral_visc_from_nodes
   use vilma_sht,             only: sht_grid, sht_grid_init
   implicit none

   integer,  parameter :: LMAX = 16
   real(wp), parameter :: YR = 3.15576e7_wp, DT = 10.0_wp*YR
   real(wp), parameter :: HI = 25.0_wp, LO = 21.0_wp, THR = 22.0_wp, CAP = 30.0_wp
   real(wp), parameter :: LID = 70.0e3_wp, TOL = 1.0e-9_wp

   type(sht_grid)    :: sht
   type(earth_model) :: e
   type(response)    :: ve
   real(wp), allocatable :: vn(:,:)
   logical :: ok
   integer :: nlid

   ok = .true.
   call sht_grid_init(sht, LMAX, nlat=2*LMAX, nphi=2*LMAX+2, mmax=LMAX)
   e = build_M3L70V01()
   e%layers(1)%rheology = RHEOL_MAXWELL
   e%layers(1)%eta      = 1.0e30_wp

   ! rule on
   call response_init_ve(ve, e, sht, DT)
   call half_field(vn)
   call response_enable_lateral_visc_from_nodes(ve, sht, vn, lid_depth=LID, &
        lid_log10max=THR, log10_cap=CAP)
   call check(.true., nlid)
   if (nlid == 0) then
      write(*,'(a)') ' FAIL: no Maxwell element lies above the lid base'
      ok = .false.
   end if
   call response_destroy(ve)

   ! rule off: the field reaches every element unchanged
   call response_init_ve(ve, e, sht, DT)
   call response_enable_lateral_visc_from_nodes(ve, sht, vn)
   call check(.false., nlid)
   call response_destroy(ve)

   call radial_fe_finalize()
   if (.not. ok) error stop 'test_visc_lid FAILED'
   write(*,'(a,i0,a)') ' PASS: lid rule caps only lid elements above the threshold (', nlid, &
        ' lid elements)'
   write(*,'(a)') ' test_visc_lid PASSED'

contains

   subroutine half_field(v)
      !! log10 η = HI on the first half of the longitudes, LO on the rest, all depths.
      real(wp), allocatable, intent(out) :: v(:,:)
      integer :: i, j
      allocate(v(sht%nphi*sht%nlat, ve%nr))
      do j = 1, sht%nlat
         do i = 1, sht%nphi
            v(i + (j-1)*sht%nphi, :) = merge(HI, LO, i <= sht%nphi/2)
         end do
      end do
   end subroutine half_field

   subroutine check(rule, nl)
      !! Effective log10 η of every Maxwell element, from its rate μ/η, against
      !! what the rule (on or off) should have produced.
      logical, intent(in)  :: rule
      integer, intent(out) :: nl
      real(wp) :: want_hi, eff
      logical  :: in_lid
      integer  :: el, i, j
      nl = 0
      do el = 1, ve%ne
         if (ve%MkPerDt(el) == 0.0_wp) cycle
         in_lid = ve%r(el) >= ve%r(ve%nr) - LID - 1.0_wp
         if (in_lid) nl = nl + 1
         want_hi = merge(CAP, HI, rule .and. in_lid)
         do j = 1, sht%nlat
            do i = 1, sht%nphi
               eff = log10(ve%mu(el)/ve%MkPerDt3(i,j,el))
               if (abs(eff - merge(want_hi, LO, i <= sht%nphi/2)) > TOL) then
                  write(*,'(a,l1,a,i0,a,l1,a,2i4,a,f8.4)') ' FAIL (rule=', rule, ') element ', el, &
                       ' lid=', in_lid, ' point', i, j, ': log10 eta ', eff
                  ok = .false.
                  return
               end if
            end do
         end do
      end do
   end subroutine check

end program test_visc_lid
