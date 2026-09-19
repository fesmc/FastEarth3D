program test_restart
   !! netCDF restart round-trip for the coupling state (yelmo-convention I/O via
   !! fe_io: a time axis lets several snapshots share one file). Run for BOTH
   !! response kinds — RESP_VE (Maxwell memory tensor + σ_n) and RESP_MODAL (the
   !! per-(l,m) modal amplitudes φ; the spectrum is rebuilt by init). Checks:
   !!
   !!   (1) direct state restore — read the last snapshot into a fresh model and
   !!       its z_bed/rsl match what was written;
   !!   (2) bit-for-bit continuation — restore the EARLIER snapshot's memory
   !!       (tau_*/φ, plus the rotation state m + channel memory rot_*) + adaptive-Δt
   !!       seed and step forward; the trajectory reproduces the uninterrupted run
   !!       exactly (the restored prognostic state is all that is needed);
   !!   (3) multi-snapshot file — two snapshots at different times coexist, and a
   !!       specific earlier time can be selected on read.
   use fe_precision,       only: wp
   use fe_constants,       only: pi, sec_per_year
   use fe_params,          only: fe_param_class
   use fe_radial_fe,       only: radial_fe_finalize
   use fe_sht,             only: sht_grid, sht_grid_init, sht_grid_destroy
   use fe_coupling,        only: solid_earth_finalize, solid_earth_update, solid_earth_init, solid_earth
   use fe_io,              only: fe_restart_write, fe_restart_read
   use ncio,               only: nc_size
   implicit none

   integer, parameter :: LMAX = 12, K1 = 3, K2 = 3
   type(sht_grid), target :: sht
   real(wp), allocatable  :: z_bed_eq(:,:), h_ice_eq(:,:), h_ice(:,:)
   logical  :: ok

   ok = .true.
   call sht_grid_init(sht, LMAX, nlat=2*LMAX, nphi=4*LMAX)
   allocate(z_bed_eq(sht%nphi,sht%nlat), h_ice_eq(sht%nphi,sht%nlat), &
            h_ice(sht%nphi,sht%nlat))
   call make_fields(z_bed_eq, h_ice_eq, h_ice)

   call roundtrip("ve",    "obj/test_restart_ve.nc",    ok)
   call roundtrip("modal", "obj/test_restart_modal.nc", ok)

   call sht_grid_destroy(sht)

   write(*,'(a)') ''
   if (ok) then
      write(*,'(a)') ' PASS: netCDF restart restores state, continues bit-for-bit,'
      write(*,'(a)') '       and stores several snapshots in one file (RESP_VE + RESP_MODAL)'
   else
      write(*,'(a)') ' FAIL: restart round-trip did not all pass'
      call radial_fe_finalize()
      error stop 1
   end if
   call radial_fe_finalize()

contains

   subroutine roundtrip(resp, file, ok)
      !! Full restart round-trip for one response kind on the shared grid/fields.
      character(len=*),  intent(in)    :: resp, file
      logical,           intent(inout) :: ok
      type(fe_param_class) :: p
      type(solid_earth)    :: a, b, c
      real(wp), allocatable :: a6_zbed(:,:), a6_rsl(:,:)
      real(wp) :: t1, t2, d_restore, d_continue, dt_yr
      integer  :: step, nt

      dt_yr           = 1.0e3_wp          ! interval per update [years]
      p%lmax = LMAX;  p%nlat = 2*LMAX;  p%nphi = 4*LMAX  ! model builds its own grid from par
      p%earth_response = resp             ! "ve" (memory tensor) or "modal" (φ amplitudes)
      p%rotation      = .true.            ! rotation on: also round-trips the polar
                                          ! motion m + both channels' memory (rot_*)

      ! === reference run A ====================================================
      a%par = p; call solid_earth_init(a, z_bed_eq, h_ice_eq)
      do step = 1, K1
         call solid_earth_update(a, h_ice, dt_yr)
      end do
      t1 = a%time
      call fe_restart_write(a, t1, filename=file, init=.true.)   ! snapshot 1 (memory @ K1)

      do step = 1, K2
         call solid_earth_update(a, h_ice, dt_yr)
      end do
      t2 = a%time
      call fe_restart_write(a, t2, filename=file, init=.false.)  ! snapshot 2 (state @ K1+K2)
      a6_zbed = a%z_bed;  a6_rsl = a%rsl

      write(*,'(a,a,a,i0,a,f6.2,a,f6.2)') ' restart [', trim(resp), ']: lmax=', LMAX, &
           '   t1=', t1*1.0e-3_wp, ' kyr   t2=', t2*1.0e-3_wp

      ! === (3) multi-snapshot file ===========================================
      nt = nc_size(file, "time")
      write(*,'(a,i0)') '   snapshots in file: ', nt
      if (nt /= 2) then
         write(*,'(a)') '   FAIL: expected two time snapshots in one file'; ok = .false.
      end if

      ! === (1) direct state restore (default = last snapshot, t2) =============
      b%par = p; call solid_earth_init(b, z_bed_eq, h_ice_eq)
      call fe_restart_read(b, file)
      d_restore = max(maxval(abs(b%z_bed - a6_zbed)), maxval(abs(b%rsl - a6_rsl)))
      write(*,'(a,es11.2)') '   (1) state restore  max|B - A|     =', d_restore
      if (d_restore > 1.0e-9_wp) then
         write(*,'(a)') '   FAIL: restored state does not match the written state'; ok = .false.
      end if

      ! === (2) bit-for-bit continuation from the earlier snapshot (t1) =========
      c%par = p; call solid_earth_init(c, z_bed_eq, h_ice_eq)
      call fe_restart_read(c, file, time=t1)                ! restore memory @ K1
      do step = 1, K2
         call solid_earth_update(c, h_ice, dt_yr)     ! same load, K2 steps
      end do
      d_continue = max(maxval(abs(c%z_bed - a6_zbed)), maxval(abs(c%rsl - a6_rsl)))
      write(*,'(a,es11.2)') '   (2) continuation   max|C - A|     =', d_continue
      if (d_continue > 1.0e-9_wp) then
         write(*,'(a)') '   FAIL: restarted trajectory diverges from the uninterrupted run'
         ok = .false.
      end if

      call solid_earth_finalize(a);  call solid_earth_finalize(b);  call solid_earth_finalize(c)
   end subroutine roundtrip

   subroutine make_fields(z_bed_eq, h_ice_eq, h_ice)
      !! Polar land cap (colat<50°, +500 m) over deep ocean (−4000 m); no
      !! reference ice; a 2 km grounded ice load on colat<30°.
      real(wp), intent(out) :: z_bed_eq(:,:), h_ice_eq(:,:), h_ice(:,:)
      integer  :: i, j
      real(wp) :: thd
      do j = 1, sht%nlat
         thd = sht%colat(j)*180.0_wp/pi
         do i = 1, sht%nphi
            z_bed_eq(i,j) = merge(500.0_wp, -4000.0_wp, thd < 50.0_wp)
            h_ice(i,j)    = merge(2000.0_wp, 0.0_wp,     thd < 30.0_wp)
         end do
      end do
      h_ice_eq = 0.0_wp
   end subroutine make_fields

end program test_restart
