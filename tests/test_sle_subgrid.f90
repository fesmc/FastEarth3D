program test_sle_subgrid
   !! Validate the subgrid (sloping-coast) ocean water load of fe_sle against an
   !! independent calculation, and contrast it with the binary coastline.
   !!
   !! With a NULL (rigid) response the SLE reduces to pure eustatic redistribution:
   !! u = N = 0, so the sea-surface change is a uniform dphi and the ocean function
   !! migrates as O = (topo0 < dphi). On a SLOPING basin the volume swept as the
   !! coastline moves depends on the bathymetry, so the mass-conserving dphi is the
   !! root of a 1-D volume balance. The subgrid load s = C*dphi − topo0*(C − C0) is
   !! exactly the water-column change max(0,dphi−topo0) − max(0,−topo0), so the
   !! subgrid solver's dphi MUST equal that volume-balance root. The binary load
   !! C*dphi ignores the swept bathymetry and gives a different (wrong) dphi.
   !!
   !! Checks: (1) subgrid eustatic == independent volume-balance root (bisection);
   !! (2) subgrid conserves ocean mass to machine precision while the coastline
   !! migrates; (3) the binary coastline gives a DIFFERENT dphi here (the term is
   !! active) and does NOT satisfy the volume balance; (4) with no migration
   !! (coast does not move) subgrid and binary agree bit-for-bit.
   use fe_precision,  only: wp
   use fe_constants,  only: rho_ice, rho_water, pi
   use fe_response,   only: response, response_init_elastic, response_init_ve, response_init_null
   use fe_sht,        only: sht_grid, sht_grid_init, sht_grid_surface_integral, sht_grid_destroy
   use fe_sle,        only: sle_solve, sle_solver, sle_result
   use fe_field,      only: spherical_cap, exp_basin
   implicit none

   real(wp), parameter :: DEG = pi/180.0_wp
   integer,  parameter :: LMAX = 48
   type(sht_grid)      :: sht
   type(response) :: resp
   type(sle_solver)    :: sle
   type(sle_result)    :: res
   real(wp), allocatable :: topo0(:,:), d_ice(:,:), ice(:,:), rsl(:,:), C(:,:)
   real(wp) :: ice_int, h_root, dphi_sub, dphi_bin, rho_ratio
   logical  :: ok

   ok = .true.
   call sht_grid_init(sht, LMAX, nlat=2*LMAX, nphi=4*LMAX)
   allocate(topo0(sht%nphi,sht%nlat), d_ice(sht%nphi,sht%nlat), &
            ice(sht%nphi,sht%nlat), rsl(sht%nphi,sht%nlat), C(sht%nphi,sht%nlat))
   rho_ratio = rho_ice/rho_water

   ! A broad, gently sloping ocean basin (centre below sea level, sloping up to a
   ! land rim) so that a sea-level change sweeps a non-trivial coastal area.
   call exp_basin(sht, 90.0_wp*DEG, 0.0_wp*DEG, 500.0_wp, 2500.0_wp, 30.0_wp*DEG, topo0)
   ! A thick ice cap on the far side (on land, away from the basin) that grows, so
   ! ocean water is locked up and sea level falls -- the coastline retreats over the
   ! shallow basin slope, where binary and subgrid disagree.
   call spherical_cap(sht, 20.0_wp*DEG, 180.0_wp*DEG, 12.0_wp*DEG, 3000.0_wp, ice)
   d_ice = ice                              ! ice-free reference => change = absolute

   sle%n_inner = 1                          ! null response: no water feedback, 1 inner is exact
   sle%n_outer = 60                         ! let the migrating coastline fully converge

   ! melt source (negative: ice grows => water removed). Ice sits on land so the
   ! grounded mask (1-C) = 1 there; ice_int = -(rho_i/rho_w) integral d_ice.
   ice_int = -rho_ratio*sht_grid_surface_integral(sht, d_ice)

   write(*,'(a)') ' --- subgrid sloping-coast ocean load (null response) ---'

   ! (1)+(2) subgrid eustatic vs independent volume-balance root
   sle%subgrid = .true.;  sle%fixed_ocean = .false.
   call sle_solve(sle, sht, resp, d_ice, ice, topo0, rsl, C, res)
   dphi_sub = res%esl
   h_root   = volume_root(ice_int)
   write(*,'(a,f12.6,a,f12.6,a,es10.2)') ' (1) subgrid dphi =', dphi_sub, &
        '   volume-root =', h_root, '   diff =', abs(dphi_sub - h_root)
   write(*,'(a,es10.2,a,f8.4,a,i3)') '     mass resid =', res%mass_resid, &
        '   ocean frac =', res%ocean_frac, '   outer iters =', res%n_outer_done
   if (abs(dphi_sub - h_root) > 1.0e-4_wp) then
      write(*,'(a)') '     FAIL: subgrid eustatic != volume-balance root';  ok = .false.
   end if
   if (res%mass_resid > 1.0e-12_wp) then
      write(*,'(a)') '     FAIL: subgrid does not conserve ocean mass';  ok = .false.
   end if

   ! (3) binary coastline: a different dphi here, and it does NOT solve the volume
   ! balance (it ignores the bathymetry swept by the moving coast).
   sle%subgrid = .false.
   call sle_solve(sle, sht, resp, d_ice, ice, topo0, rsl, C, res)
   dphi_bin = res%esl
   write(*,'(a,f12.6,a,f12.6)') ' (3) binary  dphi =', dphi_bin, &
        '   (subgrid - binary) =', dphi_sub - dphi_bin
   if (abs(dphi_bin - dphi_sub) < 1.0e-2_wp) then
      write(*,'(a)') '     FAIL: subgrid term inactive (coast not migrating here?)';  ok = .false.
   end if
   if (abs(volume_excess(dphi_bin, ice_int)) < abs(volume_excess(dphi_sub, ice_int))) then
      write(*,'(a)') '     FAIL: binary fits the volume balance better than subgrid';  ok = .false.
   end if

   ! (4) no migration => subgrid and binary agree bit-for-bit. A tiny load keeps the
   ! coastline fixed (the sea-level change does not cross any grid point's topo0).
   call no_migration_check(ok)

   write(*,'(a)') ''
   if (ok) then
      write(*,'(a)') ' PASS: subgrid ocean load matches the sloping-coast volume balance'
   else
      write(*,'(a)') ' FAIL: subgrid ocean-load validation did not all pass'
      call sht_grid_destroy(sht);  error stop 1
   end if
   call sht_grid_destroy(sht)

contains

   real(wp) function volume_excess(h, target) result(e)
      !! rho_w-scaled water-volume change for a uniform sea-level shift h on the
      !! sloping basin, minus the target (the melt source). Zero at the balanced h.
      !! Water column change per point = max(0,h-topo0) - max(0,-topo0) (the cell
      !! fills from its bed once submerged; drains to its bed once emerged).
      real(wp), intent(in) :: h, target
      real(wp), allocatable :: col(:,:)
      allocate(col(sht%nphi,sht%nlat))
      col = max(0.0_wp, h - topo0) - max(0.0_wp, -topo0)
      e = sht_grid_surface_integral(sht, col) - target
   end function volume_excess

   real(wp) function volume_root(target) result(h)
      !! Bisect the monotone volume_excess(h) for the mass-conserving sea level.
      real(wp), intent(in) :: target
      real(wp) :: a, b, fa, fm, m
      integer  :: it
      a = -4000.0_wp;  b = 0.0_wp           ! sea level falls (ice grows): h in [-4000,0]
      fa = volume_excess(a, target)
      do it = 1, 200
         m = 0.5_wp*(a + b);  fm = volume_excess(m, target)
         if (fa*fm <= 0.0_wp) then;  b = m;  else;  a = m;  fa = fm;  end if
      end do
      h = 0.5_wp*(a + b)
   end function volume_root

   subroutine no_migration_check(ok)
      !! With a load too small to move the coastline, C == C0 everywhere, so the
      !! subgrid term vanishes and subgrid == binary to machine precision.
      logical, intent(inout) :: ok
      real(wp) :: dsub, dbin
      ice = 0.0_wp
      call spherical_cap(sht, 20.0_wp*DEG, 180.0_wp*DEG, 12.0_wp*DEG, 1.0e-3_wp, ice)
      d_ice = ice
      sle%subgrid = .true.
      call sle_solve(sle, sht, resp, d_ice, ice, topo0, rsl, C, res);  dsub = res%esl
      sle%subgrid = .false.
      call sle_solve(sle, sht, resp, d_ice, ice, topo0, rsl, C, res);  dbin = res%esl
      write(*,'(a,es12.4)') ' (4) no-migration |subgrid - binary| dphi =', abs(dsub - dbin)
      if (abs(dsub - dbin) > 1.0e-10_wp) then
         write(*,'(a)') '     FAIL: subgrid != binary when the coast does not move';  ok = .false.
      end if
   end subroutine no_migration_check

end program test_sle_subgrid
