program test_mesh
   !! Validate the radial mesh over the whole sphere of M3-L70-V01: it must span
   !! [0, r_earth] exactly (Martinec 2000 meshes ⟨0,a⟩ through the centre), be
   !! strictly ascending, place a node on every material interface — including
   !! the CMB (no element straddles a discontinuity) — and honour the element-size
   !! targets (5/10/40 km by depth).
   use fe_precision,       only: wp
   use fe_earth_structure, only: earth_model, build_M3L70V01
   use fe_radial_fe,       only: radial_mesh_build, radial_mesh
   implicit none

   type(earth_model) :: em
   type(radial_mesh) :: mesh
   integer  :: i, k
   logical  :: ok, found
   real(wp) :: dr, depth_mid, rmid, target, hit
   real(wp), parameter :: km = 1.0e3_wp
   real(wp), parameter :: tolnode = 1.0_wp     ! 1 m
   ! Interior interfaces that must be nodes [m] (incl. the CMB at 3480 km).
   real(wp), parameter :: interfaces(4) = [6301.0_wp, 5951.0_wp, 5701.0_wp, &
                                           3480.0_wp]*km

   ok = .true.
   em = build_M3L70V01()
   call radial_mesh_build(mesh, em)
   print '(a,i0,a,i0,a)', ' mesh: ', mesh%nr, ' nodes, ', mesh%ne, ' elements'
   print '(a,f8.1,a,f8.1,a)', ' span: ', mesh%r(1)/km, ' km -> ', &
        mesh%r(mesh%nr)/km, ' km'

   ! --- Node count consistency and bounds -------------------------------------
   if (mesh%nr /= mesh%ne + 1) ok = .false.
   if (abs(mesh%r(1) - 0.0_wp)  > tolnode) ok = .false.
   if (abs(mesh%r(mesh%nr) - em%r_earth) > tolnode) ok = .false.

   ! --- Strictly ascending ----------------------------------------------------
   do i = 2, mesh%nr
      if (mesh%r(i) <= mesh%r(i-1)) ok = .false.
   end do

   ! --- Every interior interface is a node ------------------------------------
   do k = 1, size(interfaces)
      found = .false.
      do i = 1, mesh%nr
         if (abs(mesh%r(i) - interfaces(k)) <= tolnode) found = .true.
      end do
      if (.not. found) then
         print '(a,f8.1,a)', ' MISSING interface node at ', interfaces(k)/km, ' km'
         ok = .false.
      end if
   end do

   ! --- Element sizes within the depth-dependent target -----------------------
   hit = 0.0_wp
   do i = 1, mesh%ne
      dr        = mesh%r(i+1) - mesh%r(i)
      rmid      = 0.5_wp*(mesh%r(i+1) + mesh%r(i))
      depth_mid = em%r_earth - rmid
      if (depth_mid <= 70.0_wp*km) then
         target = 5.0_wp*km
      else if (depth_mid <= 670.0_wp*km) then
         target = 10.0_wp*km
      else
         target = 40.0_wp*km
      end if
      if (dr > target + tolnode) then
         print '(a,i0,a,f7.2,a,f7.2,a)', ' element ', i, ' size ', dr/km, &
              ' km exceeds target ', target/km, ' km'
         ok = .false.
      end if
      ! Check the element's layer assignment contains its midpoint.
      if (rmid < em%layers(mesh%elem_layer(i))%r_bot - tolnode .or. &
          rmid > em%layers(mesh%elem_layer(i))%r_top + tolnode) ok = .false.
      hit = max(hit, dr)
   end do
   print '(a,f7.2,a)', ' max element size = ', hit/km, ' km'

   if (ok) then
      print '(a)', ' PASS: radial mesh spans the shell, captures interfaces, honours spacing'
   else
      print '(a)', ' FAIL: radial mesh consistency check'
      error stop 1
   end if
end program test_mesh
