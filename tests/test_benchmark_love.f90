program test_benchmark_love
   !! Rung-2/3 validation against the GIA-benchmark normal-mode loading Love
   !! table data/benchmarks/love_M3-L70-V01/mod_M3-L70-V01 (see
   !! data/benchmarks/PROVENANCE.md). This table is the authoritative M3-L70-V01
   !! reference (TABOO/ALMA normal-mode output; independently reproduced with
   !! TABOO NV=3/CODE=7). It carries, per degree, the ELASTIC (t=0) and FLUID
   !! (t->inf) loading Love numbers h, l, k.
   !!
   !! (1) FLUID limit. The t->inf relaxed state of the Maxwell model equals the
   !!     ELASTIC solve of the model with every Maxwell layer fluidized (mu=0).
   !!     We assert this matches the table's fluid columns to <1% at all degrees.
   !!     This is a strong, layered check of self-gravity, the R_k interface
   !!     buoyancy, the inviscid core, incompressibility AND the l normalization
   !!     (so it pins l's sign/scale — previously the one open Love convention).
   !!
   !! (2) ELASTIC. Our elastic loading Love numbers match the table to <1% at
   !!     every degree 2..48 (residual is P1 mesh discretization, same size as the
   !!     fluid-limit residual). This was the long-standing low-degree discrepancy
   !!     (~50% too soft at j=2): a single transposed index in the self-gravity
   !!     potential-gradient force (eq 65/81, the U-F coupling i2) — now fixed.
   !!     See doc/formulation.md "Elastic low-degree discrepancy (FIXED)".
   use fe_precision,       only: wp
   use fe_earth_structure, only: earth_model, build_M3L70V01, RHEOL_FLUID
   use fe_radial_fe,       only: radial_operator_solve, radial_operator_assemble, radial_mesh_build, radial_mesh, radial_operator, loading_love, &
                                 radial_fe_finalize
   implicit none
   character(*), parameter :: REF = 'data/benchmarks/love_M3-L70-V01/mod_M3-L70-V01'
   integer, parameter :: NMAX = 256
   real(wp) :: he(NMAX), le(NMAX), ke(NMAX), hf(NMAX), lf(NMAX), kf(NMAX)
   logical  :: ok, okread
   integer  :: j, it
   real(wp) :: u, v, f, h, l, k, rr
   type(earth_model)     :: e, ef
   type(radial_mesh)     :: m, mf
   type(radial_operator) :: op

   ok = .true.
   call read_ref(REF, he, le, ke, hf, lf, kf, okread)
   if (.not. okread) then
      write(*,'(2a)') ' FAIL: cannot read benchmark table ', REF
      error stop 1
   end if

   ! --- (1) FLUID limit: fluidized-mantle M3 vs table fluid columns -----------
   write(*,'(a)') ' (1) M3-L70-V01 fluid (relaxed) limit vs benchmark table'
   write(*,'(a)') '      j      h_ours     h_ref     l_ours     l_ref     k_ours     k_ref'
   ef = build_M3L70V01()
   ef%layers(2)%mu = 0.0_wp;  ef%layers(2)%rheology = RHEOL_FLUID
   ef%layers(3)%mu = 0.0_wp;  ef%layers(3)%rheology = RHEOL_FLUID
   ef%layers(4)%mu = 0.0_wp;  ef%layers(4)%rheology = RHEOL_FLUID
   call radial_mesh_build(mf, ef)
   do j = 2, 8
      call radial_operator_assemble(op, ef, mf, j)
      call radial_operator_solve(op, 1.0_wp, u, v, f, iters=it, resid=rr)
      call loading_love(ef, j, 1.0_wp, u, v, f, h, l, k)
      write(*,'(i7,6f11.5)') j, h, hf(j), l, lf(j), k, kf(j)
      if (reldiff(h, hf(j)) > 1.0e-2_wp) then
         write(*,'(a)') '      FAIL: fluid-limit h off the benchmark (>1%)';  ok = .false.
      end if
      if (reldiff(l, lf(j)) > 1.0e-2_wp) then
         write(*,'(a)') '      FAIL: fluid-limit l off the benchmark (>1%)';  ok = .false.
      end if
      if (reldiff(k, kf(j)) > 1.0e-2_wp) then
         write(*,'(a)') '      FAIL: fluid-limit k off the benchmark (>1%)';  ok = .false.
      end if
   end do

   ! --- (2) ELASTIC: match the table at every degree (<1%) --------------------
   write(*,'(a)') ''
   write(*,'(a)') ' (2) M3-L70-V01 elastic loading Love numbers vs benchmark table'
   write(*,'(a)') '      j      h_ours     h_ref    dh%      k_ours     k_ref    dk%'
   e = build_M3L70V01();  call radial_mesh_build(m, e)
   do j = 2, 48
      call radial_operator_assemble(op, e, m, j)
      call radial_operator_solve(op, 1.0_wp, u, v, f, iters=it, resid=rr)
      call loading_love(e, j, 1.0_wp, u, v, f, h, l, k)
      if (j <= 8 .or. mod(j,8) == 0) &
         write(*,'(i7,2f11.5,f8.1,2f11.5,f8.1)') j, h, he(j), &
              100.0_wp*(h-he(j))/abs(he(j)), k, ke(j), 100.0_wp*(k-ke(j))/abs(ke(j))
      if (reldiff(h, he(j)) > 1.0e-2_wp .or. reldiff(k, ke(j)) > 1.0e-2_wp) then
         write(*,'(a,i0)') '      FAIL: elastic Love numbers off the benchmark (>1%) at j=', j
         ok = .false.
      end if
   end do

   write(*,'(a)') ''
   if (ok) then
      write(*,'(a)') ' PASS: elastic AND fluid M3-L70-V01 Love numbers match the benchmark (<1%)'
   else
      write(*,'(a)') ' FAIL: benchmark Love-number validation did not all pass'
      call radial_fe_finalize();  error stop 1
   end if
   call radial_fe_finalize()

contains

   pure real(wp) function reldiff(a, b) result(d)
      real(wp), intent(in) :: a, b
      d = 0.0_wp
      if (abs(b) > 0.0_wp) d = abs(a-b)/abs(b)
   end function reldiff

   subroutine read_ref(fname, he, le, ke, hf, lf, kf, okread)
      !! Parse the normal-mode table: 5 earth-model lines, then per degree a
      !! header "n nmodes k h l" (elastic), nmodes mode lines, and a fluid line
      !! "n nmodes k h l" (t->inf). Columns after the index are k, h, l.
      character(len=*), intent(in)  :: fname
      real(wp),         intent(out) :: he(:), le(:), ke(:), hf(:), lf(:), kf(:)
      logical,          intent(out) :: okread
      integer  :: u, i, ni, nm, mm, ios
      real(wp) :: kv, hv, lv, dum
      okread = .false.
      open(newunit=u, file=fname, status='old', action='read', iostat=ios)
      if (ios /= 0) return
      do i = 1, 5
         read(u,*,iostat=ios)                       ! skip earth-model lines
         if (ios /= 0) then;  close(u);  return;  end if
      end do
      do
         read(u,*,iostat=ios) ni, nm, kv, hv, lv     ! elastic header
         if (ios /= 0) exit
         if (ni < 1 .or. ni > size(he)) exit
         ke(ni) = kv;  he(ni) = hv;  le(ni) = lv
         do mm = 1, nm
            read(u,*,iostat=ios) dum                 ! skip one mode line (record)
            if (ios /= 0) then;  close(u);  return;  end if
         end do
         read(u,*,iostat=ios) ni, nm, kv, hv, lv     ! fluid line
         if (ios /= 0) then;  close(u);  return;  end if
         kf(ni) = kv;  hf(ni) = hv;  lf(ni) = lv
      end do
      close(u)
      okread = .true.
   end subroutine read_ref

end program test_benchmark_love
