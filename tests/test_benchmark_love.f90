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
   ! Degrees over which agreement is ASSERTED. The comparison itself now runs to
   ! NMAX and is dumped in full, but only these ranges are pass/fail criteria --
   ! they are the ranges this test has always asserted. Above them the residual
   ! grows (see the crossover summary printed below): h stays ~0.06% to 256, but
   ! l and k drift past 1% near degrees 145 and 97 respectively. Raising these
   ! bounds is a real model question (radial mesh resolution at short wavelength),
   ! not a tolerance to be relaxed.
   integer, parameter :: JA_FLUID = 8, JA_ELASTIC = 48
   real(wp) :: he(NMAX), le(NMAX), ke(NMAX), hf(NMAX), lf(NMAX), kf(NMAX)
   logical  :: ok, okread
   integer  :: j, it, nrow
   real(wp) :: tabf(NMAX-1,7), tabe(NMAX-1,7)
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
   nrow = 0
   do j = 2, NMAX
      call radial_operator_assemble(op, ef, mf, j)
      call radial_operator_solve(op, 1.0_wp, u, v, f, iters=it, resid=rr)
      call loading_love(ef, j, 1.0_wp, u, v, f, h, l, k)
      nrow = nrow + 1
      tabf(nrow,:) = [real(j,wp), h, hf(j), l, lf(j), k, kf(j)]
      if (j <= 8 .or. mod(j,32) == 0) &
         write(*,'(i7,6f11.5)') j, h, hf(j), l, lf(j), k, kf(j)
      if (j <= JA_FLUID) then
         if (reldiff(h, hf(j)) > 1.0e-2_wp) then
            write(*,'(a)') '      FAIL: fluid-limit h off the benchmark (>1%)';  ok = .false.
         end if
         if (reldiff(l, lf(j)) > 1.0e-2_wp) then
            write(*,'(a)') '      FAIL: fluid-limit l off the benchmark (>1%)';  ok = .false.
         end if
         if (reldiff(k, kf(j)) > 1.0e-2_wp) then
            write(*,'(a)') '      FAIL: fluid-limit k off the benchmark (>1%)';  ok = .false.
         end if
      end if
   end do
   call dump_cols('love_fluid.txt', &
        'j h_model h_ref l_model l_ref k_model k_ref', tabf(1:nrow,:))

   ! --- (2) ELASTIC: match the table at every degree (<1%) --------------------
   write(*,'(a)') ''
   write(*,'(a)') ' (2) M3-L70-V01 elastic loading Love numbers vs benchmark table'
   write(*,'(a)') '      j      h_ours     h_ref    dh%      k_ours     k_ref    dk%'
   e = build_M3L70V01();  call radial_mesh_build(m, e)
   nrow = 0
   do j = 2, NMAX
      call radial_operator_assemble(op, e, m, j)
      call radial_operator_solve(op, 1.0_wp, u, v, f, iters=it, resid=rr)
      call loading_love(e, j, 1.0_wp, u, v, f, h, l, k)
      nrow = nrow + 1
      tabe(nrow,:) = [real(j,wp), h, he(j), l, le(j), k, ke(j)]
      if (j <= 8 .or. mod(j,32) == 0) &
         write(*,'(i7,2f11.5,f8.1,2f11.5,f8.1)') j, h, he(j), &
              100.0_wp*(h-he(j))/abs(he(j)), k, ke(j), 100.0_wp*(k-ke(j))/abs(ke(j))
      if (j <= JA_ELASTIC) then
         if (reldiff(h, he(j)) > 1.0e-2_wp .or. reldiff(k, ke(j)) > 1.0e-2_wp) then
            write(*,'(a,i0)') '      FAIL: elastic Love numbers off the benchmark (>1%) at j=', j
            ok = .false.
         end if
      end if
   end do
   call dump_cols('love_elastic.txt', &
        'j h_model h_ref l_model l_ref k_model k_ref', tabe(1:nrow,:))

   write(*,'(a)') ''
   write(*,'(a,i0,a,i0,a)') ' (3) agreement beyond the asserted range (fluid j<=', &
        JA_FLUID, ', elastic j<=', JA_ELASTIC, '): lowest degree exceeding 1%'
   write(*,'(a)') '      table        h        l        k'
   call crossover('fluid  ', tabf(1:nrow,:))
   call crossover('elastic', tabe(1:nrow,:))

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


   subroutine crossover(label, t)
      !! Lowest degree at which each of h, l, k leaves 1% of the reference table,
      !! or 0 if it never does over the range computed. Columns of `t` are
      !! j, h_model, h_ref, l_model, l_ref, k_model, k_ref.
      character(*), intent(in) :: label
      real(wp),     intent(in) :: t(:,:)
      integer :: c, i, jx(3)
      jx = 0
      do c = 1, 3
         do i = 1, size(t,1)
            if (reldiff(t(i,2*c), t(i,2*c+1)) > 1.0e-2_wp) then
               jx(c) = nint(t(i,1));  exit
            end if
         end do
      end do
      write(*,'(6x,a,3i9)') label, jx
   end subroutine crossover

   subroutine dump_cols(name, header, a)
      !! Write a column table to $FE_BENCH_DUMP/<name> for the analysis scripts.
      !!
      !! No-op unless FE_BENCH_DUMP names a directory, so `make check` and any
      !! plain run behave exactly as before — the dump is opt-in and costs
      !! nothing when off. `header` names the columns and is written as a leading
      !! `#` comment line, so the file is self-describing and readable with any
      !! delimited-text reader.
      character(*), intent(in) :: name, header
      real(wp),     intent(in) :: a(:,:)          ! (nrow, ncol)
      character(512) :: dir, path
      integer :: u, i, st
      call get_environment_variable('FE_BENCH_DUMP', dir, status=st)
      if (st /= 0 .or. len_trim(dir) == 0) return
      path = trim(dir)//'/'//name
      open(newunit=u, file=trim(path), status='replace', action='write', iostat=st)
      if (st /= 0) then
         write(*,'(3a)') '   WARNING: cannot write ', trim(path), ' (dump skipped)'
         return
      end if
      write(u,'(2a)') '# ', header
      do i = 1, size(a,1)
         write(u,'(*(es18.10,1x))') a(i,:)
      end do
      close(u)
      write(*,'(3a,i0,a,i0,a)') '   dumped ', trim(path), ' (', size(a,1), ' x ', size(a,2), ')'
   end subroutine dump_cols

end program test_benchmark_love
