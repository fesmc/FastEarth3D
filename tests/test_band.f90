program test_band
   !! Validate fe_band (pivoted banded LU) against a dense Gaussian-elimination
   !! reference, including a saddle-point-style matrix with ZERO diagonal entries
   !! (which forces pivoting — the case the radial operator's pressure block hits).
   use fe_precision, only: wp
   use fe_band,      only: band_lu, band_build, band_solve, band_destroy
   implicit none
   logical :: ok
   ok = .true.

   call case_random_band(50, 4, 3, ok)
   call case_zero_diagonal(ok)
   call case_saddle(ok)

   write(*,'(a)') ''
   if (ok) then
      write(*,'(a)') ' PASS: banded LU matches dense solve (incl. pivoting on zero diagonals)'
   else
      write(*,'(a)') ' FAIL: banded LU validation did not all pass'
      error stop 1
   end if

contains

   subroutine case_random_band(n, kl, ku, ok)
      integer, intent(in) :: n, kl, ku
      logical, intent(inout) :: ok
      real(wp) :: A(n,n), b(n), x(n), xref(n), r
      integer  :: i, j, nz, t
      integer,  allocatable :: rows(:), cols(:)
      real(wp), allocatable :: vals(:)
      type(band_lu) :: blu
      logical :: okf
      ! deterministic pseudo-random banded matrix, diagonally dominant
      A = 0.0_wp
      do j = 1, n
         do i = max(1,j-ku), min(n,j+kl)
            A(i,j) = real(mod(7*i+3*j, 11) - 5, wp)
         end do
         A(j,j) = A(j,j) + real(2*(kl+ku)+5, wp)      ! dominant diagonal
      end do
      do i = 1, n;  b(i) = real(mod(13*i, 17) - 8, wp);  end do
      ! COO
      nz = 0
      do j = 1, n; do i = max(1,j-ku), min(n,j+kl); nz = nz + 1; end do; end do
      allocate(rows(nz), cols(nz), vals(nz));  t = 0
      do j = 1, n; do i = max(1,j-ku), min(n,j+kl)
         t = t+1; rows(t)=i; cols(t)=j; vals(t)=A(i,j)
      end do; end do
      call band_build(blu, n, nz, rows, cols, vals, okf)
      call band_solve(blu, b, x)
      call dense_solve(n, A, b, xref)
      r = maxval(abs(x - xref)) / max(maxval(abs(xref)), 1.0_wp)
      write(*,'(a,i0,a,i0,a,i0,a,es10.2)') ' (1) random band n=',n,' kl=',kl,' ku=',ku, &
           '  rel err = ', r
      if (.not. okf .or. r > 1.0e-12_wp) then;  write(*,'(a)') '     FAIL';  ok=.false.;  end if
      call band_destroy(blu);  deallocate(rows,cols,vals)
   end subroutine case_random_band

   subroutine case_zero_diagonal(ok)
      !! A 4x4 banded matrix with a zero on the diagonal (needs a row swap).
      logical, intent(inout) :: ok
      integer, parameter :: n=4
      real(wp) :: A(n,n), b(n), x(n), xref(n), r
      integer :: rows(10), cols(10), i, j, nz
      real(wp) :: vals(10)
      type(band_lu) :: blu;  logical :: okf
      A = 0.0_wp
      A(1,:) = [0.0_wp, 2.0_wp, 0.0_wp, 0.0_wp]   ! zero (1,1) -> must pivot
      A(2,:) = [3.0_wp, 1.0_wp, 4.0_wp, 0.0_wp]
      A(3,:) = [0.0_wp, 5.0_wp, 1.0_wp, 6.0_wp]
      A(4,:) = [0.0_wp, 0.0_wp, 7.0_wp, 2.0_wp]
      b = [1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp]
      nz = 0
      do j=1,n; do i=1,n; if (A(i,j)/=0.0_wp) then; nz=nz+1; rows(nz)=i; cols(nz)=j; vals(nz)=A(i,j); end if; end do; end do
      call band_build(blu, n, nz, rows(1:nz), cols(1:nz), vals(1:nz), okf)
      call band_solve(blu, b, x)
      call dense_solve(n, A, b, xref)
      r = maxval(abs(x - xref)) / max(maxval(abs(xref)), 1.0_wp)
      write(*,'(a,es10.2)') ' (2) zero-diagonal (pivot required)  rel err = ', r
      if (.not. okf .or. r > 1.0e-12_wp) then;  write(*,'(a)') '     FAIL';  ok=.false.;  end if
      call band_destroy(blu)
   end subroutine case_zero_diagonal

   subroutine case_saddle(ok)
      !! Small saddle point [K B; B^T 0] with zero (2,2) block, banded-ish.
      logical, intent(inout) :: ok
      integer, parameter :: n=6
      real(wp) :: A(n,n), b(n), x(n), xref(n), r
      integer :: rows(36), cols(36), i, j, nz
      real(wp) :: vals(36)
      type(band_lu) :: blu;  logical :: okf
      A = 0.0_wp
      do i=1,4; A(i,i)=4.0_wp; end do
      A(1,2)=1; A(2,1)=1; A(3,4)=1; A(4,3)=1
      ! couplings to the two "pressure" dofs 5,6 (zero diagonal there)
      A(1,5)=1; A(5,1)=1; A(2,5)=1; A(5,2)=1
      A(3,6)=1; A(6,3)=1; A(4,6)=1; A(6,4)=1
      do i=1,n; b(i)=real(i,wp); end do
      nz=0
      do j=1,n; do i=1,n; if (A(i,j)/=0.0_wp) then; nz=nz+1; rows(nz)=i; cols(nz)=j; vals(nz)=A(i,j); end if; end do; end do
      call band_build(blu, n, nz, rows(1:nz), cols(1:nz), vals(1:nz), okf)
      call band_solve(blu, b, x)
      call dense_solve(n, A, b, xref)
      r = maxval(abs(x - xref)) / max(maxval(abs(xref)), 1.0_wp)
      write(*,'(a,es10.2)') ' (3) saddle point (zero 2,2 block)   rel err = ', r
      if (.not. okf .or. r > 1.0e-12_wp) then;  write(*,'(a)') '     FAIL';  ok=.false.;  end if
      call band_destroy(blu)
   end subroutine case_saddle

   subroutine dense_solve(n, Ain, bin, x)
      !! Reference dense Gaussian elimination with partial pivoting.
      integer,  intent(in)  :: n
      real(wp), intent(in)  :: Ain(n,n), bin(n)
      real(wp), intent(out) :: x(n)
      real(wp) :: A(n,n), b(n), f, p
      integer  :: i, j, k, kp
      A = Ain;  b = bin
      do k = 1, n-1
         kp = k;  p = abs(A(k,k))
         do i = k+1, n;  if (abs(A(i,k))>p) then; p=abs(A(i,k)); kp=i; end if;  end do
         if (kp/=k) then;  A([k,kp],:)=A([kp,k],:);  f=b(k);b(k)=b(kp);b(kp)=f;  end if
         do i = k+1, n
            f = A(i,k)/A(k,k)
            A(i,k:n) = A(i,k:n) - f*A(k,k:n)
            b(i) = b(i) - f*b(k)
         end do
      end do
      do i = n, 1, -1
         x(i) = (b(i) - dot_product(A(i,i+1:n), x(i+1:n))) / A(i,i)
      end do
   end subroutine dense_solve

end program test_band
