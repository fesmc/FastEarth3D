module fe_band
   !! Banded LU with partial pivoting — a small, dependency-free, thread-safe direct
   !! solver for the per-degree radial operator. The operator is a saddle-point
   !! system with a ZERO pressure (Π,Π) block, so a factorization MUST pivot; it is
   !! banded (half-bandwidth ~6 for the node-interleaved P1/P0 layout), so a band LU
   !! factors in O(n·kl²) and solves in O(n·(kl+ku)) — far cheaper and far lighter on
   !! cache than GMRES+ILU, and (unlike LIS) re-entrant, so many degrees can be
   !! solved concurrently from different threads.
   !!
   !! Storage and algorithm follow LAPACK's unblocked dgbtf2 / dgbtrs exactly:
   !! A(i,j) lives at ab(kv+1+i-j, j) with kv = kl+ku, ldab = 2*kl+ku+1 (the top kl
   !! rows are fill-in workspace created by pivoting).
   use fe_precision, only: wp
   implicit none
   private
   public :: band_lu
   public :: band_build, band_solve, band_destroy

   type :: band_lu
      integer :: n = 0, kl = 0, ku = 0, kv = 0, ldab = 0
      real(wp), allocatable :: ab(:,:)     !! (ldab, n) factored band storage
      integer,  allocatable :: ipiv(:)     !! (n) row pivots
      logical :: ready = .false.
   end type band_lu

contains

   subroutine band_build(self, n, nz, rows, cols, vals, ok)
      !! Build band storage from a COO matrix (entries may repeat — summed) and
      !! LU-factor it in place with partial pivoting.
      type(band_lu), intent(inout) :: self
      integer,  intent(in)  :: n, nz
      integer,  intent(in)  :: rows(:), cols(:)
      real(wp), intent(in)  :: vals(:)
      logical,  intent(out) :: ok
      integer :: t, i, j, kl, ku

      call band_destroy(self)
      ! bandwidths from the sparsity pattern
      kl = 0;  ku = 0
      do t = 1, nz
         if (vals(t) == 0.0_wp) cycle
         i = rows(t);  j = cols(t)
         kl = max(kl, i - j)
         ku = max(ku, j - i)
      end do
      self%n = n;  self%kl = kl;  self%ku = ku;  self%kv = kl + ku
      self%ldab = 2*kl + ku + 1
      allocate(self%ab(self%ldab, n), self%ipiv(n))
      self%ab = 0.0_wp
      ! scatter A(i,j) -> ab(kv+1+i-j, j)
      do t = 1, nz
         if (vals(t) == 0.0_wp) cycle
         i = rows(t);  j = cols(t)
         self%ab(self%kv + 1 + i - j, j) = self%ab(self%kv + 1 + i - j, j) + vals(t)
      end do
      call band_factor(self%n, self%kl, self%ku, self%kv, self%ldab, self%ab, self%ipiv, ok)
      self%ready = ok
   end subroutine band_build

   subroutine band_factor(n, kl, ku, kv, ldab, ab, ipiv, ok)
      !! Unblocked banded LU with partial pivoting (LAPACK dgbtf2).
      integer,  intent(in)    :: n, kl, ku, kv, ldab
      real(wp), intent(inout) :: ab(ldab, n)
      integer,  intent(out)   :: ipiv(n)
      logical,  intent(out)   :: ok
      integer  :: j, jp, ju, km, i, jj, p
      real(wp) :: piv, tmp
      ok = .true.
      ju = 1
      do j = 1, n
         ! zero the fill-in column that pivoting may create
         if (j + kv <= n) ab(1:kl, j+kv) = 0.0_wp
         km = min(kl, n - j)                 ! # subdiagonal entries in column j
         ! pivot = arg max |ab(kv+1 .. kv+1+km, j)|
         jp = 0;  piv = abs(ab(kv+1, j))
         do p = 1, km
            if (abs(ab(kv+1+p, j)) > piv) then;  piv = abs(ab(kv+1+p, j));  jp = p;  end if
         end do
         ipiv(j) = j + jp
         if (ab(kv+1+jp, j) == 0.0_wp) then;  ok = .false.;  return;  end if
         ju = max(ju, min(j + ku + jp, n))
         ! swap rows j and j+jp across columns j..ju (band storage)
         if (jp /= 0) then
            do jj = j, ju
               tmp                  = ab(kv+1+(j)-jj, jj)
               ab(kv+1+(j)-jj, jj)  = ab(kv+1+(j+jp)-jj, jj)
               ab(kv+1+(j+jp)-jj, jj) = tmp
            end do
         end if
         if (km > 0) then
            ! scale subdiagonal (the L multipliers) by 1/pivot
            piv = 1.0_wp / ab(kv+1, j)
            do p = 1, km
               ab(kv+1+p, j) = ab(kv+1+p, j) * piv
            end do
            ! rank-1 update of the trailing band: A(i,jj) -= L(i,j)*U(j,jj)
            do jj = j+1, ju
               tmp = ab(kv+1+(j)-jj, jj)        ! U(j,jj)
               if (tmp /= 0.0_wp) then
                  do p = 1, km
                     i = j + p
                     ab(kv+1+i-jj, jj) = ab(kv+1+i-jj, jj) - ab(kv+1+p, j)*tmp
                  end do
               end if
            end do
         end if
      end do
   end subroutine band_factor

   subroutine band_solve(self, b, x)
      !! Solve A x = b for the factored system (LAPACK dgbtrs, no transpose).
      type(band_lu), intent(in)  :: self
      real(wp),       intent(in)  :: b(:)
      real(wp),       intent(out) :: x(:)
      integer  :: j, lm, l, i, kv
      real(wp) :: tmp
      kv = self%kv
      x = b
      ! forward: P and L  (L y = P b)
      do j = 1, self%n - 1
         lm = min(self%kl, self%n - j)
         l  = self%ipiv(j)
         if (l /= j) then;  tmp = x(l);  x(l) = x(j);  x(j) = tmp;  end if
         do i = 1, lm
            x(j+i) = x(j+i) - self%ab(kv+1+i, j)*x(j)
         end do
      end do
      ! back: U x = y
      do j = self%n, 1, -1
         x(j) = x(j) / self%ab(kv+1, j)
         lm = min(kv, j-1)
         do i = 1, lm
            x(j-i) = x(j-i) - self%ab(kv+1-i, j)*x(j)
         end do
      end do
   end subroutine band_solve

   subroutine band_destroy(self)
      type(band_lu), intent(inout) :: self
      if (allocated(self%ab))   deallocate(self%ab)
      if (allocated(self%ipiv)) deallocate(self%ipiv)
      self%n = 0;  self%kl = 0;  self%ku = 0;  self%kv = 0;  self%ldab = 0
      self%ready = .false.
   end subroutine band_destroy

end module fe_band
