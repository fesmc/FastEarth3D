program test_integrals
   !! Verify the closed-form element integrals (fe_radial_integrals, Martinec 2000
   !! Appendix C) against 5-point Gauss-Legendre numerical quadrature of their
   !! defining integrands. 5-point GL is exact for the polynomial integrands
   !! (I1-I6, K1-K3) and highly accurate for the smooth 1/r integrand (I7).
   use fe_precision, only: wp
   use fe_radial_integrals
   implicit none

   ! A representative element away from r=0 (so I7's 1/r is non-singular).
   real(wp), parameter :: rk = 5.00e6_wp, rk1 = 5.04e6_wp
   real(wp), parameter :: tol = 1.0e-7_wp
   logical :: ok
   real(wp) :: worst

   ok = .true.; worst = 0.0_wp

   call chk_mat('I1', elem_i1(rk,rk1), num_mat(1,1,2))
   call chk_mat('I2', elem_i2(rk,rk1), num_mat(0,1,2))
   call chk_mat('I3', elem_i3(rk,rk1), num_mat(0,1,1))
   call chk_mat('I4', elem_i4(rk,rk1), num_mat(0,0,2))
   call chk_mat('I5', elem_i5(rk,rk1), num_mat(0,0,1))
   call chk_mat('I6', elem_i6(rk,rk1), num_mat(0,0,0))
   call chk_mat('I7', elem_i7(rk,rk1), num_mat(0,0,-1))
   call chk_vec('K1', elem_k1(rk,rk1), num_vec(1,2))
   call chk_vec('K2', elem_k2(rk,rk1), num_vec(0,1))
   call chk_vec('K3', elem_k3(rk,rk1), num_vec(0,2))

   print '(a,es12.4)', ' worst relative error = ', worst
   if (ok) then
      print '(a)', ' PASS: Appendix C element integrals match quadrature'
   else
      print '(a)', ' FAIL: an element integral disagrees with quadrature'
      error stop 1
   end if

contains

   real(wp) function psi(a, r) result(p)
      integer,  intent(in) :: a
      real(wp), intent(in) :: r
      if (a == 1) then
         p = (rk1 - r)/(rk1 - rk)
      else
         p = (r - rk)/(rk1 - rk)
      end if
   end function psi

   real(wp) function dpsi(a) result(p)
      integer, intent(in) :: a
      if (a == 1) then
         p = -1.0_wp/(rk1 - rk)
      else
         p =  1.0_wp/(rk1 - rk)
      end if
   end function dpsi

   function num_mat(da, db, p) result(m)
      !! 2x2 matrix of ∫ ψ_a^(da) ψ_b^(db) r^p dr by 5-point Gauss-Legendre.
      integer, intent(in) :: da, db, p
      real(wp) :: m(2,2)
      real(wp) :: x(5), w(5), rm, rh, r, wq, fa, fb
      integer  :: a, b, q
      call gl5(x, w)
      rm = 0.5_wp*(rk + rk1); rh = 0.5_wp*(rk1 - rk)
      m = 0.0_wp
      do a = 1, 2
         do b = 1, 2
            do q = 1, 5
               r  = rm + rh*x(q)
               wq = rh*w(q)
               fa = merge(dpsi(a), psi(a, r), da == 1)
               fb = merge(dpsi(b), psi(b, r), db == 1)
               m(a,b) = m(a,b) + wq*fa*fb*r**p
            end do
         end do
      end do
   end function num_mat

   function num_vec(da, p) result(v)
      !! 2-vector of ∫ ψ_a^(da) r^p dr by 5-point Gauss-Legendre.
      integer, intent(in) :: da, p
      real(wp) :: v(2)
      real(wp) :: x(5), w(5), rm, rh, r, wq, fa
      integer  :: a, q
      call gl5(x, w)
      rm = 0.5_wp*(rk + rk1); rh = 0.5_wp*(rk1 - rk)
      v = 0.0_wp
      do a = 1, 2
         do q = 1, 5
            r  = rm + rh*x(q)
            wq = rh*w(q)
            fa = merge(dpsi(a), psi(a, r), da == 1)
            v(a) = v(a) + wq*fa*r**p
         end do
      end do
   end function num_vec

   subroutine gl5(x, w)
      real(wp), intent(out) :: x(5), w(5)
      x = [-0.9061798459386640_wp, -0.5384693101056831_wp, 0.0_wp, &
            0.5384693101056831_wp,  0.9061798459386640_wp]
      w = [ 0.2369268850561891_wp,  0.4786286704993665_wp, &
            0.5688888888888889_wp,  0.4786286704993665_wp, &
            0.2369268850561891_wp]
   end subroutine gl5

   subroutine chk_mat(name, a, b)
      character(*), intent(in) :: name
      real(wp),     intent(in) :: a(2,2), b(2,2)
      real(wp) :: rel
      rel = maxval(abs(a - b))/max(maxval(abs(b)), tiny(1.0_wp))
      worst = max(worst, rel)
      if (rel >= tol) then
         print '(a,a,a,es12.4)', ' FAIL ', name, '  rel err = ', rel
         ok = .false.
      end if
   end subroutine chk_mat

   subroutine chk_vec(name, a, b)
      character(*), intent(in) :: name
      real(wp),     intent(in) :: a(2), b(2)
      real(wp) :: rel
      rel = maxval(abs(a - b))/max(maxval(abs(b)), tiny(1.0_wp))
      worst = max(worst, rel)
      if (rel >= tol) then
         print '(a,a,a,es12.4)', ' FAIL ', name, '  rel err = ', rel
         ok = .false.
      end if
   end subroutine chk_vec

end program test_integrals
