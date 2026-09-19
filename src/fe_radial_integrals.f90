module fe_radial_integrals
   !! Closed-form definite integrals of the piecewise-linear (P1) radial base
   !! functions over a single finite element [r_k, r_{k+1}], from Martinec (2000)
   !! Appendix C (eqs C1-C10). These are the building blocks of the per-degree
   !! stiffness/coupling matrices assembled in fe_radial_fe.
   !!
   !! On the element, with h = r_{k+1} - r_k,
   !!     ψ_k(r)   = (r_{k+1} - r)/h        (local node 1)
   !!     ψ_{k+1}(r) = (r - r_k)/h          (local node 2)
   !! Each `elem_iN` returns the 2x2 element matrix M(a,b), a,b ∈ {1,2} indexing
   !! the local nodes (k, k+1); each `elem_kN` returns the 2-vector v(a).
   !!
   !!   I1: ∫ ψ'_a ψ'_b r² dr     I2: ∫ ψ_a ψ'_b r² dr     I3: ∫ ψ_a ψ'_b r dr
   !!   I4: ∫ ψ_a ψ_b r² dr       I5: ∫ ψ_a ψ_b r dr       I6: ∫ ψ_a ψ_b dr
   !!   I7: ∫ (1/r) ψ_a ψ_b dr  (singular at r=0; only used with coefficient R_k,
   !!                            which is 0 for the innermost element — eq 77)
   !!   K1: ∫ ψ'_a r² dr          K2: ∫ ψ_a r dr           K3: ∫ ψ_a r² dr
   use fe_precision, only: wp
   implicit none
   private
   public :: elem_i1, elem_i2, elem_i3, elem_i4, elem_i5, elem_i6, elem_i7
   public :: elem_k1, elem_k2, elem_k3

contains

   pure function elem_i1(rk, rk1) result(m)   ! ∫ ψ'_a ψ'_b r² dr   (C1)
      real(wp), intent(in) :: rk, rk1
      real(wp) :: m(2,2), h, d
      h = rk1 - rk
      d = (rk1*rk1 + rk1*rk + rk*rk)/(3.0_wp*h)
      m(1,1) =  d; m(1,2) = -d
      m(2,1) = -d; m(2,2) =  d
   end function elem_i1

   pure function elem_i2(rk, rk1) result(m)   ! ∫ ψ_a ψ'_b r² dr   (C2)
      real(wp), intent(in) :: rk, rk1
      real(wp) :: m(2,2)
      m(1,1) = -(rk1*rk1 + 2.0_wp*rk1*rk + 3.0_wp*rk*rk)/12.0_wp
      m(1,2) = -m(1,1)
      m(2,2) =  (3.0_wp*rk1*rk1 + 2.0_wp*rk1*rk + rk*rk)/12.0_wp
      m(2,1) = -m(2,2)
   end function elem_i2

   pure function elem_i3(rk, rk1) result(m)   ! ∫ ψ_a ψ'_b r dr    (C3)
      real(wp), intent(in) :: rk, rk1
      real(wp) :: m(2,2)
      m(1,1) = -(rk1 + 2.0_wp*rk)/6.0_wp
      m(1,2) = -m(1,1)
      m(2,2) =  (2.0_wp*rk1 + rk)/6.0_wp
      m(2,1) = -m(2,2)
   end function elem_i3

   pure function elem_i4(rk, rk1) result(m)   ! ∫ ψ_a ψ_b r² dr    (C4)
      real(wp), intent(in) :: rk, rk1
      real(wp) :: m(2,2), h
      h = rk1 - rk
      m(1,1) = h/30.0_wp*(rk1*rk1 + 3.0_wp*rk1*rk + 6.0_wp*rk*rk)
      m(1,2) = h/60.0_wp*(3.0_wp*rk1*rk1 + 4.0_wp*rk1*rk + 3.0_wp*rk*rk)
      m(2,1) = m(1,2)
      m(2,2) = h/30.0_wp*(6.0_wp*rk1*rk1 + 3.0_wp*rk1*rk + rk*rk)
   end function elem_i4

   pure function elem_i5(rk, rk1) result(m)   ! ∫ ψ_a ψ_b r dr     (C5)
      real(wp), intent(in) :: rk, rk1
      real(wp) :: m(2,2), h
      h = rk1 - rk
      m(1,1) = h/12.0_wp*(rk1 + 3.0_wp*rk)
      m(1,2) = h/12.0_wp*(rk1 + rk)
      m(2,1) = m(1,2)
      m(2,2) = h/12.0_wp*(3.0_wp*rk1 + rk)
   end function elem_i5

   pure function elem_i6(rk, rk1) result(m)   ! ∫ ψ_a ψ_b dr       (C6)
      real(wp), intent(in) :: rk, rk1
      real(wp) :: m(2,2), h
      h = rk1 - rk
      m(1,1) = h/3.0_wp; m(1,2) = h/6.0_wp
      m(2,1) = h/6.0_wp; m(2,2) = h/3.0_wp
   end function elem_i6

   pure function elem_i7(rk, rk1) result(m)   ! ∫ (1/r) ψ_a ψ_b dr (C7)
      real(wp), intent(in) :: rk, rk1
      real(wp) :: m(2,2), h, lg
      h  = rk1 - rk
      lg = log(rk1/rk)                          ! singular if rk = 0
      m(1,1) = (rk1*rk1/h*lg - 2.0_wp*rk1 + 0.5_wp*(rk1 + rk))/h
      m(1,2) = (-rk1*rk/h*lg + 0.5_wp*(rk1 + rk))/h
      m(2,1) = m(1,2)
      m(2,2) = (rk*rk/h*lg - 2.0_wp*rk + 0.5_wp*(rk1 + rk))/h
   end function elem_i7

   pure function elem_k1(rk, rk1) result(v)   ! ∫ ψ'_a r² dr       (C8)
      real(wp), intent(in) :: rk, rk1
      real(wp) :: v(2)
      v(1) = -(rk1*rk1 + rk1*rk + rk*rk)/3.0_wp
      v(2) = -v(1)
   end function elem_k1

   pure function elem_k2(rk, rk1) result(v)   ! ∫ ψ_a r dr         (C9)
      real(wp), intent(in) :: rk, rk1
      real(wp) :: v(2), h
      h = rk1 - rk
      v(1) = h/6.0_wp*(rk1 + 2.0_wp*rk)
      v(2) = h/6.0_wp*(2.0_wp*rk1 + rk)
   end function elem_k2

   pure function elem_k3(rk, rk1) result(v)   ! ∫ ψ_a r² dr        (C10)
      real(wp), intent(in) :: rk, rk1
      real(wp) :: v(2), h
      h = rk1 - rk
      v(1) = h/12.0_wp*(rk1*rk1 + 2.0_wp*rk1*rk + 3.0_wp*rk*rk)
      v(2) = h/12.0_wp*(3.0_wp*rk1*rk1 + 2.0_wp*rk1*rk + rk*rk)
   end function elem_k3

end module fe_radial_integrals
