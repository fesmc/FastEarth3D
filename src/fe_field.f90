module fe_field
   !! Analytic field generators on the Gauss-Legendre grid: spherical ice caps and
   !! exponential basins, plus the angular-distance helper. Shared by the
   !! Martinec-2018 SLE benchmark driver and (later) fe_coupling, which both need
   !! to place idealized loads/topography on the model grid.
   !!
   !! All fields are returned in the model's spatial layout (nphi, nlat), indexed
   !! by the grid's cached longitudes sht%lon and colatitudes sht%colat [rad].
   use fe_precision, only: wp
   use fe_sht,       only: sht_grid
   implicit none
   private

   public :: angular_distance, spherical_cap, exp_basin

contains

   elemental real(wp) function angular_distance(colat, lon, colat_c, lon_c) result(d)
      !! Great-circle angular distance [rad] between a grid point (colat,lon) and a
      !! centre (colat_c,lon_c): cos d = cosθ cosθc + sinθ sinθc cos(φ-φc).
      real(wp), intent(in) :: colat, lon, colat_c, lon_c
      real(wp) :: cosd
      cosd = cos(colat)*cos(colat_c) + sin(colat)*sin(colat_c)*cos(lon - lon_c)
      d = acos(max(-1.0_wp, min(1.0_wp, cosd)))
   end function angular_distance

   subroutine spherical_cap(sht, colat_c, lon_c, alpha, h0, field)
      !! Spherical-cap height field (Martinec/giapy convention):
      !!   h(δ) = h0 · sqrt[(cos δ − cos α)/(1 − cos α)]   for δ ≤ α, else 0,
      !! where δ is the angular distance from the cap centre. Returns a height [m]
      !! (multiply by ρ_ice for a surface mass load).
      type(sht_grid), intent(in)  :: sht
      real(wp),       intent(in)  :: colat_c, lon_c   !! cap centre [rad]
      real(wp),       intent(in)  :: alpha            !! angular radius [rad]
      real(wp),       intent(in)  :: h0               !! central height [m]
      real(wp),       intent(out) :: field(:,:)       !! (nphi, nlat) [m]
      real(wp) :: ca, d
      integer  :: i, j
      ca = cos(alpha)
      do j = 1, sht%nlat
         do i = 1, sht%nphi
            d = angular_distance(sht%colat(j), sht%lon(i), colat_c, lon_c)
            if (d <= alpha) then
               field(i,j) = h0*sqrt(max((cos(d) - ca)/(1.0_wp - ca), 0.0_wp))
            else
               field(i,j) = 0.0_wp
            end if
         end do
      end do
   end subroutine spherical_cap

   subroutine exp_basin(sht, colat_b, lon_b, bmax, b0, sigma, field)
      !! Circular exponential basin topography (Martinec/giapy convention):
      !!   B(δ) = bmax − b0 · exp(−δ²/(2 σ²)),
      !! δ the angular distance from the basin centre. Returns solid-surface
      !! elevation [m] (ocean where B < 0).
      type(sht_grid), intent(in)  :: sht
      real(wp),       intent(in)  :: colat_b, lon_b   !! basin centre [rad]
      real(wp),       intent(in)  :: bmax, b0         !! plateau / depth params [m]
      real(wp),       intent(in)  :: sigma            !! angular decay rate [rad]
      real(wp),       intent(out) :: field(:,:)       !! (nphi, nlat) [m]
      real(wp) :: d
      integer  :: i, j
      do j = 1, sht%nlat
         do i = 1, sht%nphi
            d = angular_distance(sht%colat(j), sht%lon(i), colat_b, lon_b)
            field(i,j) = bmax - b0*exp(-d*d/(2.0_wp*sigma*sigma))
         end do
      end do
   end subroutine exp_basin

end module fe_field
