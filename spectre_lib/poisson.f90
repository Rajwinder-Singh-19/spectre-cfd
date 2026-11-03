module poisson
   implicit none

contains
   subroutine Thomas_TDMA_Solver(Ny_g, a_c, b_c, c_c, d_c, x)
      use, intrinsic :: iso_c_binding
      implicit none
      integer, intent(in) :: Ny_g
      complex(C_DOUBLE_COMPLEX), intent(in), dimension(2:Ny_g)   :: a_c
      complex(C_DOUBLE_COMPLEX), intent(in), dimension(1:Ny_g)   :: b_c
      complex(C_DOUBLE_COMPLEX), intent(in), dimension(1:Ny_g - 1) :: c_c
      complex(C_DOUBLE_COMPLEX), intent(in), dimension(1:Ny_g)   :: d_c
      complex(C_DOUBLE_COMPLEX), intent(out), dimension(1:Ny_g)   :: x

      complex(C_DOUBLE_COMPLEX), dimension(2:Ny_g)   :: a
      complex(C_DOUBLE_COMPLEX), dimension(1:Ny_g)   :: b
      complex(C_DOUBLE_COMPLEX), dimension(1:Ny_g - 1) :: c
      complex(C_DOUBLE_COMPLEX), dimension(1:Ny_g)   :: d

      integer :: i

      a = a_c
      b = b_c
      c = c_c
      d = d_c

      ! Forward elimination
      c(1) = c(1)/b(1)
      d(1) = d(1)/b(1)

      do i = 2, Ny_g - 1
         c(i) = c(i)/(b(i) - a(i)*c(i - 1))
      end do

      do i = 2, Ny_g
         d(i) = (d(i) - a(i)*d(i - 1))/(b(i) - a(i)*c(i - 1))
      end do

      ! Back substitution
      x(Ny_g) = d(Ny_g)
      do i = Ny_g - 1, 1, -1
         x(i) = d(i) - c(i)*x(i + 1)
      end do
   end subroutine Thomas_TDMA_Solver

   subroutine Viscous_Poisson_Assembler(Ny_g, kxi, dy, dt, nu, a, b, c)
      use, intrinsic :: iso_c_binding
      use turbulence
      implicit none
      integer, intent(in)::Ny_g
      real(c_double), intent(in)::nu
      real(c_double), intent(in)::dy
      real(c_double), intent(in)::dt
      real(c_double), intent(in)::kxi
      complex(C_DOUBLE_COMPLEX), intent(out), dimension(2:Ny_g)   :: a
      complex(C_DOUBLE_COMPLEX), intent(out), dimension(1:Ny_g)   :: b
      complex(C_DOUBLE_COMPLEX), intent(out), dimension(1:Ny_g - 1) :: c

      real(c_double)::inv_dy2
      inv_dy2 = 1.0_c_double/(dy*dy)

      a = cmplx(-nu*dt*inv_dy2, 0.0_C_DOUBLE, C_DOUBLE)
      b = cmplx(1 + 2*nu*dt*inv_dy2 + nu*dt*kxi*kxi, 0.0_C_DOUBLE, C_DOUBLE)
      c = cmplx(-nu*dt*inv_dy2, 0.0_C_DOUBLE, C_DOUBLE)
   end subroutine Viscous_Poisson_Assembler

   subroutine Pressure_Poisson_Assembler(Ny_g, kxi, dy, a, b, c)
      use, intrinsic :: iso_c_binding
      implicit none
      integer, intent(in)::Ny_g
      real(c_double), intent(in)::dy
      real(c_double), intent(in)::kxi
      complex(C_DOUBLE_COMPLEX), intent(out), dimension(2:Ny_g)   :: a
      complex(C_DOUBLE_COMPLEX), intent(out), dimension(1:Ny_g)   :: b
      complex(C_DOUBLE_COMPLEX), intent(out), dimension(1:Ny_g - 1) :: c

      real(c_double)::inv_dy2
      inv_dy2 = 1.0_C_DOUBLE/(dy*dy)
      b = -cmplx(kxi*kxi + 2.0_C_DOUBLE*inv_dy2, 0.0_C_DOUBLE, C_DOUBLE)
      c = cmplx(inv_dy2, 0.0_C_DOUBLE, C_DOUBLE)
      a = cmplx(inv_dy2, 0.0_C_DOUBLE, C_DOUBLE)

      c(1) = cmplx(2.0_C_DOUBLE*inv_dy2, 0.0_C_DOUBLE, C_DOUBLE)
      a(Ny_g) = cmplx(2.0_C_DOUBLE*inv_dy2, 0.0_C_DOUBLE, C_DOUBLE)

   end subroutine Pressure_Poisson_Assembler
end module poisson
