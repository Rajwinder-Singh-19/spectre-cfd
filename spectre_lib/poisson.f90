module poisson
   use, intrinsic :: iso_c_binding
   implicit none

   ! Pressure Poisson matrices, assembled ONCE before the time loop.
   ! The pressure matrix depends only on kx and dy (not on dt or the viscosity), so it
   ! is identical at every time step. Column i holds the tridiagonal matrix of Fourier
   ! mode i.
   ! Every rank solves its own share of the modes (see comm.f90,
   ! Transpose_To_Modes), so every rank fills and stores only the columns of its own modes.
   complex(C_DOUBLE_COMPLEX), allocatable, save :: p_mat_a(:, :)   ! sub-diagonal   (2:Ny_g, Nx)
   complex(C_DOUBLE_COMPLEX), allocatable, save :: p_mat_b(:, :)   ! diagonal       (1:Ny_g, Nx)
   complex(C_DOUBLE_COMPLEX), allocatable, save :: p_mat_c(:, :)   ! super-diagonal (1:Ny_g-1, Nx)

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

   ! [TEST3 CHANGE] ---------------------------------------------------------------
   ! Assemble the pressure Poisson matrix of every Fourier mode once and store it.
   ! Called once from Fractional_Step_Solver before the time loop. The coefficients
   ! stay complex, exactly as in the original code. The kx = 0 gauge condition on the
   ! matrix (p(1) = 0) is applied here; the matching right-hand-side entry stays in
   ! Pressure_Poisson_Step because the right-hand side changes every step.
   subroutine Pressure_Matrix_Precompute(Ny_g, Nx, kx, dy, curr_rank_mpi, num_ranks_mpi)
      use domain
      implicit none
      integer, intent(in) :: Ny_g, Nx, curr_rank_mpi, num_ranks_mpi
      real(c_double), intent(in) :: kx(Nx)
      real(c_double), intent(in) :: dy
      integer :: i, nm, m_start, m_end

      ! modes owned by this rank (same split as in Transpose_To_Modes)
      call Decompose_Domain(Nx, num_ranks_mpi, curr_rank_mpi, nm, m_start, m_end)

      allocate (p_mat_a(2:Ny_g, Nx))
      allocate (p_mat_b(1:Ny_g, Nx))
      allocate (p_mat_c(1:Ny_g - 1, Nx))

      do i = m_start, m_end
         call Pressure_Poisson_Assembler(Ny_g, kx(i), dy, &
                                         p_mat_a(:, i), p_mat_b(:, i), p_mat_c(:, i))

         if (abs(kx(i)) < 1.0e-12_c_double) then
            ! Fix the pressure gauge: p(1) = 0.
            p_mat_b(1, i) = cmplx(1.0_c_double, 0.0_c_double, kind=c_double)
            p_mat_c(1, i) = cmplx(0.0_c_double, 0.0_c_double, kind=c_double)
         end if
      end do
   end subroutine Pressure_Matrix_Precompute
end module poisson
