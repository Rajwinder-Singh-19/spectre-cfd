module navier
   implicit none

contains
   subroutine Explicit_Convection_Step(u, v, Ny, Nx, Lx, kx, dy, dt, &
                                       forward_fft_plan, backward_fft_plan, &
                                       curr_rank_mpi, num_ranks_mpi, ierr_mpi, u_star, v_star)
      use, intrinsic :: iso_c_binding
      use domain
      use field
      use fft
      use comm
      implicit none

      integer, intent(in) :: Ny, Nx
      real(c_double), intent(in) :: kx(Nx)
      real(c_double), intent(in) :: dy, dt, Lx
      integer, intent(in) :: curr_rank_mpi, num_ranks_mpi
      integer, intent(inout) :: ierr_mpi
      complex(c_double_complex), intent(in), dimension(Ny + 2, Nx) :: u, v
      complex(c_double_complex), intent(out), dimension(Ny + 2, Nx) :: u_star, v_star

      complex(c_double_complex), allocatable, dimension(:, :) :: du_dx, du_dy
      complex(c_double_complex), allocatable, dimension(:, :) :: dv_dx, dv_dy

      complex(c_double_complex), allocatable, dimension(:, :) :: du_dx_r, du_dy_r
      complex(c_double_complex), allocatable, dimension(:, :) :: dv_dx_r, dv_dy_r

      complex(c_double_complex), allocatable, dimension(:, :) :: u_r
      complex(c_double_complex), allocatable, dimension(:, :) :: v_r

      complex(c_double_complex), allocatable, dimension(:, :) :: conv_u_r
      complex(c_double_complex), allocatable, dimension(:, :) :: conv_v_r

      complex(c_double_complex), allocatable :: u_hat(:, :), v_hat(:, :)
      real(c_double) :: inv_Nx

      type(c_ptr), intent(in)::forward_fft_plan
      type(c_ptr), intent(in)::backward_fft_plan

      allocate (du_dx(Ny + 2, Nx), du_dy(Ny + 2, Nx))
      allocate (dv_dx(Ny + 2, Nx), dv_dy(Ny + 2, Nx))
      allocate (du_dx_r(Ny + 2, Nx), du_dy_r(Ny + 2, Nx))
      allocate (dv_dx_r(Ny + 2, Nx), dv_dy_r(Ny + 2, Nx))

      allocate (u_r(Ny + 2, Nx))
      allocate (v_r(Ny + 2, Nx))

      allocate (conv_u_r(Ny + 2, Nx))
      allocate (conv_v_r(Ny + 2, Nx))
      allocate (u_hat(Ny + 2, Nx), v_hat(Ny + 2, Nx))

      u_hat = u
      v_hat = v
      inv_Nx = 1.0_c_double/real(Nx, c_double)

      ! 1. Dealias and compute spectral spatial derivatives
      call Dealias_Stronger(u_hat, Ny, Nx, Lx, kx)
      call Dealias_Stronger(v_hat, Ny, Nx, Lx, kx)

      call Dx_Spectral(u_hat, Ny, Nx, kx, curr_rank_mpi, num_ranks_mpi, du_dx)
      call Dx_Spectral(v_hat, Ny, Nx, kx, curr_rank_mpi, num_ranks_mpi, dv_dx)

      call Exchange_Halos(u_hat, Ny, nx, curr_rank_mpi, num_ranks_mpi, ierr_mpi)
      call Dy_Center(u_hat, Ny, Nx, dy, curr_rank_mpi, num_ranks_mpi, du_dy)
      call Exchange_Halos(v_hat, Ny, nx, curr_rank_mpi, num_ranks_mpi, ierr_mpi)
      call Dy_Center(v_hat, Ny, Nx, dy, curr_rank_mpi, num_ranks_mpi, dv_dy)

      ! 2. Transform fields into physical real space
      call Fourier_Transform(u_hat, Ny, Nx, backward_fft_plan, u_r)
      call Fourier_Transform(v_hat, Ny, Nx, backward_fft_plan, v_r)
      call Fourier_Transform(du_dx, Ny, Nx, backward_fft_plan, du_dx_r)
      call Fourier_Transform(dv_dx, Ny, Nx, backward_fft_plan, dv_dx_r)
      call Fourier_Transform(du_dy, Ny, Nx, backward_fft_plan, du_dy_r)
      call Fourier_Transform(dv_dy, Ny, Nx, backward_fft_plan, dv_dy_r)

      ! 3. Standardize unnormalized backward transforms to true physical values
      u_r = u_r*inv_Nx
      v_r = v_r*inv_Nx
      du_dx_r = du_dx_r*inv_Nx
      dv_dx_r = dv_dx_r*inv_Nx
      du_dy_r = du_dy_r*inv_Nx
      dv_dy_r = dv_dy_r*inv_Nx

      ! 4. Evaluate non-linear products in real space
      conv_u_r = cmplx(real(u_r, c_double)*real(du_dx_r, c_double) + real(v_r, c_double)*real(du_dy_r, c_double), 0.0_c_double, c_double_complex)
      conv_v_r = cmplx(real(u_r, c_double)*real(dv_dx_r, c_double) + real(v_r, c_double)*real(dv_dy_r, c_double), 0.0_c_double, c_double_complex)

      call Fourier_Transform(conv_u_r, Ny, Nx, forward_fft_plan, u_star)
      call Fourier_Transform(conv_v_r, Ny, Nx, forward_fft_plan, v_star)

      call Dealias_Stronger(u_star, Ny, Nx, Lx, kx)
      call Dealias_Stronger(v_star, Ny, Nx, Lx, kx)

      u_star = u_hat - dt*u_star
      v_star = v_hat - dt*v_star

      if (curr_rank_mpi == 0) then
         u_star(2, :) = cmplx(0.0_c_double, 0.0_c_double, c_double)
         v_star(2, :) = cmplx(0.0_c_double, 0.0_c_double, c_double)
      end if

      if (curr_rank_mpi == num_ranks_mpi - 1) then
         u_star(Ny + 1, :) = cmplx(0.0_c_double, 0.0_c_double, c_double)
         v_star(Ny + 1, :) = cmplx(0.0_c_double, 0.0_c_double, c_double)
      end if

      deallocate (du_dx, du_dy, dv_dx, dv_dy)
      deallocate (du_dx_r, du_dy_r, dv_dx_r, dv_dy_r)
      deallocate (u_r, v_r, conv_u_r, conv_v_r)
      deallocate (u_hat, v_hat)

   end subroutine Explicit_Convection_Step

   subroutine Implicit_Viscous_Step(u_star, v_star, Ny, Ny_g, Nx, kx, dy, dt, nu &
                                    , curr_rank_mpi, num_ranks_mpi, ierr_mpi, &
                                    u_star_star, v_star_star)

      use, intrinsic :: iso_c_binding
      use comm
      use poisson
      implicit none

      integer, intent(in) :: Ny, Nx, Ny_g
      real(c_double), intent(in) :: kx(Nx)
      real(c_double), intent(in) :: dy, dt, nu(Nx)
      integer, intent(in) :: curr_rank_mpi, num_ranks_mpi
      integer, intent(inout) :: ierr_mpi

      complex(c_double_complex), intent(inout) :: u_star(Ny + 2, Nx)
      complex(c_double_complex), intent(inout) :: v_star(Ny + 2, Nx)

      complex(c_double_complex), intent(out) :: u_star_star(Ny + 2, Nx)
      complex(c_double_complex), intent(out) :: v_star_star(Ny + 2, Nx)

      ! Tridiagonal matrix coefficients
      complex(c_double_complex), allocatable :: a(:)
      complex(c_double_complex), allocatable :: b(:)
      complex(c_double_complex), allocatable :: c(:)

      ! Gathered right-hand sides
      complex(c_double_complex), allocatable :: u_rhs(:)
      complex(c_double_complex), allocatable :: v_rhs(:)

      ! Globally solved columns
      complex(c_double_complex), allocatable :: u_glob(:)
      complex(c_double_complex), allocatable :: v_glob(:)

      integer :: i

      !--------------------------------------------------
      ! Allocate local arrays
      !--------------------------------------------------

      allocate (a(2:Ny_g))
      allocate (b(1:Ny_g))
      allocate (c(1:Ny_g - 1))

      allocate (u_rhs(1:Ny_g))
      allocate (v_rhs(1:Ny_g))

      allocate (u_glob(1:Ny_g))
      allocate (v_glob(1:Ny_g))

      !--------------------------------------------------
      ! Exchange halos
      !--------------------------------------------------

      call Exchange_Halos(u_star, Ny, Nx, curr_rank_mpi, &
                          num_ranks_mpi, ierr_mpi)

      call Exchange_Halos(v_star, Ny, Nx, curr_rank_mpi, &
                          num_ranks_mpi, ierr_mpi)

      !--------------------------------------------------
      ! Solve each Fourier mode independently
      !--------------------------------------------------

      do i = 1, Nx

         call Gather_Column_To_0(u_star, i, Ny, Ny_g, Nx, &
                                 curr_rank_mpi, num_ranks_mpi, &
                                 u_rhs, ierr_mpi)

         call Gather_Column_To_0(v_star, i, Ny, Ny_g, Nx, &
                                 curr_rank_mpi, num_ranks_mpi, &
                                 v_rhs, ierr_mpi)

         if (curr_rank_mpi == 0) then

            ! Assemble viscous matrix
            call Viscous_Poisson_Assembler(Ny_g, kx(i), &
                                           dy, dt, nu(i), a, b, c)

            ! Dirichlet boundary conditions
            u_rhs(1) = cmplx(0.0_c_double, 0.0_c_double, kind=c_double)
            u_rhs(Ny_g) = cmplx(0.0_c_double, 0.0_c_double, kind=c_double)

            v_rhs(1) = cmplx(0.0_c_double, 0.0_c_double, kind=c_double)
            v_rhs(Ny_g) = cmplx(0.0_c_double, 0.0_c_double, kind=c_double)

            b(1) = cmplx(1.0_c_double, 0.0_c_double, kind=c_double)
            c(1) = cmplx(0.0_c_double, 0.0_c_double, kind=c_double)

            b(Ny_g) = cmplx(1.0_c_double, 0.0_c_double, kind=c_double)
            a(Ny_g) = cmplx(0.0_c_double, 0.0_c_double, kind=c_double)

            ! Solve for u and v
            call Thomas_TDMA_Solver(Ny_g, a, b, c, u_rhs, u_glob)
            call Thomas_TDMA_Solver(Ny_g, a, b, c, v_rhs, v_glob)

         end if

         ! Distribute solutions
         call Scatter_Column_From_0(u_glob, Ny_g, u_star_star, i, &
                                    Ny, Nx, curr_rank_mpi, &
                                    num_ranks_mpi, ierr_mpi)

         call Scatter_Column_From_0(v_glob, Ny_g, v_star_star, i, &
                                    Ny, Nx, curr_rank_mpi, &
                                    num_ranks_mpi, ierr_mpi)

      end do

      if (curr_rank_mpi == 0) then
         u_star_star(2, :) = cmplx(0.0_c_double, 0.0_c_double, c_double)
         v_star_star(2, :) = cmplx(0.0_c_double, 0.0_c_double, c_double)
      end if

      if (curr_rank_mpi == num_ranks_mpi - 1) then
         u_star_star(Ny + 1, :) = cmplx(0.0_c_double, 0.0_c_double, c_double)
         v_star_star(Ny + 1, :) = cmplx(0.0_c_double, 0.0_c_double, c_double)
      end if

      !--------------------------------------------------
      ! Exchange output halos
      !--------------------------------------------------

      call Exchange_Halos(u_star_star, Ny, Nx, curr_rank_mpi, &
                          num_ranks_mpi, ierr_mpi)

      call Exchange_Halos(v_star_star, Ny, Nx, curr_rank_mpi, &
                          num_ranks_mpi, ierr_mpi)

      !--------------------------------------------------
      ! Deallocate
      !--------------------------------------------------

      deallocate (a, b, c)
      deallocate (u_rhs, v_rhs)
      deallocate (u_glob, v_glob)

   end subroutine Implicit_Viscous_Step

   subroutine Body_Force_Step(u_star, Ny, Nx, dpdx, dt, curr_rank_mpi, num_ranks_mpi, ierr_mpi)
      use, intrinsic::iso_c_binding
      use domain
      use comm
      implicit none

      integer, intent(in) :: Ny, Nx, curr_rank_mpi, num_ranks_mpi
      integer, intent(inout) :: ierr_mpi
      real(C_DOUBLE), intent(in) :: dpdx, dt
      complex(C_DOUBLE_COMPLEX), dimension(Ny + 2, Nx), intent(inout) :: u_star
      integer :: i, start_i, end_i

      call Get_Interior_Indices(Ny, curr_rank_mpi, num_ranks_mpi, start_i, end_i)

      do i = start_i, end_i
         u_star(i, 1) = u_star(i, 1) - dt*real(Nx, c_double)*complex(dpdx, 0.0_C_DOUBLE)
      end do

      call Exchange_Halos(u_star, Ny, Nx, curr_rank_mpi, num_ranks_mpi, ierr_mpi)

      if (curr_rank_mpi == 0) then
         u_star(2, :) = cmplx(0.0_c_double, 0.0_c_double, c_double)
      end if

      if (curr_rank_mpi == num_ranks_mpi - 1) then
         u_star(Ny + 1, :) = cmplx(0.0_c_double, 0.0_c_double, c_double)
      end if
   end subroutine

   subroutine Pressure_Poisson_Step(u_star_star, v_star_star, Ny, Ny_g, Nx, &
                                    kx, dy, dt, curr_rank_mpi, num_ranks_mpi, &
                                    ierr_mpi, p)

      use, intrinsic :: iso_c_binding
      use comm
      use poisson
      use field
      implicit none

      integer, intent(in) :: Ny, Ny_g, Nx
      real(c_double), intent(in) :: kx(Nx)
      real(c_double), intent(in) :: dy, dt

      integer, intent(in) :: curr_rank_mpi, num_ranks_mpi
      integer, intent(inout) :: ierr_mpi

      complex(c_double_complex), intent(inout) :: u_star_star(Ny + 2, Nx)
      complex(c_double_complex), intent(inout) :: v_star_star(Ny + 2, Nx)

      complex(c_double_complex), intent(out) :: p(Ny + 2, Nx)

      ! Derivatives
      complex(c_double_complex), allocatable :: div_uv(:, :)

      ! Tridiagonal coefficients
      complex(c_double_complex), allocatable :: a(:)
      complex(c_double_complex), allocatable :: b(:)
      complex(c_double_complex), allocatable :: c(:)

      ! Gathered RHS and pressure solution
      complex(c_double_complex), allocatable :: rhs(:)
      complex(c_double_complex), allocatable :: p_glob(:)
      integer :: i

      !--------------------------------------------------
      ! Allocate arrays
      !--------------------------------------------------

      allocate (div_uv(Ny + 2, Nx))

      allocate (a(2:Ny_g))
      allocate (b(1:Ny_g))
      allocate (c(1:Ny_g - 1))

      allocate (rhs(Ny_g))
      allocate (p_glob(Ny_g))

      !--------------------------------------------------
      ! Calculate velocity divergence
      !--------------------------------------------------

      call Exchange_Halos(u_star_star, Ny, Nx, curr_rank_mpi, &
                          num_ranks_mpi, ierr_mpi)

      call Exchange_Halos(v_star_star, Ny, Nx, curr_rank_mpi, &
                          num_ranks_mpi, ierr_mpi)
      call Divergence(u_star_star, v_star_star, Ny, Nx, kx, dy, curr_rank_mpi, num_ranks_mpi, div_uv)

      call Exchange_Halos(div_uv, Ny, Nx, curr_rank_mpi, num_ranks_mpi, ierr_mpi)

      !--------------------------------------------------
      ! Solve each Fourier mode
      !--------------------------------------------------

      do i = 1, Nx

         ! Gather divergence of velocity
         ! RHS = (du/dx + dv/dy)/dt

         call Gather_Column_To_0(div_uv, i, Ny, Ny_g, Nx, &
                                 curr_rank_mpi, num_ranks_mpi, &
                                 rhs, ierr_mpi)

         if (curr_rank_mpi == 0) then

            rhs = rhs/dt

            ! Assemble pressure Poisson matrix
            call Pressure_Poisson_Assembler(Ny_g, kx(i), dy, a, b, c)

            !--------------------------------------------------
            ! Neumann boundary conditions
            !--------------------------------------------------

            if (abs(kx(i)) < 1.0e-12_c_double) then

               ! Fix the pressure gauge: p(1) = 0.
               b(1) = cmplx(1.0_c_double, 0.0_c_double, &
                            kind=c_double)

               c(1) = cmplx(0.0_c_double, 0.0_c_double, &
                            kind=c_double)

               rhs(1) = cmplx(0.0_c_double, 0.0_c_double, &
                              kind=c_double)

            end if

            ! Solve pressure Poisson equation

            call Thomas_TDMA_Solver(Ny_g, a, b, c, rhs, &
                                    p_glob)

         end if

         ! Scatter pressure solution
         call Scatter_Column_From_0(p_glob, Ny_g, p, i, Ny, Nx, &
                                    curr_rank_mpi, num_ranks_mpi, ierr_mpi)

      end do

      !--------------------------------------------------
      ! Exchange pressure halos
      !--------------------------------------------------

      call Exchange_Halos(p, Ny, Nx, curr_rank_mpi, &
                          num_ranks_mpi, ierr_mpi)

      !--------------------------------------------------
      ! Cleanup
      !--------------------------------------------------

      deallocate (div_uv)
      deallocate (a, b, c)
      deallocate (rhs, p_glob)

   end subroutine Pressure_Poisson_Step

   subroutine Velocity_Correction_Step(u_star_star, v_star_star, p, &
                                       Ny, Nx, kx, nu, dy, dt, &
                                       curr_rank_mpi, num_ranks_mpi, ierr_mpi, &
                                       u_new, v_new)

      use, intrinsic :: iso_c_binding
      use field
      use domain
      use comm
      implicit none

      integer, intent(in) :: Ny, Nx
      integer, intent(in) :: curr_rank_mpi, num_ranks_mpi
      integer, intent(inout) :: ierr_mpi

      real(c_double), intent(in) :: kx(Nx)
      real(c_double), intent(in)::nu(Nx)
      real(c_double), intent(in) :: dy, dt

      complex(c_double_complex), intent(in) :: &
         u_star_star(Ny + 2, Nx), v_star_star(Ny + 2, Nx)

      complex(c_double_complex), intent(inout) :: p(Ny + 2, Nx)

      complex(c_double_complex), intent(out) :: &
         u_new(Ny + 2, Nx), v_new(Ny + 2, Nx)

      complex(c_double_complex), allocatable :: dp_dx(:, :)
      complex(c_double_complex), allocatable :: dp_dy(:, :)

      integer::i
      real(c_double), dimension(Nx)::a
      real(c_double)::inv_d2

      inv_d2 = 1.0_c_double/(dy*dy)
      !--------------------------------------------------
      ! Allocate arrays
      !--------------------------------------------------

      allocate (dp_dx(Ny + 2, Nx))
      allocate (dp_dy(Ny + 2, Nx))

      dp_dx = cmplx(0.0_c_double, 0.0_c_double, kind=c_double)
      dp_dy = cmplx(0.0_c_double, 0.0_c_double, kind=c_double)

      !--------------------------------------------------
      ! Exchange pressure halos
      !--------------------------------------------------

      call Exchange_Halos(p, Ny, Nx, &
                          curr_rank_mpi, num_ranks_mpi, ierr_mpi)

      !--------------------------------------------------
      ! Calculate pressure derivatives
      !--------------------------------------------------

      ! Spectral derivative in x

      call Dx_Spectral(p, Ny, Nx, kx, &
                       curr_rank_mpi, num_ranks_mpi, dp_dx)

      ! Finite difference derivative in y

      call Dy_Center(p, Ny, Nx, dy, &
                     curr_rank_mpi, num_ranks_mpi, dp_dy)

      !--------------------------------------------------
      ! Velocity correction
      !--------------------------------------------------

      do i = 1, Nx
         a(i) = 1.0_c_double + nu(i)*dt*(2.0_c_double*inv_d2 + kx(i)*kx(i))
      end do

      do i = 1, Nx
         u_new(:, i) = u_star_star(:, i) - (dt)*dp_dx(:, i)
         v_new(:, i) = v_star_star(:, i) - (dt)*dp_dy(:, i)
      end do

      call Exchange_Halos(u_new, Ny, Nx, curr_rank_mpi, num_ranks_mpi, ierr_mpi)
      call Exchange_Halos(v_new, Ny, Nx, curr_rank_mpi, num_ranks_mpi, ierr_mpi)
      !--------------------------------------------------
      ! Cleanup
      !--------------------------------------------------

      deallocate (dp_dx, dp_dy)

   end subroutine Velocity_Correction_Step

end module navier
