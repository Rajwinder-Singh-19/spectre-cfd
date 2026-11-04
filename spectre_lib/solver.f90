module solver
   implicit none

contains
   subroutine Fractional_Step_Solver(Nx, Ny_g, Re, epsilon_svv, Cs, cfl_target, n_steps, dpdx, eps_pert, k_mode_seed, curr_rank_mpi, num_ranks_mpi, ierr_mpi)

      use, intrinsic :: iso_c_binding
      use mpi
      use fft
      use comm
      use domain
      use field
      use navier
      use poisson
      use turbulence
      use diagnostics
      use vti_writer
      implicit none

      include 'fftw3.f03'

      integer, intent(in) :: curr_rank_mpi, num_ranks_mpi
      integer, intent(inout) :: ierr_mpi

      !--------------------------------------------------
      ! Simulation parameters
      !--------------------------------------------------

      real(c_double), parameter :: PI = acos(-1.0_c_double)
      real(c_double), parameter :: Lx = 2.0_c_double*PI
      real(c_double), parameter :: Ly = 2.0_c_double

      integer, intent(in) :: Nx
      integer, intent(in) :: Ny_g

      real(c_double), intent(in) :: Re
      real(c_double), intent(in):: epsilon_svv
      real(c_double), intent(in) :: Cs
      real(c_double), intent(in)::cfl_target

      real(c_double), intent(in) :: eps_pert
      integer, intent(in) :: k_mode_seed

      integer, intent(in) :: n_steps
      real(c_double) :: dpdx
      real(c_double) :: nu
      !--------------------------------------------------
      ! Domain decomposition
      !--------------------------------------------------

      integer :: Ny
      integer :: start_row, end_row

      !--------------------------------------------------
      ! Time integration
      !--------------------------------------------------

      integer :: n

      real(c_double) :: dy
      real(c_double) :: dt
      real(c_double) :: time_s

      !--------------------------------------------------
      ! Spectral quantities
      !--------------------------------------------------

      real(c_double), allocatable :: kx(:)
      real(c_double), allocatable :: nu_sgs(:)
      real(c_double), allocatable :: nu_svv(:)
      real(c_double), allocatable :: nu_eff(:)

      real(c_double) :: kc
      real(c_double) :: k_nyq
      real(c_double) :: Q_k

      !--------------------------------------------------
      ! Velocity and pressure fields
      !--------------------------------------------------

      complex(c_double_complex), allocatable :: u(:, :), v(:, :), v_face(:, :)
      complex(c_double_complex), allocatable :: p(:, :)

      complex(c_double_complex), allocatable :: u_star(:, :)
      complex(c_double_complex), allocatable :: v_star(:, :)

      complex(c_double_complex), allocatable :: u_star_star(:, :)
      complex(c_double_complex), allocatable :: v_star_star(:, :)

      integer::i

      !--------------------------------------------------
      ! FFTW plans
      !--------------------------------------------------

      type(c_ptr) :: forward_fft_plan
      type(c_ptr) :: backward_fft_plan

      !--------------------------------------------------
      ! Domain setup
      !--------------------------------------------------

      dy = Ly/real(Ny_g, c_double)

      call Decompose_Domain(Ny_g, num_ranks_mpi, curr_rank_mpi, &
                            Ny, start_row, end_row)

      allocate (kx(Nx))

      call Wavelengths(Nx, Lx, kx)

      !--------------------------------------------------
      ! Allocate fields
      !--------------------------------------------------

      allocate (u(Ny + 2, Nx), v(Ny + 2, Nx), v_face(Ny + 2, Nx))
      allocate (p(Ny + 2, Nx))

      allocate (u_star(Ny + 2, Nx), v_star(Ny + 2, Nx))
      allocate (u_star_star(Ny + 2, Nx), v_star_star(Ny + 2, Nx))

      allocate (nu_sgs(Nx), nu_svv(Nx), nu_eff(Nx))

      !--------------------------------------------------
      ! Create FFTW plans
      !--------------------------------------------------

      call Forward_Plan(Ny, Nx, forward_fft_plan)
      call Backward_Plan(Ny, Nx, backward_fft_plan)

      !--------------------------------------------------
      ! Initialise simulation
      !--------------------------------------------------

      ! call Initialise_Flow(u, v, Ny, Ny_g, Nx, Lx, Ly, forward_fft_plan, &
      ! curr_rank_mpi, num_ranks_mpi, &
      ! start_row, 3, ierr_mpi)

      call Exchange_Halos(u, Ny, Nx, curr_rank_mpi, &
                          num_ranks_mpi, ierr_mpi)

      call Exchange_Halos(v, Ny, Nx, curr_rank_mpi, &
                          num_ranks_mpi, ierr_mpi)

      time_s = 0.0_c_double

      dt = 1e-4_c_double

      nu = 1.0_c_double/Re
      !--------------------------------------------------
      ! Calculate spectral parameters
      !--------------------------------------------------

      k_nyq = maxval(abs(kx))

      kc = (2.0_c_double/3.0_c_double)*k_nyq

      !--------------------------------------------------
      ! Main time_s integration loop
      !--------------------------------------------------

      p = cmplx(0.0_c_double, 0.0_C_DOUBLE, C_DOUBLE)
      u = cmplx(0.0_c_double, 0.0_C_DOUBLE, C_DOUBLE)
      v = cmplx(0.0_c_double, 0.0_C_DOUBLE, C_DOUBLE)

      call Seed_Perturbation(u, v, Ny, Nx, Ly, dy, &
                             kx, start_row, curr_rank_mpi, &
                             num_ranks_mpi, ierr_mpi, eps_pert, k_mode_seed)

      call MPI_Barrier(MPI_COMM_WORLD, ierr_mpi)

      n = 0
      call Write_VTI(u, v, p, n, Ny, Nx, Ny_g, dy, Lx, &
                     curr_rank_mpi, num_ranks_mpi, &
                     backward_fft_plan)
      do n = 1, n_steps

         !-----------------------------------------------
         ! 1. Calculate Chollet-Lesieur viscosity
         !-----------------------------------------------

         call Chollet_Lesieur_Viscosity(u, v, Ny, Nx, dy, kx, kc, Cs, &
                                        ierr_mpi, nu_sgs)

         !-----------------------------------------------
         ! 2. Calculate SVV viscosity
         !-----------------------------------------------

         do i = 1, Nx
            call Kolmogorov_Sharp_SVV(kx(i), epsilon_svv, kc, k_nyq, Q_k, nu_svv(i))
         end do

         !-----------------------------------------------
         ! 3. Combine viscosities
         !-----------------------------------------------

         nu_eff = nu + nu_sgs + nu_svv

         !-----------------------------------------------
         ! 4. Explicit convection
         !-----------------------------------------------

         call Explicit_Convection_Step( &
            u, v, Ny, Nx, Lx, kx, dy, dt, &
            forward_fft_plan, backward_fft_plan, &
            curr_rank_mpi, num_ranks_mpi, ierr_mpi, &
            u_star, v_star)

         !call Mass_Flux_Body_Force(u, Ny, Ny_g, Nx, dy, Ly, nu, 2.0_c_double*Lx, dpdx, curr_rank_mpi, num_ranks_mpi, ierr_mpi)

         call Body_Force_Step(u_star, Ny, Nx, dpdx, dt, curr_rank_mpi, num_ranks_mpi, ierr_mpi)

         !-----------------------------------------------
         ! 5. Implicit viscous diffusion
         !-----------------------------------------------

         call Implicit_Viscous_Step( &
            u_star, v_star, Ny, Ny_g, Nx, kx, dy, dt, &
            nu_eff, curr_rank_mpi, num_ranks_mpi, &
            ierr_mpi, u_star_star, v_star_star)

         !call Rhie_Chow_Face_Velocity(v_star_star, p, Ny, Nx, kx, nu_eff, dy, dt, curr_rank_mpi, num_ranks_mpi, ierr_mpi, v_face)

         !-----------------------------------------------
         ! 6. Pressure Poisson equation
         !-----------------------------------------------

         call Pressure_Poisson_Step(u_star_star, v_star_star, &
                                    Ny, Ny_g, Nx, kx, dy, dt, &
                                    curr_rank_mpi, num_ranks_mpi, ierr_mpi, p)

         !-----------------------------------------------
         ! 7. Velocity correction
         !-----------------------------------------------

         call Velocity_Correction_Step( &
            u_star_star, v_star_star, p, &
            Ny, Nx, kx, nu_eff, dy, dt, &
            curr_rank_mpi, num_ranks_mpi, ierr_mpi, u, v)

         !call Filter_Y(u, v, Ny, Nx, curr_rank_mpi, num_ranks_mpi, ierr_mpi)
         call Dealias_Stronger(u, Ny, Nx, Lx, kx)
         call Dealias_Stronger(v, Ny, Nx, Lx, kx)

         time_s = time_s + dt

         !-----------------------------------------------
         ! 8. Diagnostics
         !-----------------------------------------------

         !if (mod(n, 100) == 0) then

         !call Solver_Diagnostics( &
         ! u, v, nu_eff, Ny, Nx, Lx, Ly, dy, dt, kx, &
         ! backward_fft_plan, curr_rank_mpi, num_ranks_mpi, &
         ! n, time_s, ierr_mpi)

         ! end if

         if (mod(n, 500) == 0 .or. n == n_steps) then

            if (curr_rank_mpi == 0) then
               write (*, '(A,F8.4)') "Simulation Time: ", time_s
            end if

            call Write_VTI(u, v, p, n, Ny, Nx, Ny_g, dy, Lx, &
                           curr_rank_mpi, num_ranks_mpi, &
                           backward_fft_plan)
         end if

         call Adaptive_Timestep(u, v, Ny, Nx, Lx, dy, cfl_target, backward_fft_plan, ierr_mpi, dt)

      end do

      !--------------------------------------------------
      ! Cleanup
      !--------------------------------------------------

      call fftw_destroy_plan(forward_fft_plan)
      call fftw_destroy_plan(backward_fft_plan)

      deallocate (kx, nu_sgs, nu_svv, nu_eff)

      deallocate (u, v, p)
      deallocate (u_star, v_star)
      deallocate (u_star_star, v_star_star)

   end subroutine Fractional_Step_Solver

   subroutine Mass_Flux_Body_Force(u, Ny, Ny_g, Nx, dy, Ly, nu, U_target, dpdx, &
                                   curr_rank_mpi, num_ranks_mpi, ierr_mpi)
      use, intrinsic :: iso_c_binding
      use mpi
      use domain
      implicit none

      integer, intent(in) :: Ny, Ny_g, Nx
      integer, intent(in) :: curr_rank_mpi, num_ranks_mpi
      integer, intent(inout) :: ierr_mpi
      real(c_double), intent(in) :: dy, Ly, nu, U_target
      real(c_double), intent(inout) :: dpdx
      complex(c_double_complex), intent(in), dimension(Ny + 2, Nx) :: u

      real(c_double) :: local_ub, global_ub, error
      real(c_double) :: local_tau_bot, global_tau_bot
      real(c_double) :: local_tau_top, global_tau_top
      real(c_double) :: tau_avg, h
      integer :: j, start_j, end_j

      call Get_Interior_Indices(Ny, curr_rank_mpi, num_ranks_mpi, start_j, end_j)

      !--------------------------------------------------
      ! 1. Compute current bulk velocity
      !--------------------------------------------------
      local_ub = 0.0_c_double
      do j = start_j, end_j
         local_ub = local_ub + real(u(j, 1), c_double)/real(Nx, c_double)
      end do

      call MPI_Allreduce(local_ub, global_ub, 1, MPI_DOUBLE_PRECISION, &
                         MPI_SUM, MPI_COMM_WORLD, ierr_mpi)
      global_ub = global_ub/real(Ny_g, c_double)

      !--------------------------------------------------
      ! 2. Compute Wall Shear Stress
      !--------------------------------------------------
      local_tau_bot = 0.0_c_double
      local_tau_top = 0.0_c_double

      ! Rank 0 handles the bottom wall (assuming j=2 is the wall, u=0)
      if (curr_rank_mpi == 0) then
         ! du/dy at wall = u(3)/dy (since u(2)=0)
         local_tau_bot = nu*(real(u(3, 1), c_double)/real(Nx, c_double))/dy
      end if

      ! Last rank handles the top wall (assuming j=Ny+1 is the wall, u=0)
      if (curr_rank_mpi == num_ranks_mpi - 1) then
         ! du/dy at wall = -u(Ny)/dy (since u(Ny+1)=0)
         local_tau_top = -nu*(real(u(Ny, 1), c_double)/real(Nx, c_double))/dy
      end if

      call MPI_Allreduce(local_tau_bot, global_tau_bot, 1, MPI_DOUBLE_PRECISION, &
                         MPI_SUM, MPI_COMM_WORLD, ierr_mpi)
      call MPI_Allreduce(local_tau_top, global_tau_top, 1, MPI_DOUBLE_PRECISION, &
                         MPI_SUM, MPI_COMM_WORLD, ierr_mpi)

      tau_avg = 0.5_c_double*(global_tau_bot + global_tau_top)
      h = Ly/2.0_c_double  ! half channel height

      !--------------------------------------------------
      ! 3. Update dpdx using momentum balance + soft relaxation
      !--------------------------------------------------
      error = U_target - global_ub

      ! Base force balances the wall friction
      ! Note: Body_Force_Step does u = u - dt*dpdx, so dpdx is negative of the force.
      dpdx = -tau_avg/h

      ! Add a very weak proportional correction to drive mass flux error to zero
      ! A relaxation time_s of ~ 50 time_s units prevents any sloshing.
      dpdx = dpdx - 0.1_c_double*error

   end subroutine Mass_Flux_Body_Force

   subroutine Adaptive_Timestep(u, v, Ny, Nx, Lx, dy, cfl_target, &
                                backward_fft_plan, &
                                ierr_mpi, dt)
      use, intrinsic :: iso_c_binding
      use mpi
      use fft
      implicit none
      integer, intent(in) :: Ny, Nx
      integer, intent(inout) :: ierr_mpi
      real(c_double), intent(in) :: cfl_target, Lx, dy
      complex(c_double_complex), intent(in), dimension(Ny + 2, Nx) :: u, v
      type(c_ptr), intent(in) :: backward_fft_plan
      real(c_double), intent(inout) :: dt

      real(c_double) :: local_u_max, global_u_max
      real(c_double) :: local_v_max, global_v_max
      real(c_double) :: dx, dt_x, dt_y
      complex(c_double_complex), allocatable :: u_r(:, :), v_r(:, :)
      integer :: i, j

      !--------------------------------------------------
      ! 1. Transform to physical space
      !--------------------------------------------------
      allocate (u_r(Ny + 2, Nx), v_r(Ny + 2, Nx))
      call Fourier_Transform(u, Ny, Nx, backward_fft_plan, u_r)
      call Fourier_Transform(v, Ny, Nx, backward_fft_plan, v_r)

      ! Normalize (same as Solver_Diagnostics)
      u_r = u_r/real(Nx, c_double)
      v_r = v_r/real(Nx, c_double)

      !--------------------------------------------------
      ! 2. Find local maximum physical velocities
      !--------------------------------------------------
      local_u_max = 0.0_c_double
      local_v_max = 0.0_c_double

      do j = 1, Nx
         do i = 2, Ny + 1
            local_u_max = max(local_u_max, abs(real(u_r(i, j), c_double)))
            local_v_max = max(local_v_max, abs(real(v_r(i, j), c_double)))
         end do
      end do

      deallocate (u_r, v_r)

      !--------------------------------------------------
      ! 3. MPI Reduction (Note: send, recv order is fixed)
      !--------------------------------------------------
      call MPI_Allreduce(local_u_max, global_u_max, 1, MPI_DOUBLE_PRECISION, &
                         MPI_MAX, MPI_COMM_WORLD, ierr_mpi)
      call MPI_Allreduce(local_v_max, global_v_max, 1, MPI_DOUBLE_PRECISION, &
                         MPI_MAX, MPI_COMM_WORLD, ierr_mpi)

      !--------------------------------------------------
      ! 4. Compute CFL-based dt
      !--------------------------------------------------
      dx = Lx/real(Nx, c_double)

      ! Prevent division by zero if flow is perfectly at rest
      global_u_max = max(global_u_max, 1.0e-12_c_double)
      global_v_max = max(global_v_max, 1.0e-12_c_double)

      ! Advective CFL limits in x and y
      dt_x = cfl_target*dx/global_u_max
      dt_y = cfl_target*dy/global_v_max

      ! Take the minimum of the x and y constraints
      dt = min(dt_x, dt_y)

      ! Optional: impose a maximum dt so it doesn't grow unbounded
      ! in case the flow momentarily stops
      dt = min(dt, 1.0e-2_c_double)

   end subroutine Adaptive_Timestep

   subroutine Filter_Y(u, v, Ny, Nx, curr_rank_mpi, num_ranks_mpi, ierr_mpi)
      use, intrinsic :: iso_c_binding
      use comm
      use domain
      implicit none
      integer, intent(in) :: Ny, Nx, curr_rank_mpi, num_ranks_mpi
      integer, intent(inout) :: ierr_mpi
      complex(c_double_complex), intent(inout) :: u(Ny + 2, Nx)
      complex(c_double_complex), intent(inout) :: v(Ny + 2, Nx)

      real(c_double), parameter :: alpha = 0.02_c_double
      integer :: j, i, start_j, end_j

      call Exchange_Halos(u, Ny, Nx, curr_rank_mpi, num_ranks_mpi, ierr_mpi)
      call Exchange_Halos(v, Ny, Nx, curr_rank_mpi, num_ranks_mpi, ierr_mpi)

      call Get_Interior_Indices(Ny, curr_rank_mpi, num_ranks_mpi, start_j, end_j)

      do i = 1, Nx
         do j = start_j, end_j
            ! 3-point smoothing filter (damps highest y-frequencies)
            u(j, i) = u(j, i) + alpha*(u(j + 1, i) - 2.0_c_double*u(j, i) + u(j - 1, i))
            v(j, i) = v(j, i) + alpha*(v(j + 1, i) - 2.0_c_double*v(j, i) + v(j - 1, i))
         end do
      end do

      ! Re-enforce physical wall boundaries after filtering
      if (curr_rank_mpi == 0) then
         u(2, :) = (0.0_c_double, 0.0_c_double)
         v(2, :) = (0.0_c_double, 0.0_c_double)
      end if
      if (curr_rank_mpi == num_ranks_mpi - 1) then
         u(Ny + 1, :) = (0.0_c_double, 0.0_c_double)
         v(Ny + 1, :) = (0.0_c_double, 0.0_c_double)
      end if
   end subroutine Filter_Y
end module solver

