module diagnostics
   implicit none

contains
   subroutine Solver_Diagnostics(u, v, nu_eff, Ny, Nx, Lx, Ly, &
                                 dy, dt, kx, backward_fft_plan, &
                                 curr_rank_mpi, num_ranks_mpi, &
                                 step, time, ierr_mpi)

      use, intrinsic :: iso_c_binding
      use mpi
      use fft
      use comm
      use field
      implicit none

      integer, intent(in) :: Ny, Nx
      integer, intent(in) :: curr_rank_mpi, num_ranks_mpi
      integer, intent(in) :: step
      integer, intent(inout) :: ierr_mpi

      real(c_double), intent(in) :: Lx, Ly, dy, dt, time
      real(c_double), intent(in) :: kx(Nx)
      real(c_double), intent(in) :: nu_eff(Nx)

      complex(c_double_complex), intent(in) :: u(Ny + 2, Nx)
      complex(c_double_complex), intent(in) :: v(Ny + 2, Nx)

      type(c_ptr), intent(in) :: backward_fft_plan

      complex(c_double_complex), allocatable :: u_r(:, :), v_r(:, :)
      complex(c_double_complex), allocatable :: dudx(:, :), dvdy(:, :)
      complex(c_double_complex), allocatable :: dvdx(:, :), dudy(:, :)

      complex(c_double_complex), allocatable :: dudx_r(:, :)
      complex(c_double_complex), allocatable :: dvdy_r(:, :)
      complex(c_double_complex), allocatable :: dvdx_r(:, :)
      complex(c_double_complex), allocatable :: dudy_r(:, :)

      complex(c_double_complex), allocatable :: div_r(:, :)

      real(c_double) :: dx

      real(c_double) :: local_energy, local_enstrophy
      real(c_double) :: local_umax, local_vmax
      real(c_double) :: local_divmax, local_cfl

      real(c_double) :: global_energy, global_enstrophy
      real(c_double) :: global_umax, global_vmax
      real(c_double) :: global_divmax, global_cfl

      real(c_double) :: local_numax, local_numin
      real(c_double) :: global_numax, global_numin

      real(c_double) :: uval, vval
      real(c_double) :: omega_val, div_val

      integer :: i, j

      dx = Lx/real(Nx, c_double)

      allocate (u_r(Ny + 2, Nx), v_r(Ny + 2, Nx))

      allocate (dudx(Ny + 2, Nx), dvdy(Ny + 2, Nx))
      allocate (dvdx(Ny + 2, Nx), dudy(Ny + 2, Nx))

      allocate (dudx_r(Ny + 2, Nx), dvdy_r(Ny + 2, Nx))
      allocate (dvdx_r(Ny + 2, Nx), dudy_r(Ny + 2, Nx))

      allocate (div_r(Ny + 2, Nx))

      !-----------------------------------------------
      ! Transform velocity into physical space
      !-----------------------------------------------

      call Fourier_Transform(u, Ny, Nx, backward_fft_plan, u_r)
      call Fourier_Transform(v, Ny, Nx, backward_fft_plan, v_r)

      u_r = u_r/real(Nx, c_double)
      v_r = v_r/real(Nx, c_double)

      !-----------------------------------------------
      ! Calculate spectral derivatives
      !-----------------------------------------------

      call Dx_Spectral(u, Ny, Nx, kx, curr_rank_mpi, &
                       num_ranks_mpi, dudx)

      call Dx_Spectral(v, Ny, Nx, kx, curr_rank_mpi, &
                       num_ranks_mpi, dvdx)

      call Dy_Center(u, Ny, Nx, dy, curr_rank_mpi, &
                     num_ranks_mpi, dudy)

      call Dy_Center(v, Ny, Nx, dy, curr_rank_mpi, &
                     num_ranks_mpi, dvdy)

      !-----------------------------------------------
      ! Transform derivatives into physical space
      !-----------------------------------------------

      call Fourier_Transform(dudx, Ny, Nx, backward_fft_plan, dudx_r)
      call Fourier_Transform(dvdx, Ny, Nx, backward_fft_plan, dvdx_r)

      call Fourier_Transform(dudy, Ny, Nx, backward_fft_plan, dudy_r)
      call Fourier_Transform(dvdy, Ny, Nx, backward_fft_plan, dvdy_r)

      dudx_r = dudx_r/real(Nx, c_double)
      dvdx_r = dvdx_r/real(Nx, c_double)

      dudy_r = dudy_r/real(Nx, c_double)
      dvdy_r = dvdy_r/real(Nx, c_double)

      !-----------------------------------------------
      ! Physical-space divergence
      !-----------------------------------------------

      div_r = dudx_r + dvdy_r

      !-----------------------------------------------
      ! Initialise diagnostics
      !-----------------------------------------------

      local_energy = 0.0_c_double
      local_enstrophy = 0.0_c_double

      local_umax = 0.0_c_double
      local_vmax = 0.0_c_double

      local_divmax = 0.0_c_double
      local_cfl = 0.0_c_double

      local_numin = minval(nu_eff)
      local_numax = maxval(nu_eff)

      !-----------------------------------------------
      ! Calculate local diagnostics
      !-----------------------------------------------

      do j = 2, Ny + 1
         do i = 1, Nx

            uval = real(u_r(j, i), c_double)
            vval = real(v_r(j, i), c_double)

            omega_val = real(dvdx_r(j, i) - dudy_r(j, i), c_double)

            div_val = abs(real(div_r(j, i), c_double))

            ! Kinetic energy
            local_energy = local_energy + &
                           0.5_c_double*(uval*uval + vval*vval)*dx*dy

            ! Enstrophy
            local_enstrophy = local_enstrophy + &
                              0.5_c_double*omega_val*omega_val*dx*dy

            ! Maximum velocities
            local_umax = max(local_umax, abs(uval))
            local_vmax = max(local_vmax, abs(vval))

            ! Maximum divergence
            local_divmax = max(local_divmax, div_val)

            ! CFL number
            local_cfl = max(local_cfl, &
                            abs(uval)*dt/dx + abs(vval)*dt/dy)

         end do
      end do

      !-----------------------------------------------
      ! Global MPI reductions
      !-----------------------------------------------

      call MPI_Allreduce(local_energy, global_energy, 1, &
                         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr_mpi)

      call MPI_Allreduce(local_enstrophy, global_enstrophy, 1, &
                         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr_mpi)

      call MPI_Allreduce(local_umax, global_umax, 1, &
                         MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr_mpi)

      call MPI_Allreduce(local_vmax, global_vmax, 1, &
                         MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr_mpi)

      call MPI_Allreduce(local_divmax, global_divmax, 1, &
                         MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr_mpi)

      call MPI_Allreduce(local_cfl, global_cfl, 1, &
                         MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr_mpi)

      call MPI_Allreduce(local_numin, global_numin, 1, &
                         MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr_mpi)

      call MPI_Allreduce(local_numax, global_numax, 1, &
                         MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr_mpi)

      !-----------------------------------------------
      ! Print diagnostics
      !-----------------------------------------------

      if (curr_rank_mpi == 0) then

         print *
         print '(A)', '=============================================='
         print '(A,I8)', ' STEP: ', step
         print '(A,ES14.6)', ' TIME: ', time
         print '(A)', '----------------------------------------------'

         print '(A,ES14.6)', ' Kinetic energy       : ', &
            global_energy/(Lx*Ly)

         print '(A,ES14.6)', ' Enstrophy            : ', &
            global_enstrophy/(Lx*Ly)

         print '(A,ES14.6)', ' Maximum |u|          : ', global_umax
         print '(A,ES14.6)', ' Maximum |v|          : ', global_vmax

         print '(A,ES14.6)', ' Maximum divergence   : ', global_divmax

         print '(A,ES14.6)', ' CFL number           : ', global_cfl

         print '(A,ES14.6)', ' Minimum nu_eff       : ', global_numin
         print '(A,ES14.6)', ' Maximum nu_eff       : ', global_numax

         print '(A)', '=============================================='

      end if

      !-----------------------------------------------
      ! Cleanup
      !-----------------------------------------------

      deallocate (u_r, v_r)

      deallocate (dudx, dvdy, dvdx, dudy)

      deallocate (dudx_r, dvdy_r, dvdx_r, dudy_r)

      deallocate (div_r)

   end subroutine Solver_Diagnostics
end module diagnostics
