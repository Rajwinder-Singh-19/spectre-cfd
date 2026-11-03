program main
   use, intrinsic::iso_c_binding
   use mpi
   use comm
   use solver
   implicit none
   include 'fftw3.f03'
   integer::curr_rank_mpi
   integer::num_ranks_mpi
   integer::ierr_mpi
   real(8) :: start_time, end_time, local_elapsed, max_elapsed
   integer :: file_unit, io_status

   integer :: Nx
   integer :: Ny_g

   real(c_double) :: Re
   real(c_double):: epsilon_svv
   real(c_double) :: Cs
   real(c_double)::cfl_target
   integer :: n_steps
   real(c_double)::dpdx
   real(c_double)::eps_pert
   integer :: k_mode_seed
   namelist /spectre_config/ Nx, Ny_g, Re, epsilon_svv, Cs, cfl_target, n_steps, dpdx, eps_pert, k_mode_seed

   call Init_MPI(curr_rank_mpi, num_ranks_mpi, ierr_mpi)

   if (curr_rank_mpi == 0) then
      open (newunit=file_unit, file="config.nml", status="old", action="read", iostat=io_status)
      if (io_status /= 0) then
         print *, "Error: Namelist file not found."
         call MPI_Abort(MPI_COMM_WORLD, io_status, ierr_mpi)
      end if

      read (file_unit, nml=spectre_config, iostat=io_status)
      if (io_status /= 0) then
         print *, "Error reading namelist parameters."
         call MPI_Abort(MPI_COMM_WORLD, io_status, ierr_mpi)
      end if
      close (file_unit)
   end if

   ! --- BROADCAST PARAMETERS TO ALL OTHER RANKS ---
   ! Integer broadcasts (1 element each)
   call MPI_Bcast(Nx, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr_mpi)
   call MPI_Bcast(Ny_g, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr_mpi)
   call MPI_Bcast(n_steps, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr_mpi)
   call MPI_Bcast(k_mode_seed, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr_mpi)

   ! Double precision broadcasts (1 element each)
   call MPI_Bcast(Re, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierr_mpi)
   call MPI_Bcast(epsilon_svv, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierr_mpi)
   call MPI_Bcast(Cs, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierr_mpi)
   call MPI_Bcast(cfl_target, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierr_mpi)
   call MPI_Bcast(dpdx, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierr_mpi)
   call MPI_Bcast(eps_pert, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierr_mpi)

   ! --- PRINT FROM RANK 0 ONLY ---
   if (curr_rank_mpi == 0) then
      print '(A)', "========================================="
      print '(A)', "     SPECTRE SIMULATION PARAMETERS       "
      print '(A)', "========================================="
      print '(A25, I6)', " Grid Size Nx:", Nx
      print '(A25, I6)', " Grid Size Ny_g:", Ny_g
      print '(A25, F8.2)', " Reynolds Number Re:", Re
      print '(A25, E12.4)', " Epsilon SVV:", epsilon_svv
      print '(A25, F8.4)', " Chollet-Lesieur Cs:", Cs
      print '(A25, F8.4)', " Target CFL:", cfl_target
      print '(A25, I8)', " Total Time Steps:", n_steps
      print '(A25, F8.4)', "Driving Gradient:", dpdx
      print '(A25, F8.4)', "Seed Amp:", eps_pert
      print '(A25, I6)', "Seed Mode:", k_mode_seed
      print '(A25, I4)', " Running on MPI Ranks:", num_ranks_mpi
      print '(A)', "========================================="
   end if

   call MPI_Barrier(MPI_COMM_WORLD, ierr_mpi)
   
   start_time = MPI_Wtime()
   
   call Fractional_Step_Solver(Nx, Ny_g, Re, epsilon_svv, Cs, cfl_target, n_steps, dpdx, eps_pert, k_mode_seed, curr_rank_mpi, num_ranks_mpi, ierr_mpi)
   
   end_time = MPI_Wtime()
   local_elapsed = end_time - start_time
   
   call MPI_Reduce(local_elapsed, max_elapsed, 1, MPI_DOUBLE_PRECISION, &
                   MPI_MAX, 0, MPI_COMM_WORLD, ierr_mpi)

   if (curr_rank_mpi == 0) then
      print '("Execution Time = ", f12.6, " seconds.")', max_elapsed
   end if
   call MPI_Finalize(ierr_mpi)

end program main

