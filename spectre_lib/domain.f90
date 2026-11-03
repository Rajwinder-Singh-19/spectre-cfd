module domain
   implicit none

contains
   subroutine Decompose_Domain(Ny_g, num_ranks_mpi, curr_rank_mpi, Ny, l2g_start, l2g_end)
      ! This subroutine calculates the number of rows owned by each processor rank based on the number of processors
      ! Each processor gets at least (Ny) / (num. of ranks) rows, all the remaining rows are allocated to all the processors one by one
      ! Local to global array map is also calculated which indicates how each each local array is placed in the global array
      implicit none
      integer, intent(in)::Ny_g
      integer, intent(in)::num_ranks_mpi
      integer, intent(in)::curr_rank_mpi
      integer, intent(out)::Ny
      integer, intent(out)::l2g_start
      integer, intent(out)::l2g_end

      integer::Ny_base, Ny_rem

      Ny_base = Ny_g/num_ranks_mpi ! Each processor gets at least this many rows
      Ny_rem = mod(Ny_g, num_ranks_mpi) ! Remaining number of extra rows which need to be divided amongst all processors

      ! 0 <= Ny_rem <= num_ranks_mpi (property of mod operation)

      if (curr_rank_mpi .LT. Ny_rem) then
         Ny = Ny_base + 1 ! All the ranks get one extra row until the rank is greater than all the extra rows (property of mod)
      else
         Ny = Ny_base
      end if

      l2g_start = curr_rank_mpi*Ny_base + min(curr_rank_mpi, Ny_rem) + 1 ! Because each rank gets 1 extra row, offset cumulatively adds up (+1 for fortran indexing)
      l2g_end = l2g_start + Ny - 1 ! Last l2g map is at a distance of Ny (-1 for fortran indexing)

   end subroutine Decompose_Domain

   subroutine Get_Neighbour_Ranks(curr_rank_mpi, num_ranks_mpi, north_rank_mpi, south_rank_mpi)
      ! This subroutine gets the ranks of the processor handling the subdomain above (north) and below (south) the domain handled by the current rank
      ! For all the interior ranks, north rank is (current rank + 1), and south rank is (current rank - 1)
      ! Boundary ranks need special treatment, which is
      ! 1. Lower boundary rank has no south rank, i.e. south rank is NULL
      ! 2. Upper boundary rank has no north rank, i.e. north rank is NULL
      use mpi
      implicit none
      integer, intent(in)::curr_rank_mpi! Rank of the processor calling the subroutine
      integer, intent(in)::num_ranks_mpi! Total number of ranks
      integer, intent(out)::north_rank_mpi ! Rank of the processor handling subdomain above the current subdomain
      integer, intent(out)::south_rank_mpi ! Rank of the processor handling subdomain below the current subdomain

      north_rank_mpi = curr_rank_mpi + 1
      south_rank_mpi = curr_rank_mpi - 1

      if (curr_rank_mpi .EQ. 0) then
         south_rank_mpi = MPI_PROC_NULL ! Lower boundary rank has no south rank
      end if

      if (curr_rank_mpi .EQ. num_ranks_mpi - 1) then
         north_rank_mpi = MPI_PROC_NULL ! Upper boundary rank has no north rank
      end if

   end subroutine Get_Neighbour_Ranks

   subroutine Get_Interior_Indices(Ny, curr_rank_mpi, num_ranks_mpi, start_interior_row, end_interior_row)
      ! This subroutine finds the row index of the first and last interior row, i.e. row which is neither a ghost row or a boundary row
      ! For all the ranks, interior rows are between 2 and Ny+1 inclusive, only the boundary processors need to be handled, which are
      ! 1. Lower boundary processor interior rows start at 3
      ! 2. Upper boundary procerssor interior rows end at Ny
      implicit none
      integer, intent(in)::Ny ! Number of local rows
      integer, intent(in)::curr_rank_mpi ! Rank of the processor calling the subroutine
      integer, intent(in)::num_ranks_mpi ! Total number of ranks
      integer, intent(out)::start_interior_row ! Interior row starting index
      integer, intent(out)::end_interior_row ! Interior row ending index

      start_interior_row = 2
      end_interior_row = Ny + 1

      if (curr_rank_mpi .EQ. 0) then
         start_interior_row = 3 ! Lower boundary processor
      end if

      if (curr_rank_mpi .EQ. num_ranks_mpi - 1) then
         end_interior_row = Ny ! Upper boundary processor
      end if

   end subroutine Get_Interior_Indices
end module domain
