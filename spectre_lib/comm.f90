module comm
   implicit none

contains
   subroutine Init_MPI(curr_rank_mpi, num_ranks_mpi, ierr_mpi)
      use mpi
      implicit none

      integer, intent(out)::curr_rank_mpi
      integer, intent(out)::num_ranks_mpi
      integer, intent(inout)::ierr_mpi

      call MPI_Init(ierr_mpi)
      call MPI_Comm_Size(MPI_COMM_WORLD, num_ranks_mpi, ierr_mpi)
      call MPI_Comm_Rank(MPI_COMM_WORLD, curr_rank_mpi, ierr_mpi)

   end subroutine Init_MPI

   subroutine Exchange_Halos(f, Ny, Nx, curr_rank_mpi, num_ranks_mpi, ierr_mpi)
      use mpi
      use, intrinsic::iso_c_binding
      use domain
      implicit none
      integer, intent(in)::Ny
      integer, intent(in)::Nx
      complex(C_DOUBLE_COMPLEX), intent(inout), dimension(Ny + 2, Nx)::f
      integer, intent(in)::curr_rank_mpi
      integer, intent(in)::num_ranks_mpi
      integer, intent(inout)::ierr_mpi

      integer::north_rank_mpi
      integer::south_rank_mpi
      integer::tag_up
      integer::tag_down

      integer::start_interior_row
      integer::end_interior_row
      integer::halo_row_type_mpi

      integer::reqs(4)

      call Get_Neighbour_Ranks(curr_rank_mpi, num_ranks_mpi, north_rank_mpi, south_rank_mpi)

      call MPI_Type_Vector(Nx, 1, Ny + 2, MPI_C_DOUBLE_COMPLEX, halo_row_type_mpi, ierr_mpi)
      call MPI_Type_Commit(halo_row_type_mpi, ierr_mpi)

      tag_up = 100
      tag_down = 200

      call Get_Interior_Indices(Ny, curr_rank_mpi, num_ranks_mpi, start_interior_row, end_interior_row)

      ! --- Post receives FIRST so incoming data has a destination ---
      call MPI_Irecv(f(Ny + 2, 1), 1, halo_row_type_mpi, north_rank_mpi, tag_down, &
                     MPI_COMM_WORLD, reqs(1), ierr_mpi)

      call MPI_Irecv(f(1, 1), 1, halo_row_type_mpi, south_rank_mpi, tag_up, &
                     MPI_COMM_WORLD, reqs(2), ierr_mpi)

      ! --- Post sends ---
      call MPI_Isend(f(end_interior_row, 1), 1, halo_row_type_mpi, north_rank_mpi, tag_up, &
                     MPI_COMM_WORLD, reqs(3), ierr_mpi)

      call MPI_Isend(f(start_interior_row, 1), 1, halo_row_type_mpi, south_rank_mpi, tag_down, &
                     MPI_COMM_WORLD, reqs(4), ierr_mpi)

      ! --- Wait for all four to complete before touching the datatype/buffers ---
      call MPI_Waitall(4, reqs, MPI_STATUSES_IGNORE, ierr_mpi)

      call MPI_Type_free(halo_row_type_mpi, ierr_mpi)
   end subroutine Exchange_Halos

   ! ============================================================
! Gather_Column_To_0: collect the column f(2:Ny+1, j) from
! each rank into x_global(1:Ny_g) on rank 0.
!
! Per mode j:
!   - Each rank has Ny local interior values for this mode.
!   - Rank 0 receives them, places each into the correct slot
!     of the global Ny_g-column based on l2g_start.
!
! Input:
!   f_local(Ny+2, Nx)  - local field
!   j                  - which mode (column) to gather
!   Ny                 - local interior rows
!   Ny_g               - global interior rows
!   Nx                 - total columns (unused, kept for signature)
!   curr_rank_mpi      - this rank
!   num_ranks_mpi      - total ranks
!
! Output (rank 0 only):
!   x_global(Ny_g)     - full column for mode j
!
! Other ranks: x_global not used.
! ============================================================
   subroutine Gather_Column_To_0(f_local, j, Ny, Ny_g, Nx, &
                                 curr_rank_mpi, num_ranks_mpi, &
                                 x_global, ierr_mpi)
      use mpi
      use, intrinsic :: iso_c_binding
      use domain
      implicit none

      integer, intent(in) :: j, Ny, Ny_g, Nx, curr_rank_mpi, num_ranks_mpi
      complex(C_DOUBLE_COMPLEX), intent(in), dimension(Ny + 2, Nx) :: f_local
      complex(C_DOUBLE_COMPLEX), intent(out), dimension(Ny_g) :: x_global
      integer, intent(inout) :: ierr_mpi

      complex(C_DOUBLE_COMPLEX), allocatable :: send_buf(:), recv_buf(:)
      integer, allocatable :: counts(:), displs(:)
      integer :: r, ny_r, l2g_start_r, l2g_end_r
      integer :: ierr, req

      ! --- counts / displs ---
      allocate (counts(num_ranks_mpi), displs(num_ranks_mpi))
      do r = 0, num_ranks_mpi - 1
         call Decompose_Domain(Ny_g, num_ranks_mpi, r, ny_r, l2g_start_r, l2g_end_r)
         counts(r + 1) = ny_r
         displs(r + 1) = l2g_start_r - 1
      end do

      ! --- pack this rank's column ---
      allocate (send_buf(Ny))
      send_buf(1:Ny) = f_local(2:Ny + 1, j)

      allocate (recv_buf(Ny_g))
      recv_buf = complex(0.0_C_DOUBLE, 0.0_C_DOUBLE)

      ! --- non-blocking gather ---
      call MPI_Igatherv(send_buf, Ny, MPI_C_DOUBLE_COMPLEX, &
                        recv_buf, counts, displs, MPI_C_DOUBLE_COMPLEX, &
                        0, MPI_COMM_WORLD, req, ierr)

      ! *** Buffers above MUST remain untouched until the wait below. ***
      ! *** This is where you would overlap with local computation.      ***

      call MPI_Wait(req, MPI_STATUS_IGNORE, ierr)

      if (ierr /= MPI_SUCCESS) ierr_mpi = ierr
      if (curr_rank_mpi == 0) x_global = recv_buf

      deallocate (counts, displs, send_buf, recv_buf)

   end subroutine Gather_Column_To_0

   ! ============================================================
! Scatter_Column_From_0: distribute x_global(1:Ny_g) from
! rank 0 back into f_local(2:Ny+1, j) on each rank.
!
! Inverse of Gather_Column_To_0.
!
! Input (rank 0 only):
!   x_global(Ny_g)     - solved column for mode j
!
! Output:
!   f_local(Ny+2, Nx)  - interior rows 2..Ny+1 for column j filled
!
! Other ranks: x_global ignored.
! ============================================================
   subroutine Scatter_Column_From_0(x_global, Ny_g, f_local, j, Ny, Nx, &
                                    curr_rank_mpi, num_ranks_mpi, ierr_mpi)
      use mpi
      use domain
      use, intrinsic :: iso_c_binding
      implicit none

      integer, intent(in) :: j, Ny_g, Ny, Nx, curr_rank_mpi, num_ranks_mpi
      complex(C_DOUBLE_COMPLEX), intent(in), dimension(Ny_g) :: x_global
      complex(C_DOUBLE_COMPLEX), intent(inout), dimension(Ny + 2, Nx) :: f_local
      integer, intent(inout) :: ierr_mpi

      complex(C_DOUBLE_COMPLEX), allocatable :: send_buf(:), recv_buf(:)
      integer, allocatable :: counts(:), displs(:)
      integer :: r, ny_r, l2g_start_r, l2g_end_r
      integer :: ierr, req

      allocate (counts(num_ranks_mpi), displs(num_ranks_mpi))
      do r = 0, num_ranks_mpi - 1
         call Decompose_Domain(Ny_g, num_ranks_mpi, r, ny_r, l2g_start_r, l2g_end_r)
         counts(r + 1) = ny_r
         displs(r + 1) = l2g_start_r - 1
      end do

      ! --- rank 0 holds the full solved column ---
      allocate (send_buf(Ny_g))
      if (curr_rank_mpi == 0) send_buf = x_global

      allocate (recv_buf(Ny))

      call MPI_Iscatterv(send_buf, counts, displs, MPI_C_DOUBLE_COMPLEX, &
                         recv_buf, Ny, MPI_C_DOUBLE_COMPLEX, &
                         0, MPI_COMM_WORLD, req, ierr)

      ! *** Overlap point: do any local work that does not touch recv_buf ***
      ! *** or f_local(2:Ny+1, j) here.                                   ***

      call MPI_Wait(req, MPI_STATUS_IGNORE, ierr)

      if (ierr /= MPI_SUCCESS) ierr_mpi = ierr

      f_local(2:Ny + 1, j) = recv_buf

      deallocate (counts, displs, send_buf, recv_buf)

   end subroutine Scatter_Column_From_0

   ! ============================================================
   ! [TEST3_1 CHANGE] Transpose_To_Modes and Transpose_From_Modes
   !
   ! Why: the tridiagonal system of a Fourier mode needs ALL Ny_g rows of
   ! that mode, but each rank only holds Ny rows of EVERY mode. The original
   ! code therefore gathered one column at a time to rank 0 and solved
   ! every mode there. Here the data are redistributed once per solve with
   ! one MPI_Alltoallv, so that each rank holds all Ny_g rows of its own
   ! share of the modes and can run the (unchanged) Thomas solver on them.
   ! A second MPI_Alltoallv sends the solutions back.
   !
   ! The modes are split over the ranks with the existing Decompose_Domain
   ! (called with Nx instead of Ny_g), exactly like the rows are split:
   ! rank r owns modes m_start..m_end, nm = m_end - m_start + 1 of them.
   !
   ! Transpose_To_Modes  : f_local(Ny+2, Nx)  ->  f_modes(Ny_g, nm)
   ! Transpose_From_Modes: f_modes(Ny_g, nm)  ->  f_local(2:Ny+1, all Nx)
   ! (ghost rows 1 and Ny+2 of f_local are not touched, as before)
   ! ============================================================
   subroutine Transpose_To_Modes(f_local, Ny, Ny_g, Nx, nm, &
                                 curr_rank_mpi, num_ranks_mpi, &
                                 f_modes, ierr_mpi)
      use mpi
      use, intrinsic :: iso_c_binding
      use domain
      implicit none

      integer, intent(in) :: Ny, Ny_g, Nx, nm, curr_rank_mpi, num_ranks_mpi
      complex(C_DOUBLE_COMPLEX), intent(in), dimension(Ny + 2, Nx) :: f_local
      complex(C_DOUBLE_COMPLEX), intent(out), dimension(Ny_g, nm) :: f_modes
      integer, intent(inout) :: ierr_mpi

      complex(C_DOUBLE_COMPLEX), allocatable :: send_buf(:), recv_buf(:)
      integer, allocatable :: scounts(:), sdispls(:), rcounts(:), rdispls(:)
      integer, allocatable :: m_start_r(:), m_end_r(:), row_start_r(:), ny_r(:)
      integer :: r, k, j, pos, nm_r, ms, me, ny_rr, gs, ge, ierr

      allocate (scounts(num_ranks_mpi), sdispls(num_ranks_mpi))
      allocate (rcounts(num_ranks_mpi), rdispls(num_ranks_mpi))
      allocate (m_start_r(num_ranks_mpi), m_end_r(num_ranks_mpi))
      allocate (row_start_r(num_ranks_mpi), ny_r(num_ranks_mpi))

      do r = 0, num_ranks_mpi - 1
         call Decompose_Domain(Nx, num_ranks_mpi, r, nm_r, ms, me)         ! modes owned by rank r
         call Decompose_Domain(Ny_g, num_ranks_mpi, r, ny_rr, gs, ge)      ! rows owned by rank r
         m_start_r(r + 1) = ms
         m_end_r(r + 1) = me
         row_start_r(r + 1) = gs
         ny_r(r + 1) = ny_rr
         scounts(r + 1) = Ny*nm_r          ! my rows of the modes that rank r owns
         rcounts(r + 1) = ny_rr*nm         ! rank r's rows of the modes that I own
      end do
      sdispls(1) = 0
      rdispls(1) = 0
      do r = 2, num_ranks_mpi
         sdispls(r) = sdispls(r - 1) + scounts(r - 1)
         rdispls(r) = rdispls(r - 1) + rcounts(r - 1)
      end do

      allocate (send_buf(sum(scounts)), recv_buf(sum(rcounts)))

      ! --- pack: for each destination rank, my interior rows of its modes, column by column ---
      pos = 0
      do r = 1, num_ranks_mpi
         do j = m_start_r(r), m_end_r(r)
            send_buf(pos + 1:pos + Ny) = f_local(2:Ny + 1, j)
            pos = pos + Ny
         end do
      end do

      call MPI_Alltoallv(send_buf, scounts, sdispls, MPI_C_DOUBLE_COMPLEX, &
                         recv_buf, rcounts, rdispls, MPI_C_DOUBLE_COMPLEX, &
                         MPI_COMM_WORLD, ierr)
      if (ierr /= MPI_SUCCESS) ierr_mpi = ierr

      ! --- unpack: rows of every rank, placed at their global row numbers ---
      pos = 0
      do r = 1, num_ranks_mpi
         do k = 1, nm
            f_modes(row_start_r(r):row_start_r(r) + ny_r(r) - 1, k) = recv_buf(pos + 1:pos + ny_r(r))
            pos = pos + ny_r(r)
         end do
      end do

      deallocate (scounts, sdispls, rcounts, rdispls, m_start_r, m_end_r, row_start_r, ny_r)
      deallocate (send_buf, recv_buf)

   end subroutine Transpose_To_Modes

   subroutine Transpose_From_Modes(f_modes, Ny, Ny_g, Nx, nm, &
                                   curr_rank_mpi, num_ranks_mpi, &
                                   f_local, ierr_mpi)
      use mpi
      use, intrinsic :: iso_c_binding
      use domain
      implicit none

      integer, intent(in) :: Ny, Ny_g, Nx, nm, curr_rank_mpi, num_ranks_mpi
      complex(C_DOUBLE_COMPLEX), intent(in), dimension(Ny_g, nm) :: f_modes
      complex(C_DOUBLE_COMPLEX), intent(inout), dimension(Ny + 2, Nx) :: f_local
      integer, intent(inout) :: ierr_mpi

      complex(C_DOUBLE_COMPLEX), allocatable :: send_buf(:), recv_buf(:)
      integer, allocatable :: scounts(:), sdispls(:), rcounts(:), rdispls(:)
      integer, allocatable :: m_start_r(:), nm_all(:), row_start_r(:), ny_r(:)
      integer :: r, k, pos, nm_r, ms, me, ny_rr, gs, ge, ierr

      allocate (scounts(num_ranks_mpi), sdispls(num_ranks_mpi))
      allocate (rcounts(num_ranks_mpi), rdispls(num_ranks_mpi))
      allocate (m_start_r(num_ranks_mpi), nm_all(num_ranks_mpi))
      allocate (row_start_r(num_ranks_mpi), ny_r(num_ranks_mpi))

      do r = 0, num_ranks_mpi - 1
         call Decompose_Domain(Nx, num_ranks_mpi, r, nm_r, ms, me)
         call Decompose_Domain(Ny_g, num_ranks_mpi, r, ny_rr, gs, ge)
         m_start_r(r + 1) = ms
         nm_all(r + 1) = nm_r
         row_start_r(r + 1) = gs
         ny_r(r + 1) = ny_rr
         scounts(r + 1) = ny_rr*nm         ! rank r's rows of the modes that I own
         rcounts(r + 1) = Ny*nm_r          ! my rows of the modes that rank r owns
      end do
      sdispls(1) = 0
      rdispls(1) = 0
      do r = 2, num_ranks_mpi
         sdispls(r) = sdispls(r - 1) + scounts(r - 1)
         rdispls(r) = rdispls(r - 1) + rcounts(r - 1)
      end do

      allocate (send_buf(sum(scounts)), recv_buf(sum(rcounts)))

      ! --- pack: for each destination rank, ITS rows of my modes, column by column ---
      pos = 0
      do r = 1, num_ranks_mpi
         do k = 1, nm
            send_buf(pos + 1:pos + ny_r(r)) = f_modes(row_start_r(r):row_start_r(r) + ny_r(r) - 1, k)
            pos = pos + ny_r(r)
         end do
      end do

      call MPI_Alltoallv(send_buf, scounts, sdispls, MPI_C_DOUBLE_COMPLEX, &
                         recv_buf, rcounts, rdispls, MPI_C_DOUBLE_COMPLEX, &
                         MPI_COMM_WORLD, ierr)
      if (ierr /= MPI_SUCCESS) ierr_mpi = ierr

      ! --- unpack: my interior rows of every mode, from the rank that solved it ---
      pos = 0
      do r = 1, num_ranks_mpi
         do k = 1, nm_all(r)
            f_local(2:Ny + 1, m_start_r(r) + k - 1) = recv_buf(pos + 1:pos + Ny)
            pos = pos + Ny
         end do
      end do

      deallocate (scounts, sdispls, rcounts, rdispls, m_start_r, nm_all, row_start_r, ny_r)
      deallocate (send_buf, recv_buf)

   end subroutine Transpose_From_Modes
end module comm
