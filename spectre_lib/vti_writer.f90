module vti_writer

contains

   subroutine Write_VTI(u_s, v_s, p_s, step, Ny, Nx, Ny_g, dy, Lx, &
                        rank, nranks, backward_fft_plan)
      use, intrinsic :: iso_c_binding
      use mpi
      use fft
      use comm
      implicit none

      integer, intent(in) :: step, Ny, Nx, Ny_g, rank, nranks
      real(C_DOUBLE), intent(in) :: dy, Lx
      type(c_ptr), intent(in) :: backward_fft_plan
      complex(C_DOUBLE_COMPLEX), dimension(Ny + 2, Nx), intent(in) :: u_s, v_s, p_s

      complex(C_DOUBLE_COMPLEX), allocatable :: u_r_local(:, :), v_r_local(:, :), p_r_local(:, :)
      complex(C_DOUBLE_COMPLEX), allocatable :: u_col(:), v_col(:), p_col(:)
      complex(C_DOUBLE_COMPLEX), allocatable :: u_r_glob(:, :), v_r_glob(:, :), p_r_glob(:, :)
      integer :: i, j, ierr, fu
      character(len=20), parameter :: folder_name = "output"
      character(len=64) :: fname

      allocate (u_r_local(Ny + 2, Nx), v_r_local(Ny + 2, Nx), p_r_local(Ny + 2, Nx))

      !-------------------------------------------------------------------
      ! Fourier transform to real space (1D along x, per-row, unnormalized)
      !-------------------------------------------------------------------
      call Fourier_Transform(u_s, Ny, Nx, backward_fft_plan, u_r_local)
      call Fourier_Transform(v_s, Ny, Nx, backward_fft_plan, v_r_local)
      call Fourier_Transform(p_s, Ny, Nx, backward_fft_plan, p_r_local)

      !-------------------------------------------------------------------
      ! Normalize — same pattern as Solver_Diagnostics (divide by Nx)
      !-------------------------------------------------------------------
      u_r_local = u_r_local/real(Nx, C_DOUBLE)
      v_r_local = v_r_local/real(Nx, C_DOUBLE)
      p_r_local = p_r_local/real(Nx, C_DOUBLE)

      allocate (u_col(Ny_g), v_col(Ny_g), p_col(Ny_g))
      if (rank == 0) then
         allocate (u_r_glob(Ny_g, Nx), v_r_glob(Ny_g, Nx), p_r_glob(Ny_g, Nx))
      else
         allocate (u_r_glob(1, 1), v_r_glob(1, 1), p_r_glob(1, 1))
      end if

      do j = 1, Nx
         call Gather_Column_To_0(u_r_local, j, Ny, Ny_g, Nx, rank, nranks, u_col, ierr)
         call Gather_Column_To_0(v_r_local, j, Ny, Ny_g, Nx, rank, nranks, v_col, ierr)
         call Gather_Column_To_0(p_r_local, j, Ny, Ny_g, Nx, rank, nranks, p_col, ierr)
         if (rank == 0) then
            u_r_glob(:, j) = u_col
            v_r_glob(:, j) = v_col
            p_r_glob(:, j) = p_col
         end if
      end do

      if (rank == 0) then
         write (fname, '(A,I6.6,A)') 'flow_field_', step, '.vti'
         open (newunit=fu, file=trim(folder_name)//"/"//trim(fname), status='replace', action='write')

         write (fu, '(A)') '<?xml version="1.0"?>'
         write (fu, '(A)') '<VTKFile type="ImageData" version="0.1" byte_order="LittleEndian">'
         write (fu, '(A,I0,A,I0,A,A)') &
            '  <ImageData WholeExtent="0 ', Nx - 1, ' 0 ', Ny_g - 1, ' 0 0"', &
            ' Origin="0 0 0"'
         write (fu, '(A,ES16.8,A,ES16.8,A)') &
            '   Spacing="', Lx/real(Nx, C_DOUBLE), ' ', dy, ' 1.0">'
         write (fu, '(A,I0,A,I0,A)') &
            '    <Piece Extent="0 ', Nx - 1, ' 0 ', Ny_g - 1, ' 0 0">'
         write (fu, '(A)') '      <PointData Vectors="velocity" Scalars="u">'

         write (fu, '(A)') '        <DataArray type="Float32" Name="u" NumberOfComponents="1" format="ascii">'
         do i = 1, Ny_g
            do j = 1, Nx
               write (fu, '(A,E16.8)') '          ', real(u_r_glob(i, j), C_DOUBLE)
            end do
         end do
         write (fu, '(A)') '        </DataArray>'

         write (fu, '(A)') '        <DataArray type="Float32" Name="v" NumberOfComponents="1" format="ascii">'
         do i = 1, Ny_g
            do j = 1, Nx
               write (fu, '(A,E16.8)') '          ', real(v_r_glob(i, j), C_DOUBLE)
            end do
         end do
         write (fu, '(A)') '        </DataArray>'

         write (fu, '(A)') '        <DataArray type="Float32" Name="p" NumberOfComponents="1" format="ascii">'
         do i = 1, Ny_g
            do j = 1, Nx
               write (fu, '(A,E16.8)') '          ', real(p_r_glob(i, j), C_DOUBLE)
            end do
         end do
         write (fu, '(A)') '        </DataArray>'

         write (fu, '(A)') '        <DataArray type="Float32" Name="velocity" NumberOfComponents="3" format="ascii">'
         do i = 1, Ny_g
            do j = 1, Nx
               write (fu, '(A,3E16.8)') '          ', &
                  real(u_r_glob(i, j), C_DOUBLE), &
                  real(v_r_glob(i, j), C_DOUBLE), &
                  0.0_C_DOUBLE
            end do
         end do
         write (fu, '(A)') '        </DataArray>'

         write (fu, '(A)') '      </PointData>'
         write (fu, '(A)') '    </Piece>'
         write (fu, '(A)') '  </ImageData>'
         write (fu, '(A)') '</VTKFile>'
         close (fu)

         write (*, '(A,A)') ' Wrote ', trim(fname)
      end if

      deallocate (u_r_local, v_r_local, p_r_local, &
                  u_col, v_col, p_col, &
                  u_r_glob, v_r_glob, p_r_glob)

   end subroutine Write_VTI

end module vti_writer
