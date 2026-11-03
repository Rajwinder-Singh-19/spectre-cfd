module field
   implicit none
   ! All input fields are in fourier space unless specified otherwise
contains
   subroutine Dx_Spectral(f, Ny, Nx, kx, curr_rank_mpi, num_ranks_mpi, df_dx)
      use, intrinsic::iso_c_binding
      use fft
      use domain
      implicit none
      integer, intent(in)::Ny
      integer, intent(in)::Nx
      integer, intent(in)::curr_rank_mpi
      integer, intent(in)::num_ranks_mpi
      real(c_double), intent(in), dimension(Nx)::kx
      complex(c_double_complex), dimension(Ny + 2, Nx), intent(in)::f
      complex(c_double_complex), dimension(Ny + 2, Nx), intent(out)::df_dx

      complex(c_double_complex), parameter::imag = complex(0.0_c_double, 1.0_c_double)
      integer::i, j, j_start, j_end

      call Get_Interior_Indices(Ny, curr_rank_mpi, num_ranks_mpi, j_start, j_end)

      df_dx = cmplx(0.0_c_double, 0.0_c_double, kind=c_double)
      do j = j_start, j_end
         do i = 1, Nx
            df_dx(j, i) = imag*kx(i)*f(j, i)
         end do
      end do

   end subroutine Dx_Spectral

   subroutine Dy_Forward(f, Ny, Nx, dy, curr_rank_mpi, num_ranks_mpi, df_dy)
      use, intrinsic::iso_c_binding
      use domain
      implicit none
      integer, intent(in)::Ny
      integer, intent(in)::Nx
      integer, intent(in)::curr_rank_mpi
      integer, intent(in)::num_ranks_mpi
      real(c_double), intent(in)::dy
      complex(c_double_complex), dimension(Ny + 2, Nx), intent(in)::f
      complex(c_double_complex), dimension(Ny + 2, Nx), intent(out)::df_dy

      integer::i, j, j_start, j_end
      real(c_double)::inv_dy

      call Get_Interior_Indices(Ny, curr_rank_mpi, num_ranks_mpi, j_start, j_end)

      inv_dy = 1.0_c_double/dy

      df_dy = cmplx(0.0_c_double, 0.0_c_double, kind=c_double)
      do j = j_start, j_end
         do i = 1, Nx
            df_dy(j, i) = (f(j + 1, i) - f(j, i))*inv_dy
         end do
      end do

   end subroutine Dy_Forward

   subroutine Dy_Backward(f, Ny, Nx, dy, curr_rank_mpi, num_ranks_mpi, df_dy)
      use, intrinsic::iso_c_binding
      use domain
      implicit none
      integer, intent(in)::Ny
      integer, intent(in)::Nx
      integer, intent(in)::curr_rank_mpi
      integer, intent(in)::num_ranks_mpi
      real(c_double), intent(in)::dy
      complex(c_double_complex), dimension(Ny + 2, Nx), intent(in)::f
      complex(c_double_complex), dimension(Ny + 2, Nx), intent(out)::df_dy

      integer::i, j, j_start, j_end
      real(c_double)::inv_dy

      call Get_Interior_Indices(Ny, curr_rank_mpi, num_ranks_mpi, j_start, j_end)

      inv_dy = 1.0_c_double/dy

      df_dy = cmplx(0.0_c_double, 0.0_c_double, kind=c_double)
      do j = j_start, j_end
         do i = 1, Nx
            df_dy(j, i) = (f(j, i) - f(j - 1, i))*inv_dy
         end do
      end do

   end subroutine Dy_Backward

   subroutine Dy_Center(f, Ny, Nx, dy, curr_rank_mpi, num_ranks_mpi, df_dy)
      use, intrinsic::iso_c_binding
      use domain
      implicit none
      integer, intent(in)::Ny
      integer, intent(in)::Nx
      integer, intent(in)::curr_rank_mpi
      integer, intent(in)::num_ranks_mpi
      real(c_double), intent(in)::dy
      complex(c_double_complex), dimension(Ny + 2, Nx), intent(in)::f
      complex(c_double_complex), dimension(Ny + 2, Nx), intent(out)::df_dy

      integer::i, j, j_start, j_end
      real(c_double)::inv_dy

      call Get_Interior_Indices(Ny, curr_rank_mpi, num_ranks_mpi, j_start, j_end)

      inv_dy = 1.0_c_double/dy

      df_dy = cmplx(0.0_c_double, 0.0_c_double, kind=c_double)
      do j = j_start, j_end
         do i = 1, Nx
            df_dy(j, i) = (f(j + 1, i) - f(j - 1, i))*inv_dy*0.5_c_double
         end do
      end do

   end subroutine Dy_Center

   subroutine Divergence(f, g, Ny, Nx, kx, dy, curr_rank_mpi, num_ranks_mpi, div_fg)
      use, intrinsic::iso_c_binding
      implicit none
      integer, intent(in)::Ny
      integer, intent(in)::Nx
      integer, intent(in)::curr_rank_mpi
      integer, intent(in)::num_ranks_mpi
      real(c_double), intent(in), dimension(Nx)::kx
      real(c_double), intent(in)::dy
      complex(c_double_complex), dimension(Ny + 2, Nx), intent(in)::f
      complex(c_double_complex), dimension(Ny + 2, Nx), intent(in)::g
      complex(c_double_complex), dimension(Ny + 2, Nx), intent(out)::div_fg
      complex(c_double_complex), allocatable, dimension(:, :)::dg_dy

      allocate (dg_dy(Ny + 2, Nx))

      div_fg = cmplx(0.0_c_double, 0.0_c_double, kind=c_double)

      call Dx_Spectral(f, Ny, Nx, kx, curr_rank_mpi, num_ranks_mpi, div_fg)
      call Dy_Center(g, Ny, Nx, dy, curr_rank_mpi, num_ranks_mpi, dg_dy)

      div_fg = div_fg + dg_dy

      deallocate (dg_dy)

   end subroutine Divergence

   subroutine Rhie_Chow_Face_Velocity(v, p, Ny, Nx, kx, nu, dy, dt, curr_rank_mpi, num_ranks_mpi, ierr_mpi, v_face)
      use, intrinsic::iso_c_binding
      use domain
      use comm
      implicit none
      integer, intent(in)::Ny
      integer, intent(in)::Nx
      real(c_double), intent(in)::dy
      real(c_double), intent(in)::dt
      real(c_double), intent(in)::kx(Nx)
      real(c_double), intent(in)::nu(Nx)
      integer, intent(in)::curr_rank_mpi, num_ranks_mpi
      integer, intent(inout)::ierr_mpi
      complex(c_double_complex), intent(inout), dimension(Ny + 2, Nx)::v
      complex(c_double_complex), intent(inout), dimension(Ny + 2, Nx)::p
      complex(c_double_complex), intent(out), dimension(Ny + 2, Nx)::v_face

      real(c_double), dimension(Nx)::a
      real(c_double)::inv_d2

      integer::start_j, end_j, i, j

      call Get_Interior_Indices(Ny, curr_rank_mpi, num_ranks_mpi, start_j, end_j)

      call Exchange_Halos(v, Ny, nx, curr_rank_mpi, num_ranks_mpi, ierr_mpi)
      call Exchange_Halos(p, Ny, nx, curr_rank_mpi, num_ranks_mpi, ierr_mpi)

      v_face = cmplx(0.0_c_double, 0.0_c_double, c_double_complex)
      inv_d2 = 1.0_c_double/(dy*dy)
      !a(1) = 1.0_c_double + nu(1)*dt*(3.0_c_double*inv_d2 + kx(i)*kx(i))
      !a(Nx) = 1.0_c_double + nu(Nx)*dt*(3.0_c_double*inv_d2 + kx(i)*kx(i))
      do i = 1, Nx
         a(i) = 1.0_c_double + nu(i)*dt*(2.0_c_double*inv_d2 + kx(i)*kx(i))
      end do

      do j = start_j, end_j
         do i = 1, Nx
            v_face(j, i) = (v(j, i) + v(j + 1, i))*0.5_c_double &
                           - (dt/a(i))*(p(j + 1, i) - p(j, i))/dy
         end do
      end do

      if (curr_rank_mpi == 0) then
         v_face(2, :) = cmplx(0.0_c_double, 0.0_c_double, c_double)
      end if

      if (curr_rank_mpi == num_ranks_mpi - 1) then
         v_face(Ny + 1, :) = cmplx(0.0_c_double, 0.0_c_double, c_double)
      end if

   end subroutine Rhie_Chow_Face_Velocity

end module field
