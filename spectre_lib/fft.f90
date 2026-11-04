module fft
   implicit none

contains
   subroutine Forward_Plan(Ny, Nx, fft_plan)
      use, intrinsic::iso_c_binding
      implicit none
      include 'fftw3.f03'

      integer, intent(in)::Ny
      integer, intent(in)::Nx
      type(c_ptr), intent(out)::fft_plan

      complex(c_double_complex), allocatable, dimension(:, :)::f_in
      complex(c_double_complex), allocatable, dimension(:, :)::f_out
      integer(c_int) :: rank
      integer(c_int), dimension(1) :: n_dims
      integer(c_int) :: howmany, idist, odist, istride, ostride

      allocate (f_in(Ny + 2, Nx))
      allocate (f_out(Ny + 2, Nx))

      rank = 1
      n_dims(1) = Nx
      howmany = Ny + 2
      istride = Ny + 2
      ostride = Ny + 2
      idist = 1
      odist = 1

      ! Batch FFT Plan
      fft_plan = fftw_plan_many_dft(rank, n_dims, howmany, &
                                    f_in, n_dims, istride, idist, &
                                    f_out, n_dims, ostride, odist, &
                                    FFTW_FORWARD, FFTW_ESTIMATE)
      deallocate (f_in)
      deallocate (f_out)
   end subroutine Forward_Plan

   subroutine Backward_Plan(Ny, Nx, fft_plan)
      use, intrinsic::iso_c_binding
      implicit none
      include 'fftw3.f03'

      integer, intent(in)::Ny
      integer, intent(in)::Nx
      type(c_ptr), intent(out)::fft_plan

      complex(c_double_complex), allocatable, dimension(:, :)::f_in
      complex(c_double_complex), allocatable, dimension(:, :)::f_out
      integer(c_int) :: rank
      integer(c_int), dimension(1) :: n_dims
      integer(c_int) :: howmany, idist, odist, istride, ostride

      allocate (f_in(Ny + 2, Nx))
      allocate (f_out(Ny + 2, Nx))

      rank = 1
      n_dims(1) = Nx
      howmany = Ny + 2
      istride = Ny + 2
      ostride = Ny + 2
      idist = 1
      odist = 1

      ! Batch FFT Plan
      fft_plan = fftw_plan_many_dft(rank, n_dims, howmany, &
                                    f_in, n_dims, istride, idist, &
                                    f_out, n_dims, ostride, odist, &
                                    FFTW_BACKWARD, FFTW_ESTIMATE)
      deallocate (f_in)
      deallocate (f_out)
   end subroutine Backward_Plan

   subroutine Fourier_Transform(f_in, Ny, Nx, fft_plan, f_out)
      use, intrinsic::iso_c_binding
      implicit none
      include 'fftw3.f03'

      integer, intent(in)::Ny
      integer, intent(in)::Nx
      complex(c_double_complex), intent(in)::f_in(Ny + 2, Nx)
      complex(c_double_complex), intent(out)::f_out(Ny + 2, Nx)
      type(c_ptr), intent(in)::fft_plan

      complex(c_double_complex), allocatable, dimension(:, :)::f_scratch

      allocate (f_scratch(Ny + 2, Nx))
      f_scratch = f_in

      call fftw_execute_dft(fft_plan, f_scratch, f_out)

      deallocate (f_scratch)

   end subroutine

   subroutine Wavelengths(Nx, Lx, kx)
      use, intrinsic::iso_c_binding
      implicit none
      integer, intent(in)::Nx
      real(c_double), intent(in)::Lx
      real(c_double), intent(out), dimension(Nx)::kx

      real(c_double), parameter::PI = acos(-1.0_C_DOUBLE)

      integer::i

      real(c_double)::factor
      factor = 2*PI/Lx
      do i = 1, Nx/2 + 1
         kx(i) = factor*real(i - 1, c_double)
      end do

      do i = Nx/2 + 2, Nx
         kx(i) = -real(Nx - i + 1, c_double)*factor
      end do

   end subroutine Wavelengths

   subroutine Dealias(f, Ny, Nx, Lx, kx)

      use, intrinsic :: iso_c_binding
      implicit none

      integer, intent(in) :: Ny, Nx

      real(c_double), intent(in) :: Lx
      real(c_double), intent(in) :: kx(Nx)

      complex(c_double_complex), intent(inout) :: f(Ny + 2, Nx)

      integer :: i
      real(c_double) :: kx_max
      real(c_double), parameter :: PI = acos(-1.0_c_double)

      ! 2/3-rule cutoff
      kx_max = (2.0_c_double*PI/Lx)*real(Nx, c_double)/3.0_c_double

      do i = 1, Nx

         if (abs(kx(i)) > kx_max) then

            f(:, i) = cmplx(0.0_c_double, 0.0_c_double, &
                            kind=c_double)

         end if

      end do

   end subroutine Dealias

   subroutine Dealias_Strong(f, Ny, Nx, Lx, kx)

      use, intrinsic :: iso_c_binding
      implicit none

      integer, intent(in) :: Ny, Nx

      real(c_double), intent(in) :: Lx
      real(c_double), intent(in) :: kx(Nx)

      complex(c_double_complex), intent(inout) :: f(Ny + 2, Nx)

      integer :: i
      real(c_double) :: kx_max
      real(c_double), parameter :: PI = acos(-1.0_c_double)

      ! 1/3-rule cutoff
      kx_max = (1.0_c_double*PI/Lx)*real(Nx, c_double)/3.0_c_double

      do i = 1, Nx

         if (abs(kx(i)) > kx_max) then

            f(:, i) = cmplx(0.0_c_double, 0.0_c_double, &
                            kind=c_double)

         end if

      end do

   end subroutine Dealias_Strong

   subroutine Dealias_Stronger(f, Ny, Nx, Lx, kx)

      use, intrinsic :: iso_c_binding
      implicit none

      integer, intent(in) :: Ny, Nx

      real(c_double), intent(in) :: Lx
      real(c_double), intent(in) :: kx(Nx)

      complex(c_double_complex), intent(inout) :: f(Ny + 2, Nx)

      integer :: i
      real(c_double) :: kx_max
      real(c_double), parameter :: PI = acos(-1.0_c_double)

      ! 1/6-rule cutoff
      kx_max = (1.0_c_double*PI/Lx)*real(Nx, c_double)/6.0_c_double

      do i = 1, Nx

         if (abs(kx(i)) > kx_max) then

            f(:, i) = cmplx(0.0_c_double, 0.0_c_double, &
                            kind=c_double)

         end if

      end do

   end subroutine Dealias_Stronger

end module fft
