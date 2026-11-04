module turbulence
   implicit none

contains
   subroutine Exponential_SVV(kxi, epsilon_svv, k_cut, k_Nyq, Q_k, nu_svv)
      use, intrinsic :: iso_c_binding
      implicit none
      real(C_DOUBLE), intent(in)  :: kxi           ! x-wavenumber for this mode
      real(C_DOUBLE), intent(in)  :: epsilon_svv   ! SVV amplitude
      real(C_DOUBLE), intent(in)  :: k_cut         ! SVV cutoff wavenumber
      real(C_DOUBLE), intent(in)  :: k_Nyq         ! Nyquist wavenumber

      real(C_DOUBLE), intent(out) :: Q_k           ! SVV kernel value at this mode
      real(C_DOUBLE), intent(out) :: nu_svv        ! total viscosity: nu + eps*Q

      real(C_DOUBLE) :: arg

      if (abs(kxi) > k_cut) then
         arg = (k_Nyq - abs(kxi))/(k_Nyq - k_cut)
         Q_k = exp(-arg*arg)
      else
         Q_k = 0.0_C_DOUBLE
      end if

      nu_svv = epsilon_svv*Q_k

   end subroutine Exponential_SVV

   subroutine Kolmogorov_Sharp_SVV(kxi, epsilon_svv, k_cut, k_Nyq, Q_k, nu_svv)
      use, intrinsic :: iso_c_binding
      implicit none
      real(C_DOUBLE), intent(in)  :: kxi           ! x-wavenumber for this mode
      real(C_DOUBLE), intent(in)  :: epsilon_svv   ! SVV amplitude
      real(C_DOUBLE), intent(in)  :: k_cut         ! SVV cutoff wavenumber
      real(C_DOUBLE), intent(in)  :: k_Nyq         ! Nyquist wavenumber

      real(C_DOUBLE), intent(out) :: Q_k           ! SVV kernel value at this mode
      real(C_DOUBLE), intent(out) :: nu_svv        ! total viscosity: nu + eps*Q

      real(C_DOUBLE) :: abs_k

      abs_k = abs(kxi)

      if (abs_k <= k_cut) then
         Q_k = 0.0_C_DOUBLE
      else
         ! Kolmogorov k^-4/3 decay, normalized to 1 at k_Nyq
         Q_k = (k_Nyq/abs_k)**(4.0_C_DOUBLE/3.0_C_DOUBLE)
      end if

      nu_svv = epsilon_svv*Q_k

   end subroutine Kolmogorov_Sharp_SVV

   subroutine Kolmogorov_Smooth_SVV(kxi, epsilon_svv, k_cut, k_Nyq, Q_k, nu_svv)
      use, intrinsic :: iso_c_binding
      implicit none
      real(C_DOUBLE), intent(in)  :: kxi           ! x-wavenumber for this mode
      real(C_DOUBLE), intent(in)  :: epsilon_svv   ! SVV amplitude
      real(C_DOUBLE), intent(in)  :: k_cut         ! SVV cutoff wavenumber
      real(C_DOUBLE), intent(in)  :: k_Nyq         ! Nyquist wavenumber

      real(C_DOUBLE), intent(out) :: Q_k           ! SVV kernel value at this mode
      real(C_DOUBLE), intent(out) :: nu_svv        ! total viscosity: nu + eps*Q

      real(C_DOUBLE) :: abs_k, power_part, taper

      abs_k = abs(kxi)

      if (abs_k <= k_cut) then
         Q_k = 0.0_C_DOUBLE
      else
         ! Kolmogorov power-law part
         power_part = (k_Nyq/abs_k)**(4.0_C_DOUBLE/3.0_C_DOUBLE)
         ! Smooth exponential taper from 0 at k_cut to 1 at k_Nyq
         taper = 1.0_C_DOUBLE - exp(-((abs_k - k_cut)/(k_Nyq - k_cut))**2)
         Q_k = power_part*taper
      end if

      nu_svv = epsilon_svv*Q_k

   end subroutine Kolmogorov_Smooth_SVV

   subroutine Plateau_SVV(kxi, epsilon_svv, k_cut, Q_k, nu_svv)
      use, intrinsic :: iso_c_binding
      implicit none
      real(C_DOUBLE), intent(in)  :: kxi, epsilon_svv, k_cut
      real(C_DOUBLE), intent(out) :: Q_k, nu_svv
      real(C_DOUBLE) :: abs_k

      abs_k = abs(kxi)
      if (abs_k <= k_cut) then
         Q_k = 0.0_C_DOUBLE
      else
         Q_k = 1.0_C_DOUBLE   ! Constant dissipation above cutoff
      end if
      nu_svv = epsilon_svv*Q_k
   end subroutine Plateau_SVV

   subroutine Power_2_SVV(kxi, epsilon_svv, k_cut, k_Nyq, Q_k, nu_svv)
      use, intrinsic :: iso_c_binding
      implicit none
      real(C_DOUBLE), intent(in)  :: kxi, epsilon_svv, k_cut, k_Nyq
      real(C_DOUBLE), intent(out) :: Q_k, nu_svv
      real(C_DOUBLE) :: abs_k, arg

      abs_k = abs(kxi)
      if (abs_k <= k_cut) then
         Q_k = 0.0_C_DOUBLE
      else
         arg = (abs_k - k_cut)/(k_Nyq - k_cut)
         Q_k = arg**2
      end if
      nu_svv = epsilon_svv*Q_k
   end subroutine Power_2_SVV

   subroutine Power_3_SVV(kxi, epsilon_svv, k_cut, k_Nyq, Q_k, nu_svv)
      use, intrinsic :: iso_c_binding
      implicit none
      real(C_DOUBLE), intent(in)  :: kxi, epsilon_svv, k_cut, k_Nyq
      real(C_DOUBLE), intent(out) :: Q_k, nu_svv
      real(C_DOUBLE) :: abs_k, arg

      abs_k = abs(kxi)
      if (abs_k <= k_cut) then
         Q_k = 0.0_C_DOUBLE
      else
         arg = (abs_k - k_cut)/(k_Nyq - k_cut)
         Q_k = arg**3
      end if
      nu_svv = epsilon_svv*Q_k
   end subroutine Power_3_SVV

   subroutine Power_4_SVV(kxi, epsilon_svv, k_cut, k_Nyq, Q_k, nu_svv)
      use, intrinsic :: iso_c_binding
      implicit none
      real(C_DOUBLE), intent(in)  :: kxi, epsilon_svv, k_cut, k_Nyq
      real(C_DOUBLE), intent(out) :: Q_k, nu_svv
      real(C_DOUBLE) :: abs_k, arg

      abs_k = abs(kxi)
      if (abs_k <= k_cut) then
         Q_k = 0.0_C_DOUBLE
      else
         arg = (abs_k - k_cut)/(k_Nyq - k_cut)
         Q_k = arg**4
      end if
      nu_svv = epsilon_svv*Q_k
   end subroutine Power_4_SVV

   subroutine Chollet_Lesieur_Viscosity(u, v, Ny, Nx, dy, kx, kc, C_nu, &
                              ierr_mpi, nu_sgs)

      use, intrinsic :: iso_c_binding
      use mpi
      implicit none

      integer, intent(in) :: Ny, Nx
      integer, intent(inout) :: ierr_mpi

      real(c_double), intent(in) :: dy, kc, C_nu
      real(c_double), intent(in) :: kx(Nx)

      complex(c_double_complex), intent(in) :: u(Ny + 2, Nx)
      complex(c_double_complex), intent(in) :: v(Ny + 2, Nx)

      real(c_double), intent(out) :: nu_sgs(Nx)

      ! Local variables
      real(c_double), allocatable :: E_local(:)
      real(c_double), allocatable :: E_k(:)

      real(c_double) :: E_kc
      real(c_double) :: x
      real(c_double) :: f_k
      real(c_double) :: abs_k

      integer :: i, j
      integer :: kc_index

      !--------------------------------------------------
      ! Allocate energy spectrum
      !--------------------------------------------------

      allocate (E_local(Nx))
      allocate (E_k(Nx))

      E_local = 0.0_c_double
      E_k = 0.0_c_double
      nu_sgs = 0.0_c_double

      !--------------------------------------------------
      ! Calculate local energy spectrum
      !--------------------------------------------------

      do i = 1, Nx
         do j = 2, Ny + 1

            E_local(i) = E_local(i) + 0.5_c_double* &
                         (abs(u(j, i))**2 + abs(v(j, i))**2)*dy

         end do
      end do

      !--------------------------------------------------
      ! Global energy spectrum
      !--------------------------------------------------

      call MPI_Allreduce(E_local, E_k, Nx, &
                         MPI_DOUBLE_PRECISION, MPI_SUM, &
                         MPI_COMM_WORLD, ierr_mpi)

      !--------------------------------------------------
      ! Find energy at cutoff wavenumber
      !--------------------------------------------------

      kc_index = minloc(abs(abs(kx) - kc), dim=1)

      E_kc = E_k(kc_index)

      !--------------------------------------------------
      ! Calculate Chollet-Lesieur eddy viscosity
      !--------------------------------------------------

      if (kc <= 0.0_c_double .or. E_kc <= 0.0_c_double) then

         nu_sgs = 0.0_c_double

      else

         do i = 1, Nx

            abs_k = abs(kx(i))

            if (abs_k <= kc) then

               nu_sgs(i) = 0.0_c_double

            else

               x = abs_k/kc

               f_k = sqrt(1.0_c_double + &
                          34.5_c_double*exp(-3.03_c_double/x))

               nu_sgs(i) = C_nu* &
                           sqrt(E_kc)/(kc**1.5_c_double)*f_k

            end if

         end do

      end if

      !--------------------------------------------------
      ! Cleanup
      !--------------------------------------------------

      deallocate (E_local, E_k)

   end subroutine Chollet_Lesieur_Viscosity

   subroutine Seed_Perturbation(u, v, Ny, Nx, Ly, dy, kx, &
                                start_row, curr_rank_mpi, num_ranks_mpi, &
                                ierr_mpi, eps_pert, k_mode_seed)
      use, intrinsic :: iso_c_binding
      use mpi
      use comm
      implicit none

      integer, intent(in) :: Ny, Nx, start_row, curr_rank_mpi, num_ranks_mpi
      integer, intent(in) :: k_mode_seed
      integer, intent(inout) :: ierr_mpi
      real(c_double), intent(in) :: Ly, dy, eps_pert
      real(c_double), intent(in) :: kx(Nx)
      complex(c_double_complex), intent(inout) :: u(Ny + 2, Nx), v(Ny + 2, Nx)

      real(c_double), parameter :: PI = acos(-1.0_c_double)
      real(c_double) :: y, fy, dfdy, kx_phys, inv2Nx
      integer :: jj, k_minus

      if (k_mode_seed < 2 .or. k_mode_seed > Nx/2) then
         if (curr_rank_mpi == 0) print *, 'Seed: k_mode_seed out of range'
         ierr_mpi = 1
         return
      end if

      kx_phys = kx(k_mode_seed)
      k_minus = Nx + 2 - k_mode_seed   ! Hermitian mirror: -k_mode
      inv2Nx = 0.5_c_double*real(Nx, c_double)   ! Nx/2  scaling

      !  psi(x,y) = eps * cos(k*x) * sin^2(pi*y/Ly)
      !  u =  dpsi/dy = eps * cos(kx) * sin(2pi y/Ly)*(pi/Ly)
      !  v = -dpsi/dx = eps * k * sin(kx) * sin^2(pi y/Ly)
      !
      !  spectral (FFTW convention):
      !      û(+k) = (Nx/2) * eps * dfdy             (real)
      !      v̂(+k) = (Nx/2) * eps * (-i) * k * fy    (imag, neg)
      !      û(-k) = conj(û(+k)), v̂(-k) = conj(v̂(+k))

      do jj = 2, Ny + 1
         y = (real(start_row + jj - 2, c_double) + 0.5_c_double)*dy
         fy = sin(PI*y/Ly)**2
         dfdy = sin(2.0_c_double*PI*y/Ly)*(PI/Ly)

         u(jj, k_mode_seed) = u(jj, k_mode_seed) &
                              + cmplx(inv2Nx*eps_pert*dfdy, &
                                      0.0_c_double, c_double)

         v(jj, k_mode_seed) = v(jj, k_mode_seed) &
                              + cmplx(0.0_c_double, &
                                      -inv2Nx*eps_pert*kx_phys*fy, c_double)

         ! Hermitian conjugate (-k)
         u(jj, k_minus) = u(jj, k_minus) &
                          + cmplx(inv2Nx*eps_pert*dfdy, &
                                  0.0_c_double, c_double)

         v(jj, k_minus) = v(jj, k_minus) &
                          + cmplx(0.0_c_double, &
                                  +inv2Nx*eps_pert*kx_phys*fy, c_double)
      end do

      call Exchange_Halos(u, Ny, Nx, curr_rank_mpi, num_ranks_mpi, ierr_mpi)
      call Exchange_Halos(v, Ny, Nx, curr_rank_mpi, num_ranks_mpi, ierr_mpi)
   end subroutine Seed_Perturbation
end module turbulence
