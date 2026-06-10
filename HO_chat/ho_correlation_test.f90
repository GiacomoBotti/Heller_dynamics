      program ho_correlation_test
      implicit none

      integer, parameter :: dp = selected_real_kind(15, 307)
      integer, parameter :: nsteps =5000

      real(dp), parameter :: m    = 1.0_dp
      real(dp), parameter :: hbar = 1.0_dp
      real(dp), parameter :: omega = 20.0_dp
      real(dp), parameter :: dt    = 0.01_dp

      real(dp) :: q0, p0, t,damp
      complex(dp) :: alpha, C, iunit
      integer :: i, unit

      ! Initial coherent state
      ! For a simple ground-state test, use q0 = 0, p0 = 0.
      q0 = 0.0_dp
      p0 = 0.0_dp

      iunit = cmplx(0.0_dp, 1.0_dp, kind=dp)

      alpha = sqrt(m*omega/(2.0_dp*hbar))*q0 + &
              iunit * p0 / sqrt(2.0_dp*m*hbar*omega)

      open(newunit=unit,file='corr_ho.dat',status='replace',&
           &action='write')
      write(unit,'(A)') '# time   ReC(t)   ImC(t)   |C(t)|^2'

      do i = 0, nsteps
         t = i * dt
         damp = exp(-0.46d0*t/nsteps/dt)
         C = exp( +abs(alpha)**2 * (1.0_dp - exp(-iunit*omega*t)) + &
              0.5_dp*iunit*omega*t )*damp
         !C = exp(+0.5_dp*iunit*omega*t)

     write(unit,'(F12.6,1X,ES24.16,1X,ES24.16,1X,ES24.16)') &
          t, real(C,dp), aimag(C), abs(C)**2
      end do

      close(unit)
      end program ho_correlation_test
