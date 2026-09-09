!----------------------------------------------------------------------!
! Code to perform Heller dynamics using an ab initio trajectory &     !
! Hessian matrices in Cartesian coordinates                           !
!----------------------------------------------------------------------!

      program heller_cartesian

      use constants
      use correlation_module
      use, intrinsic :: iso_c_binding
      implicit none

      include 'fftw3.f03'

      real(C_DOUBLE), allocatable :: frequency(:),omega(:)
      complex(C_DOUBLE_COMPLEX), allocatable :: Coft(:), fft(:)
      type(C_PTR) :: plan
      integer :: i,j,k
      integer :: nat,ncart,steps,padding,calculation,fftsteps,rotrasl
      real*8 :: Epot0,L,dt,logdetA0,time,S,Eref,eta,Etot
      real*8, parameter :: width_floor_factor = 1.d-4
      real*8, parameter :: hessian_neg_tol = 1.d-7
      real*8 :: eigmax,eig_tol,eig_floor,eigval
      integer, allocatable, dimension(:) :: mask,mask_cart
      real*8, allocatable, dimension(:) :: xeq,veq,xm,x,v
      real*8, allocatable, dimension(:) :: p0,q0
      real*8, allocatable, dimension(:,:) :: hesseq,hessian,A0
      real*8, allocatable, dimension(:) :: ww
      real*8, allocatable, dimension(:,:) :: cnorm,Ad
      complex*16 :: trace,deltaph,ph0,pht,Ct
      complex*16, allocatable, dimension(:,:) :: At,Z,Y,invZ
      character(len=2), allocatable, dimension(:) :: symb
      character(len=50) :: geom_eq,hess_eq,traj,hess,velocity,output
      character(len=50) :: fourierout,powerout

      namelist /input_files/ geom_eq,hess_eq,traj,hess,velocity,output,&
                            &fourierout,powerout
      namelist /options/ calculation,eta,rotrasl
      namelist /trajectory/ steps,dt,padding,mask

      call execute_command_line('cat banner.txt')

      write(*,*) "+---------------------------------------------------+"
      write(*,*) "|            CARTESIAN MAIN CODE EXECUTION          |"
      write(*,*) "|              FWGP ONLY, CALCULATION IS            |"
      write(*,*) "|                      IGNORED                      |"
      write(*,*) "+---------------------------------------------------+"

!_____Input reading_____________________________________________________

      open(unit=111,file='input',status='old',action='read')

      geom_eq="equilibrium_geometry.xyz"
      hess_eq="Hessian_flat.out"
      traj="parsed_log_traj.xyz"
      hess="Hessian_traj.out"
      velocity="velocity.xyz"
      output="correlation.dat"
      fourierout="fourier.dat"
      powerout="power.dat"

      read(111,input_files)

      calculation=0
      eta = 0.d0
      rotrasl = 6       ! Retained for input compatibility; unused in Cartesian dynamics.

      read(111,options)

      write(*,*) "Equilibrium geometry from: ", geom_eq
      write(*,*) "Equilibrium Hessian from:  ", hess_eq
      write(*,*) "Trajectory from:           ", traj
      write(*,*) "Hessian matrices from:     ", hess
      write(*,*) "Initial velocities from:   ", velocity

!_____Masses and equilibrium geometry___________________________________

      open(unit=112,file=geom_eq,status="old",action="read")

      read(112,*) nat
      ncart = 3*nat

      allocate(symb(nat))
      allocate(xm(ncart),xeq(ncart),veq(ncart),hesseq(ncart,ncart))
      allocate(x(ncart),v(ncart),hessian(ncart,ncart),mask(nat),&
               &mask_cart(ncart))
      allocate(q0(ncart),p0(ncart),A0(ncart,ncart))
      allocate(At(ncart,ncart),Z(ncart,ncart),Y(ncart,ncart),&
               &invZ(ncart,ncart))
      allocate(ww(ncart),cnorm(ncart,ncart),Ad(ncart,ncart))

      read(112,*)
      do i = 1, nat
        read (112,*) symb(i), xeq(3*i-2:3*i)
        SELECT CASE (symb(i))
          CASE ('H')
            xm(3*i-2:3*i)=1837.15d0
          CASE ('D')
            xm(3*i-2:3*i)=3671.48d0
            symb(i) = 'H1'
          CASE ('C')
            xm(3*i-2:3*i)=21874.66d0
          CASE ('N')
            xm(3*i-2:3*i)=25526.06d0
          CASE ('O')
            xm(3*i-2:3*i)=29156.95d0
          CASE ('S')
            xm(3*i-2:3*i)=58422.43d0
          CASE ('P')
            xm(3*i-2:3*i)=7348.6d0
          CASE ('I')
            xm(3*i-2:3*i)=97368.95d0
          CASE default
            write(*,*) 'One or more atoms in your molecule is not',&
                      ' present in the database'
            stop
        END SELECT
      end do

      write(*,*) "@---------------------------------------------------@"
      write(*,*) "Equilibrium Geometry:"
      do i = 1,nat
         write(*,*) symb(i), xeq(3*i-2:3*i), xm(3*i-2)
      end do

!_____Equilibrium Cartesian Hessian and coherent width__________________

      open(unit=113,file=hess_eq,status="old",action="read")

      read(113,*)
      read(113,*)
      do i = 1,ncart
         do j = 1,i
            read(113,*) hesseq(j,i)
            hesseq(i,j) = hesseq(j,i)
         end do
      end do

!     Mass-weight the Cartesian Hessian.
      do i = 1,ncart
         do j = 1,ncart
            hesseq(i,j) = hesseq(i,j)/dsqrt(xm(i)*xm(j))
         end do
      end do

!     Diagonalize H in Cartesian space only to construct sqrt(H).
!     A0 = C diag(sqrt(lambda)) C^T.
      call diagonalizer(ncart,hesseq,ww,cnorm)

      eigmax = maxval(abs(ww))
      eig_tol = hessian_neg_tol!*eigmax
      eig_floor = (width_floor_factor**2)*eigmax

      if (eigmax <= 0.d0) then
        stop 'Equilibrium Hessian has no positive eigenvalues.'
      end if

      A0(:,:) = 0.d0
      Ad(:,:) = 0.d0
      logdetA0 = 0.d0

      do i = 1,ncart
        eigval = ww(i)

      !  Treat small negative eigenvalues as numerical noise.
      if (eigval < 0.d0 .and. abs(eigval) <= eig_tol) then
      write(*,*) 'Setting small negative Hessian eigenvalue to zero:', &
                 i, eigval
        eigval = 0.d0
      elseif (eigval < -eig_tol) then
      write(*,*) 'Warning: significant negative Hessian eigenvalue:', &
                 i, eigval
      stop 'Equilibrium Hessian is not positive semidefinite.'
      end if

     !  Zero modes still need a tiny finite width because the Gaussian
     !  matrices must remain invertible.
      if (eigval <= eig_floor) then
         eigval = eig_floor
      end if

      Ad(i,i) = dsqrt(eigval)
      logdetA0 = logdetA0 + log(Ad(i,i))
      end do

      A0 = matmul(cnorm,matmul(Ad,transpose(cnorm)))

      write(*,*) "@---------------------------------------------------@"
      write(*,*) "Cartesian Hessian / coherent Gaussian width"
      write(*,*) "Dimension: ", ncart
      write(*,*) "sqrt(H) diagonal eigenvalues:"
      do i = 1,ncart
         write(*,*) "Cartesian eigenvalue[",i,"]", Ad(i,i)
      end do
      write(*,*) "Zero-mode width floor: ", dsqrt(eig_floor)
      write(*,*) "@---------------------------------------------------@"

!_____Read trajectory and compute autocorrelation_______________________

      open(unit=114,file=traj,status="old",action="read")
      open(unit=115,file=hess,status="old",action="read")
      open(unit=116,file=velocity,status="old",action="read")
      open(unit=200,file=output,status="unknown",action="write")

      write(200,'(A12,3A20)') '# Time', 'Re[C(t)]', &
                              'Im[C(t)]', '|C(t)|**2'

!     Defaults retained for compatibility with the original code.
      steps = 2500
      dt = 8.2682749151502d0
      padding = 0
      mask(:) = 1
      mask_cart(:) = 1
      Epot0 = 0.d0
      read(111,trajectory)
      fftsteps = steps+padding

!     Expand the atom mask into Cartesian coordinates.
      do i = 1,nat
         mask_cart(3*i-2:3*i) = mask(i)
      end do

      allocate(Coft(fftsteps),fft(fftsteps),frequency(fftsteps),&
              &omega(fftsteps))

      write(*,*) "@---------------------------------------------------@"
      write(*,*) "Reading", steps, "steps of dynamics"
      write(*,*) "of",dt,"Dau each"
      write(*,*) "@---------------------------------------------------@"
      write(*,*) "Final damping: ", exp(-eta*steps*dt)

!_____Initial conditions________________________________________________

      read(116,*)
      read(116,*)
      do i = 1,nat
        read(116,*) symb(i), veq(3*i-2:3*i)
      end do

!     Convert positions to Bohr; velocities are already in atomic units.
      x(:) = xeq(:)/bohr_radius
      v(:) = veq(:)

!     Mass-weight Cartesian coordinates and momenta.
      x(:) = x(:)*dsqrt(xm(:))
      v(:) = v(:)*dsqrt(xm(:))

      q0 = x
      p0 = mask_cart(:)*v(:)
      !p0 = v(:)

      Y = iu*A0
      Z = (0.d0,0.d0)
      do i = 1,ncart
         Z(i,i) = (1.d0,0.d0)
      end do

      Coft(:) = cmplx(0.d0,0.d0,kind=C_DOUBLE)
      At = cmplx(A0)

!     A0 is positive definite after the zero-mode regularization.
      ph0 = -iu*0.25d0*(logdetA0 - ncart*log(pi))
      pht = ph0
      S = 0.d0
      time = 0.d0

!     Same reference-energy convention as the original implementation,
!     now evaluated in the full Cartesian space.
      trace = (0.d0,0.d0)
      do i = 1,ncart
         trace = trace + A0(i,i)
      end do
      Eref = 0.5d0*real(trace)
      Etot = dot_product(p0,p0)/2.d0 + Epot0

      write(*,*) "Coherent Cartesian ZPE-like reference: ",&
                 & Eref, Eref*Ha2cmm1

      call correlation(ncart,time,q0,p0,A0,x,p0,At,ph0,ph0,Ct)
      Coft(1) = cmplx(Ct,kind=C_DOUBLE)

!_____Propagation______________________________________________________

      do k = 1,steps
        time = time + dt

        read(114,*)
        read(114,*)
        do i = 1,nat
          read(114,*) symb(i), x(3*i-2:3*i), v(3*i-2:3*i)
        end do

!       Convert to AU.
        x(:) = x(:)/bohr_radius
        v(:) = v(:)

!       Mass-weight Cartesian coordinates and momenta.
        x(:) = x(:)*dsqrt(xm(:))
        v(:) = v(:)*dsqrt(xm(:))

!       Atom mask: active atoms follow the trajectory; inactive atoms
!       remain at their t=0 Cartesian center and momentum.
        x(:) = x(:)*mask_cart(:) + q0(:)*(1-mask_cart(:))
        v(:) = v(:)*mask_cart(:) + p0(:)*(1-mask_cart(:))

!       Evolves A.
!       FROZEN GAUSSIAN
        At = cmplx(A0)

!       Second half of the action.
        L = dot_product(v,v) - Etot
        trace = (0.d0,0.d0)
        do i = 1,ncart
          trace = trace + At(i,i)
        end do
        pht = pht + (L - 0.5d0*trace)*dt

        call correlation(ncart,time,q0,p0,A0,x,v,At,ph0,pht,Ct)
        Coft(k) = cmplx(Ct*exp(-eta*time+iu*0.75d0*Eref*time)&
                  &,kind=C_DOUBLE)
      end do

      write(*,*) "@---------------------------------------------------@"
      write(*,*) "Final Snapshot (mass-weighted!):"
      do i = 1,nat
         write(*,*) symb(i), x(3*i-2:3*i), v(3*i-2:3*i)
      end do
      write(*,*) "@---------------------------------------------------@"
      if (calculation.ge.2) then
         write(*,*) "Final Hessian (1,1) entry:"
         write(*,*) hessian(1,1)
      write(*,*) "@---------------------------------------------------@"
      end if

!_____Fourier___________________________________________________________

      do i = 1, fftsteps
         if (i <= fftsteps/2 + 1) then
           frequency(i)=dble(i-1)/(dble(fftsteps)*dt)
         else
           frequency(i)=dble(i-1-fftsteps)/(dble(fftsteps)*dt)
         end if
      end do

      omega(:) = 2*pi*frequency(:) + 0.75d0*Eref

      plan=fftw_plan_dft_1d(fftsteps,Coft,fft,&
                            &FFTW_BACKWARD,FFTW_ESTIMATE)
      call fftw_execute_dft(plan,Coft,fft)
      call fftw_destroy_plan(plan)

      fft(:) = dt*fft(:)

      open(unit=222,file=fourierout,status='replace',action='write')

      write(222,'(A12,3A20)') '# Freq','Ang Freq','Re[FFT]','Im[FFT]'

      do i = 1, fftsteps
        write(222,'(F12.5,3ES20.10)') frequency(i), omega(i),&
                     &dreal(fft(i)), aimag(fft(i))
      end do

      close(222)

      open(unit=333,file=powerout,status='replace',action='write')

      write(333,'(A12,3A20)') '# Freq','Ang Freq [au]',&
                              &'Ang Freq [cm**-1]','|FFT|**2'

      do i = 1, fftsteps/2+1
        write(333,'(F12.5,3ES20.10)') frequency(i), omega(i), &
          &omega(i)*219474.6313705,(dreal(fft(i))**2 +aimag(fft(i))**2)
      end do

      close(333)

!_____Closing and deallocating__________________________________________

      deallocate(xm,xeq,veq,hesseq,symb,ww,cnorm,Ad,x,v,hessian,q0,p0,&
                 A0,At,Z,Y,invZ,frequency,omega,Coft,fft,mask,mask_cart)

      close(111)
      close(112)
      close(113)
      close(114)
      close(115)
      close(116)
      close(200)

      end program heller_cartesian
