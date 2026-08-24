!----------------------------------------------------------------------!
! Code to perform Heller dynamics using the ab initio trajectory &     !
! Hessian matrices in Dragonball format                                !
!----------------------------------------------------------------------!

      program heller

      use constants
      use correlation_module
      use, intrinsic :: iso_c_binding
      implicit none

      include 'fftw3.f03'

      real(C_DOUBLE), allocatable :: frequency(:),omega(:)
      complex(C_DOUBLE_COMPLEX), allocatable :: Coft(:), fft(:)
      type(C_PTR) :: plan
      integer :: i,j,k,rotrasl
      integer :: nat,ncart,nvib,steps,padding,calculation,fftsteps
      real*8 :: Epot0,L,dt,detA0,time,S,Eref,eta,Etot
      integer, allocatable, dimension(:) :: mask
      real*8, allocatable, dimension(:) :: xeq,veq,xm,ww,x,v
      real*8, allocatable, dimension(:) :: qrt,prt,qvib,pvib 
      real*8, allocatable, dimension(:) :: p0,q0 
      real*8, allocatable, dimension(:,:) :: hesseq,cnorm,tmp,hessian
      real*8, allocatable, dimension(:,:) :: Hrt,Hvib,A0,HA,invA0
      complex*16 :: dotph,trace,deltaph,ph0,pht,detZi,traceHA,Ct
      complex*16, allocatable, dimension(:,:) :: At,Z,Y,invZ,dotAt
      character(len=2), allocatable, dimension(:) :: symb
      character(len=50) :: geom_eq,hess_eq,traj,hess,velocity,output
      character(len=50) :: fourierout,powerout 

      namelist /input_files/ geom_eq,hess_eq,traj,hess,velocity,output,&
                            &fourierout,powerout
      namelist /options/ calculation,eta,rotrasl
      namelist /trajectory/ steps,dt,padding,mask

      call execute_command_line('cat banner.txt')

      write(*,*) "+---------------------------------------------------+"
      write(*,*) "|               MAIN CODE EXECUTION                 |"
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
      rotrasl = 6
 
      read(111,options)

      write(*,*) "Equilibrium geometry from: ", geom_eq
      write(*,*) "Equilibrium Hessian from:  ", hess_eq
      write(*,*) "Trajectory from:           ", traj
      write(*,*) "Hessian matrices from:     ", hess
      write(*,*) "Initial velocities from:   ", velocity

      if(calculation.eq.0) then
         write(*,*) "Calculation:     ", "Frozen"
      elseif(calculation.eq.1) then
         write(*,*) "Calculation:     ", "Single Hessian"
      else
         write(*,*) "Calculation:     ", "Thawed"
      end if


!_____Masses and equilibrium geometry___________________________________

      open(unit=112,file=geom_eq,status="old",action="read")

      read(112,*) nat

      ncart = 3*nat
      nvib = ncart - rotrasl 

      allocate(symb(nat))                                           
      allocate(xm(ncart),xeq(ncart),veq(ncart),hesseq(ncart,ncart))    
      allocate(x(ncart),v(ncart),hessian(ncart,ncart),mask(nvib))    
      allocate(qrt(ncart),prt(ncart),Hrt(ncart,ncart),dotAt(nvib,nvib))
      allocate(qvib(nvib),pvib(nvib),Hvib(nvib,nvib),invA0(nvib,nvib))
      allocate(q0(nvib),p0(nvib),A0(nvib,nvib),At(nvib,nvib))
      allocate(Z(nvib,nvib),Y(nvib,nvib),invZ(nvib,nvib),HA(nvib,nvib))
      allocate(ww(ncart),cnorm(ncart,ncart),tmp(ncart,ncart))    

      read(112,*)
      do i = 1, nat
        read (112,*) symb(i), xeq(3*i-2:3*i)
        SELECT CASE (symb(i))
          CASE ('H')
            xm(3*i-2:3*i)=1837.15d0
          CASE ('D')  !Deuterium
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
            xm(3*i-2:3*i)=97368.95
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

!_____Equilibrium Hessian and Normal modes______________________________

      open(unit=113,file=hess_eq,status="old",action="read")

      read(113,*)
      read(113,*)
      do i = 1,ncart
         do j = 1,i
            read(113,*) hesseq(j,i)
            hesseq(i,j) = hesseq(j,i)
         end do
      end do

      do i = 1,ncart
         do j = 1,ncart
            hesseq(i,j) = hesseq(i,j)/dsqrt(xm(i)*xm(j))
         end do
      end do

      call diagonalizer(ncart,hesseq,ww,cnorm)

      do j = 1, nvib                    ! vibrational modes index
         do i = 1, ncart                ! eigenvector index
            tmp(i,j) = cnorm(i,j+ncart-nvib)
         enddo
      enddo

      do j = nvib+1, ncart
         do i = 1, ncart
            tmp(i,j) = cnorm(i,j-nvib)
         enddo
      enddo
      do i = 1, ncart
         do j = 1, ncart
            cnorm(i,j) = tmp(i,j)
         enddo
      enddo
      
      tmp = matmul(hesseq,cnorm)
      Hrt = matmul(transpose(cnorm),tmp)
      write(*,*) "@---------------------------------------------------@"
      write(*,*) "Harmonic frequencies"
      do i = 1,nvib
         write(*,*) "Mode[",i,"]", dsqrt(ww(rotrasl+i))*Ha2cmm1
      end do
      write(*,*) "Harmonic zpe: ",&
                 & sum(dsqrt(ww(rotrasl+1:ncart)))*Ha2cmm1/2.d0
      write(*,*) "@---------------------------------------------------@"
      write(*,*) "Harmonic frequencies from NM Hessian"
      do i = 1,ncart
         write(*,*) "Mode[",i,"]", dsqrt(abs(Hrt(i,i)))*Ha2cmm1
      end do
      write(*,*) "This value should be zero:", Hrt(2,3)
 


!_____Read trajectory and compute autocorrelation_______________________
 
      open(unit=114,file=traj,status="old",action="read")
      open(unit=115,file=hess,status="old",action="read")
      open(unit=116,file=velocity,status="old",action="read")
      open(unit=200,file=output,status="unknown",action="write")

      write(200,'(A12,3A20)') '# Time', 'Re[C(t)]', &
                              'Im[C(t)]', '|C(t)|**2'

      steps = 2500
      dt = 8.2682749151502d0 
      padding = 0
      mask(:) = 1
      Epot0 = 0.d0
      read(111,trajectory)
      fftsteps = steps+padding

      allocate(Coft(fftsteps),fft(fftsteps),frequency(fftsteps)&
               &,omega(fftsteps))

      write(*,*) "@---------------------------------------------------@"
      write(*,*) "Reading", steps, "steps of dynamics"
      write(*,*) "of",dt,"Dau each"
      write(*,*) "@---------------------------------------------------@"
      write(*,*) "Final damping: ", exp(-eta*steps*dt)

      ! INITIAL CONDITIONS
      ! Read initial velocity
      read(116,*)
      read(116,*)
      do i = 1,nat
        read(116,*) symb(i), veq(3*i-2:3*i)
      end do
      ! Convert to AU
      x(:) = xeq(:)/bohr_radius
      v(:) = veq(:)!*FROMangTOau_vel Velocities are already in atomic units
      ! Mass scale
      x(:) = x(:)*dsqrt(xm(:))
      v(:) = v(:)*dsqrt(xm(:))
      ! Normal modes
      qrt = matmul(transpose(cnorm),x)
      prt = matmul(transpose(cnorm),v)
      ! Vibrational only
      qvib(:) = qrt(1:nvib)
      pvib(:) = prt(1:nvib)
      Hvib(:,:) = Hrt(1:nvib,1:nvib)

      q0 = qvib !-qvib
      p0 = pvib !-pvib

      A0(:,:) = 0.d0
      invA0(:,:) = 0.d0
      detA0 = 1.d0
      do i = 1,nvib
        A0(i,i) = dsqrt(ww(rotrasl+i))
        Z(i,i) = (1.d0,0.d0)!*mask(i)
        detA0 = detA0*A0(i,i)
        invA0(i,i) = 1/dsqrt(ww(rotrasl+i))
      end do 
      Y(:,:) = iu*A0(:,:)!*Z(:,:)
 
      ! COMPUTE C(0)
      Coft(:) = cmplx(0.d0,0.d0,kind=C_DOUBLE)
      At = cmplx(A0)
      !write(502,*) time, real(At(1,1)), aimag(At(1,1))
      !ph0 = 0.d0!-iu*0.25d0*log(detA0/pi**nvib)
      ph0 = -iu*0.25d0*log(detA0/pi**nvib)
      pht = ph0
      S = 0.d0
      time = 0.d0
      ! REFERENCE ENERGY 
      trace = (0.d0,0.d0)
      HA = matmul(Hvib,invA0)
      traceHA = (0.d0,0.d0)
      do i = 1,nvib
        trace = trace + A0(i,i)
        traceHA = traceHA + HA(i,i)
      end do

      !Eref= +0.25d0*trace + dot_product(p0,p0)/2.d0 +0.25*traceHA
      Eref = 0.5d0*trace
      Etot = dot_product(p0,p0)/2.d0 + Epot0 ! Epot of reference
 
      write(*,*) "Harmonic ZPE: ", Eref, Eref*Ha2cmm1
      
      call correlation(nvib,time,q0,p0,A0,qvib,pvib,At,ph0,ph0,detA0,Ct)
      Coft(1) = cmplx(Ct,kind=C_DOUBLE)
  
      do k = 1,steps
        time = time + dt
        read(114,*) 
        read(114,*) 
        do i = 1,nat
        read(114,*) symb(i), x(3*i-2:3*i), v(3*i-2:3*i)
        end do
        ! Convert to AU
        x(:) = x(:)/bohr_radius
        v(:) = v(:)!/toautime!*FROMangTOau_vel Velocities are already in atomic units
        ! Mass scale
        x(:) = x(:)*dsqrt(xm(:))
        v(:) = v(:)*dsqrt(xm(:))
        ! Normal modes
        qrt = matmul(transpose(cnorm),x)
        prt = matmul(transpose(cnorm),v)
        ! Vibrational only
        qvib(:) = qrt(1:nvib)*mask(:) + (1-mask(:))*q0(:) 
        pvib(:) = prt(1:nvib)*mask(:) + (1-mask(:))*p0(:)
        ! Evolves A
        ! FROZEN GAUSSIAN
        if(calculation.eq.0) then
           At = cmplx(A0)
        ! SINGLE HESSIAN
        elseif(calculation.eq.1) then
           !dotAt = +iu*matmul(At,At)+iu*Hvib
           !At = At + dotAt*dt
           Y = Y - matmul(Hvib,Z)*dt
           Z = Z + Y*dt
           invZ = invgen(nvib,Z)
           At = - iu*matmul(Y,invZ)
        ! THAWED GAUSSIAN
        else
           read(115,*)
           read(115,*)
           do i = 1,ncart
              do j = 1,i
                 read(115,*) hessian(j,i)
                 hessian(i,j) = hessian(j,i)
              end do
           end do
           do i = 1,ncart
              do j = 1,ncart
                 hessian(i,j) = hessian(i,j)/dsqrt(xm(i)*xm(j))
              end do
           end do
           tmp = matmul(hessian,cnorm)
           Hrt = matmul(transpose(cnorm),tmp)
           Hvib(:,:) = Hrt(1:nvib,1:nvib)
           Y = Y - matmul(Hvib,Z)*dt
           Z = Z + Y*dt
           invZ = invgen(nvib,Z)
           At = - iu*matmul(Y,invZ)
        end if
        !write(502,*) time, real(At(1,1)), aimag(At(1,1))
        ! Second half of the action 
        !L = dot_product(pvib,pvib)/2.d0 - Epot 
        L = dot_product(pvib,pvib) - Etot
        trace = (0.d0,0.d0)
        do i = 1,nvib
          trace = trace + At(i,i)
        end do
        pht = pht + (L -0.5d0*trace)*dt 
        !write(*,*) L
        !write(505,*) time, real(pht), aimag(pht)
      call correlation(nvib,time,q0,p0,A0,qvib,pvib,At,ph0,pht,detA0,Ct)
        Coft(k) = cmplx(Ct*exp(-eta*time) ,kind=C_DOUBLE)
      end do !k

      write(*,*) "@---------------------------------------------------@"
      write(*,*) "Final Snapshot (mass-weighted!):"
      do i = 1,nat
         write(*,*) symb(i), x(3*i-2:3*i), v(3*i-2:3*i)
      end do
      write(*,*) "@---------------------------------------------------@"
      write(*,*) "Final Hessian (1,1) entry:"
      write(*,*) hessian(1,1)
      write(*,*) "@---------------------------------------------------@"
      !write(*,*) "Final Energy entry:"
      !write(*,*) Epot
      
      ! PRINT PADDING TO HAVE BIGGER OUTPUTSSSSSS
 
      !do k = 1,padding
      !  time = time + dt
      !  write(200,*) time, 0.d0, 0.d0, 0.d0
      !end do

!_____Fourier___________________________________________________________
 
      do i = 1, fftsteps
         if (i <= fftsteps/2 + 1) then
           frequency(i)=dble(i-1)/(dble(fftsteps)*dt)
         else
           frequency(i)=dble(i-1-fftsteps)/(dble(fftsteps)*dt)
         end if
      end do

      omega(:) = 2*pi*frequency(:)

      plan=fftw_plan_dft_1d(fftsteps,Coft,fft,&
                            &FFTW_BACKWARD,FFTW_ESTIMATE)
      call fftw_execute_dft(plan,Coft,fft)
      call fftw_destroy_plan(plan)

      fft(:) = dt*fft(:)!/dble(n)

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

      deallocate(xm,xeq,veq,hesseq,symb,ww,cnorm,tmp,x,v,hessian,q0,p0&
                 ,A0,At,Z,Y,invZ,HA,invA0,frequency,omega,Coft,fft&
                 ,dotAt)

      close(111)
      close(112)
      close(113)
      close(114)
      close(115)
      close(116)
      close(200)
      close(222)
      close(333)

      end program heller
