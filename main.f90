!----------------------------------------------------------------------!
! Code to perform Heller dynamics using the ab initio trajectory &     !
! Hessian matrices in Dragonball format                                !
!----------------------------------------------------------------------!

      program heller

      use constants

      implicit none

      integer :: i,j,k
      integer :: nat,ncart,nvib,steps
      real*8 :: Epot,L,dt
      real*8, allocatable, dimension(:) :: xeq,veq,xm,ww,x,v
      real*8, allocatable, dimension(:) :: qrt,prt,qvib,pvib 
      real*8, allocatable, dimension(:,:) :: hesseq,cnorm,tmp,hessian
      real*8, allocatable, dimension(:,:) :: Hrt,Hvib 
      character(len=2), allocatable, dimension(:) :: symb
      character(len=50) :: geom_eq,hess_eq,traj,hess,energy,output

      namelist /input_files/ geom_eq,hess_eq,traj,hess,energy,output
      namelist /trajectory/ steps,dt

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
      energy="energies.dat"
      output="correlation.dat"

      read(111,input_files)

      write(*,*) "Equilibrium geometry from: ", geom_eq
      write(*,*) "Equilibrium Hessian from:  ", hess_eq
      write(*,*) "Trajectory from:           ", traj
      write(*,*) "Hessian matrices from:     ", hess
      write(*,*) "Potential energy from:     ", energy

!_____Masses and equilibrium geometry___________________________________

      open(unit=112,file=geom_eq,status="old",action="read")

      read(112,*) nat

      ncart = 3*nat
      nvib = ncart - 6

      allocate(symb(nat))                                           
      allocate(xm(ncart),xeq(ncart),veq(ncart),hesseq(ncart,ncart))    
      allocate(x(ncart),v(ncart),hessian(ncart,ncart))    
      allocate(qrt(ncart),prt(ncart),Hrt(ncart,ncart))    
      allocate(qvib(nvib),pvib(nvib),Hvib(nvib,nvib))
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
          CASE ('SP')
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
         write(*,*) "Mode[",i,"]", dsqrt(ww(6+i))*Ha2cmm1
      end do
      write(*,*) "@---------------------------------------------------@"
      write(*,*) "Harmonic frequencies from NM Hessian"
      do i = 1,ncart
         write(*,*) "Mode[",i,"]", dsqrt(abs(Hrt(i,i)))*Ha2cmm1
      end do
      write(*,*) "This value should be zero:", Hrt(2,3)
 


!_____Read trajectory and compute autocorrelation_______________________
 
      open(unit=114,file=traj,status="old",action="read")
      open(unit=115,file=hess,status="old",action="read")
      open(unit=116,file=energy,status="old",action="read")

      steps = 2500
      dt = 8.2682749151502d0 
      read(111,trajectory)

      write(*,*) "@---------------------------------------------------@"
      write(*,*) "Reading", steps, "steps of dynamics"
      write(*,*) "of",dt,"Dau each"

      ! INITIAL CONDITIONS
      read(113,*) 
      read(114,*) 
      do i = 1,nat
      read(114,*) symb(i), x(3*i-2:3*i), v(3*i-2:3*i)
      end do
      read(115,*)
      read(115,*)
      do i = 1,ncart
         do j = 1,i
            read(115,*) hessian(j,i)
             hessian(i,j) = hessian(j,i)
         end do
      end do
      read(116,*) Epot 
      ! Convert to AU
      x(:) = x(:)/bohr_radius
      v(:) = v(:)!*FROMangTOau_vel Velocities are already in atomic units
      ! Mass scale
      x(:) = x(:)*dsqrt(xm(:))
      v(:) = v(:)*dsqrt(xm(:))
      ! Normal modes
      qrt = matmul(transpose(cnorm),x)
      prt = matmul(transpose(cnorm),v)
      tmp = matmul(hessian,cnorm)
      Hrt = matmul(transpose(cnorm),tmp)
      ! Vibrational only
      qvib(:) = qrt(1:nvib)
      pvib(:) = prt(1:nvib)
      Hvib(:,:) = Hrt(1:nvib,1:nvib)

      q0 = qvib
      p0 = pvib

      do i = 1,nvib
        A0 = dsqrt(ww(6+i))
      end do 

      ! COMPUTE C(0)

      do k = 2,steps
        read(114,*) 
        read(114,*) 
        do i = 1,nat
        read(114,*) symb(i), x(3*i-2:3*i), v(3*i-2:3*i)
        end do
        read(115,*)
        read(115,*)
        do i = 1,ncart
           do j = 1,i
              read(115,*) hessian(j,i)
               hessian(i,j) = hessian(j,i)
           end do
        end do
        read(116,*) Epot 
        ! Convert to AU
        x(:) = x(:)/bohr_radius
        v(:) = v(:)!*FROMangTOau_vel Velocities are already in atomic units
        ! Mass scale
        x(:) = x(:)*dsqrt(xm(:))
        v(:) = v(:)*dsqrt(xm(:))
        ! Normal modes
        qrt = matmul(transpose(cnorm),x)
        prt = matmul(transpose(cnorm),v)
        tmp = matmul(hessian,cnorm)
        Hrt = matmul(transpose(cnorm),tmp)
        ! Vibrational only
        qvib(:) = qrt(1:nvib)
        pvib(:) = prt(1:nvib)
        Hvib(:,:) = Hrt(1:nvib,1:nvib)
        ! Evolves Delta gamma
        L = dot_product(pvib,pvib) - V
        ! Evolves width



      end do !k

      write(*,*) "@---------------------------------------------------@"
      write(*,*) "Final Snapshot:"
      do i = 1,nat
         write(*,*) symb(i), x(3*i-2:3*i), v(3*i-2:3*i)
      end do
      write(*,*) "@---------------------------------------------------@"
      write(*,*) "Final Hessian (1,1) entry:"
      write(*,*) hessian(1,1)
      write(*,*) "@---------------------------------------------------@"
      write(*,*) "Final Energy entry:"
      write(*,*) Epot
      


!_____Closing and deallocating__________________________________________

      deallocate(xm,xeq,veq,hesseq,symb,ww,cnorm,tmp,x,v,hessian)

      close(111)
      close(112)
      close(113)
      close(114)
      close(115)
      close(116)

      end program heller
