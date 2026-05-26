!----------------------------------------------------------------------!
! Code to perform Heller dynamics using the ab initio trajectory &     !
! Hessian matrices in Dragonball format                                !
!----------------------------------------------------------------------!

      program heller

      use constants

      implicit none

      character(len=50) :: geom_eq,hess_eq,traj,hess

      namelist /input_files/ geom_eq,hess_eq,traj,hess

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

      read(111,input_files)

      write(*,*) "Equilibrium geometry from: ", geom_eq
      write(*,*) "Equilibrium Hessian from:  ", hess_eq
      write(*,*) "Trajectory from:           ", traj
      write(*,*) "Hessian matrices from:     ", hess


      end program heller
