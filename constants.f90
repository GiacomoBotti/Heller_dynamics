!moule containing teh constant definitions

MODULE constants

     implicit none
     save

     real*8, parameter :: bohr_radius = 0.52917721067 
     real*8, parameter :: Ha2cmm1 = 219474.625d0
     real*8, parameter :: Ha2K = 315775.1293573255
     real*8, parameter :: NWchem_toauvel = 1.d0/21.87676d0
     real*8, parameter :: pi = 4.d0*datan(1.d0)
     real*8, parameter :: Ha2eV = 27.2114d0
     real*8, parameter :: toautime = 4.1341d4
     real*8, parameter :: toauvel = 1.d0!/21.87676d0

     complex*16,parameter :: iu = complex(0.d0,1.d0)
     
     logical, parameter :: cnorm_tascivr = .true.

END MODULE constants
