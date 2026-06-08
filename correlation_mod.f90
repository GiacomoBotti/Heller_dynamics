!---------------------------------------------------------------------!
! MODULE CONTAINING ALL I NEED TO COMPUTE THE CORRELATION FUNCTION    !
!---------------------------------------------------------------------!

      module correlation_module

      use constants

      implicit none

      private
      public :: correlation,invgen,det_cmplx

      contains 

!.....Inversion for a general COMPLEX matrix............................
! Returns the inverse of a matrix calculated by finding the LU
! decomposition.  Depends on LAPACK.
      function invgen(npar,A) result(Ainv)
        implicit none
        integer, intent(in) :: npar
        complex*16, dimension(npar,npar), intent(in) :: A
        complex*16, dimension(npar,npar) :: Ainv

        complex*16, dimension(npar) :: work  ! work array for LAPACK
        integer, dimension(npar) :: ipiv   ! pivot indices
        integer :: n, info

        ! External procedures defined in LAPACK
        external ZGETRF
        external ZGETRI

        ! Store A in Ainv to prevent it from being overwritten by LAPACK
        Ainv = A
        n = npar

        ! DGETRF computes an LU factorization of a general M-by-N matrix A
        ! using partial pivoting with row interchanges.
        call ZGETRF(n, n, Ainv, n, ipiv, info)

        if (info /= 0) then
          !write(*,*) "DGETRF info : ",info
          stop 'Invgen Matrix is numerically singular!'
        end if

        ! DGETRI computes the inverse of a matrix using the LU factorization
        ! computed by DGETRF.
        call ZGETRI(n, Ainv, n, ipiv, work, n, info)

        if (info /= 0) then
          !write(*,*) "DGETRI info : ",info
          stop 'Invgen Matrix inversion failed!'
        end if
      end function 

!.....Diagonalization with LAPACK (complex).............................

      function det_real(nd,Amat) result(detAmat) 
      ! nd: dimension of the matrices
      ! Amat: real matrix matrix (precision)
      ! detAmat: determinant

       integer, intent(in) :: nd
       real*8, dimension(nd,nd), intent(in) :: Amat
 
       real*8 :: detAmat 


       integer, dimension(nd) :: ipiv   ! pivot indices
       real*8 :: detL,detP
       real*8 :: detU
       real*8, dimension(nd) :: work  ! work array for LAPACK
       real*8, dimension(nd,nd) :: Awork
       integer :: i,n, info

       external DGETRF

       n = nd
       ! Store A in Ainv to prevent it from being overwritten by LAPACK
       Awork = Amat

       call DGETRF(n,n,Awork,n,ipiv,info) 

       if (info /= 0) then
         !write(*,*) "DGETRF info : ",info
         stop 'diag_complx Matrix is numerically singular!'
       end if

       ! Determinants of the decomposition 
       detU = 1.d0
       detL = 1.d0
       detP = 1.d0

       do i = 1,nd
         detU = detU*Awork(i,i)
         if(ipiv(i).ne.i) then
            detP = - detP
         end if
       end do
 
       ! Total determinant

       detAmat = detP*detU*detL

      end function

!.....Diagonalization with LAPACK (complex).............................

      function det_cmplx(nd,Amat) result(detAmat) 
      ! nd: dimension of the matrices
      ! Amat: complex matrix matrix (precision)
      ! detAmat: determinant

       integer, intent(in) :: nd
       complex*16, dimension(nd,nd), intent(in) :: Amat
 
       complex*16 :: detAmat 


       integer, dimension(nd) :: ipiv   ! pivot indices
       real*8 :: detL,detP
       complex*16 :: detU
       complex*16, dimension(nd) :: work  ! work array for LAPACK
       complex*16, dimension(nd,nd) :: Awork
       integer :: i,n, info

       external ZGETRF

       n = nd
       ! Store A in Ainv to prevent it from being overwritten by LAPACK
       Awork = Amat

       call ZGETRF(n,n,Awork,n,ipiv,info) 

       if (info /= 0) then
         !write(*,*) "DGETRF info : ",info
         stop 'diag_complx Matrix is numerically singular!'
       end if

       ! Determinants of the decomposition 
       detU = complex(1.d0,0.d0)
       detL = 1.d0
       detP = 1.d0

       do i = 1,nd
         detU = detU*Awork(i,i)
         if(ipiv(i).ne.i) then
            detP = - detP
         end if
       end do
 
       ! Total determinant

       detAmat = detP*detU*detL

      end function

!_____Correlation function______________________________________________

      subroutine correlation(nd,time,q0,p0,A0,qt,pt,At,ph0,pht,detA0)
      ! < 0 | t >
      ! nd: system dimensions
      ! time: simulation time
      ! q0, p0, A0: initial gaussian center, momentum and width
      ! qt, pt, At: instantaneous gaussian center, momentum and width
      ! deltaph: phase difference
      ! detA0 :: determinant of A0
       integer, intent(in) :: nd
       real*8, intent(in) :: time,detA0
       complex*16, intent(in) :: ph0,pht
       real*8, dimension(nd), intent(in) :: q0,p0,qt,pt
       real*8, dimension(nd,nd) :: A0
       complex*16, dimension(nd,nd), intent(in) :: At

       integer :: i
       real*8 :: p0q0,ptqt,N0,Nt,Adet,eta
       complex*16 :: q0A0q0,qtAtqt,bWb,Wdet,corr,c,gint,Dph
       complex*16, dimension(nd) :: A0q0,Atqt,Wb,bvec
       complex*16, dimension(nd,nd) :: W,invW

        eta = 0.d0
        N0 = 1.d0! (detA0/pi**nd)**(1.d0/4.d0)
        Adet = det_cmplx(nd,At)
        Nt = 1.d0!(Adet/pi**nd)**(1.d0/4.d0)

        W = (At + transpose(A0))
        invW = invgen(nd,W)
        Wdet = det_cmplx(nd,W)

        Gint = zsqrt((2.d0*pi)**nd/Wdet)

        A0q0 = matmul(transpose(A0),q0)
        q0A0q0 = dot_product(q0,A0q0)
        Atqt = matmul(At,qt)
        qtAtqt = dot_product(qt,Atqt)

        p0q0 = dot_product(p0,q0)
        ptqt = dot_product(pt,qt)

        Dph = iu*(pht-conjg(ph0)) 
        c = iu*(p0q0-ptqt) -0.5d0*qtAtqt-0.5d0*q0A0q0 +Dph
        bvec = -iu*(p0-pt) + A0q0 + Atqt

        Wb = matmul(invW,bvec)
        bWb = dot_product(dconjg(bvec),Wb)

        !write(*,*) bWb

        corr = Gint*Nt*N0*cdexp(0.5d0*bWb + c - eta*time)
        !corr = Gint*cdexp(0.5d0*bWb + c - eta*time)

       write(200,*) time,real(corr),aimag(corr),dreal(corr*dconjg(corr))
       write(201,*) time, pt(1) ,0.d0
       write(202,*) time, pt(2) ,0.d0
       write(203,*) time, pt(3) ,0.d0
       write(204,*) time, real(exp(Dph)),aimag(exp(Dph)),&
                    dreal(exp(Dph)*conjg(exp(Dph)))
       
      end subroutine

      end module
