!---------------------------------------------------------------------!
! MODULE CONTAINING ALL I NEED TO COMPUTE THE CORRELATION FUNCTION    !
!---------------------------------------------------------------------!

      module correlation_module

      use constants

      implicit none

      private
      public :: correlation,invgen,det_cmplx,logdet_cmplx

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

!.....Log determinant of a general complex matrix.........................
      function logdet_cmplx(nd,Amat) result(logdet)
      ! Returns log(det(A)) without explicitly forming det(A).
       integer, intent(in) :: nd
       complex*16, dimension(nd,nd), intent(in) :: Amat
       complex*16 :: logdet
       integer, dimension(nd) :: ipiv
       complex*16, dimension(nd,nd) :: Awork
       complex*16 :: diag
       integer :: i,n,info,nperm
       external ZGETRF

       n = nd
       Awork = Amat
       call ZGETRF(n,n,Awork,n,ipiv,info)
       if (info /= 0) then
         stop 'logdet_cmplx Matrix is numerically singular!'
       end if

       logdet = (0.d0,0.d0)
       nperm = 0
       do i = 1,nd
         diag = Awork(i,i)
         if (abs(diag) == 0.d0) then
           stop 'logdet_cmplx zero diagonal in LU factorization!'
         end if
         logdet = logdet + log(diag)
         if (ipiv(i) /= i) nperm = nperm + 1
       end do

       if (mod(nperm,2) /= 0) logdet = logdet + iu*pi
      end function

      subroutine correlation(nd,time,q0,p0,A0,qt,pt,At,ph0,pht,Ct)
      ! < 0 | t >
      ! Gaussian normalization is evaluated in logarithmic form.
       integer, intent(in) :: nd
       real*8, intent(in) :: time
       complex*16, intent(in) :: ph0,pht
       real*8, dimension(nd), intent(in) :: q0,p0,qt,pt
       real*8, dimension(nd,nd), intent(in) :: A0
       complex*16, dimension(nd,nd), intent(in) :: At

       real*8 :: p0q0,ptqt
       complex*16 :: q0A0q0,qtAtqt,bWb,Wlogdet,logcorr,c,Dph,Ct
       complex*16, dimension(nd) :: A0q0,Atqt,Wb,bvec
       complex*16, dimension(nd,nd) :: W,invW

       W = At + transpose(A0)
       invW = invgen(nd,W)
       Wlogdet = logdet_cmplx(nd,W)

       A0q0 = matmul(transpose(A0),q0)
       q0A0q0 = dot_product(q0,A0q0)
       Atqt = matmul(At,qt)
       qtAtqt = dot_product(qt,Atqt)

       p0q0 = dot_product(p0,q0)
       ptqt = dot_product(pt,qt)

       Dph = iu*(pht-conjg(ph0))
       c = iu*(p0q0-ptqt) - 0.5d0*qtAtqt - 0.5d0*q0A0q0 + Dph
       bvec = -iu*(p0-pt) + A0q0 + Atqt

       Wb = matmul(invW,bvec)
       bWb = dot_product(dconjg(bvec),Wb)

       ! log(Gint) = 1/2 [ nd*log(2*pi) - log(det(W)) ].
       ! Combine all exponential factors before evaluating exp().
       logcorr = 0.5d0*(nd*log(2.d0*pi) - Wlogdet) &
               + 0.5d0*bWb + c

       Ct = cdexp(logcorr)

       write(200,'(F12.5,3ES20.10)') time,real(Ct), &
                  aimag(Ct),dreal(Ct*dconjg(Ct))
      end subroutine

      end module
