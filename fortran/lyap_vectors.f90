! lyap_vectors.f90
!
!> Module for computation of Lyapunov exponents and vectors
!
!> @copyright                                                               
!> 2016 Sebastian Schubert.
!
!---------------------------------------------------------------------------!
!                                                                           
!>  @remark                                                                 
!>  This module contains the necessary tools to perform the Bennettin
!>  steps to compute the lyapunov exponents. (Ginelli for CLV will be added
!>  later) 
!                                                                           
!---------------------------------------------------------------------------


MODULE lyap_vectors

  !-----------------------------------------------------!
  !                                                     !
  ! Preamble and variables declaration                  !
  !                                                     !
  !-----------------------------------------------------!

  USE params, only: ndim,dt,tw
  IMPLICIT NONE
  
  PUBLIC :: bennettin_step,loclyap,lyapunov,ensemble,prop,init_lyap,multiply_prop,init_one
 
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: loclyap
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: lyapunov
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: ensemble !> ensemble contains the QR decompsoition of the ensemble
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: prop
  
  PRIVATE
  INTEGER :: lwork
  REAL(kind=8), DIMENSION(:), ALLOCATABLE :: work
  REAL(kind=8), DIMENSION(:), ALLOCATABLE :: work2
  REAL(kind=8), DIMENSION(:), ALLOCATABLE :: tau
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: prop_buf

  !-----------------------------------------------------!
  !                                                     !
  ! End of preamble                                     !
  !                                                     !
  !-----------------------------------------------------!


CONTAINS
 !-----------------------------------------------------!
 !                                                     !
 ! Function declarations                               !
 !                                                     !
 !-----------------------------------------------------!

 !> Initialize a matrix A of dimension (ndim,ndim) as a unit matrix
 SUBROUTINE init_one(A)
    REAL(KIND=8), dimension(ndim,ndim),INTENT(INOUT) :: A
    INTEGER :: i

    A=0.0d0
    DO i=1,ndim
      A(i,i)=1.0d0
    END DO

  END SUBROUTINE init_one

  !> initialize Lyapunov computation (possibly also vectors in later version)
  !> initializes also a random orthogonal matrix for the matrix ensemble. 
  SUBROUTINE init_lyap
    INTEGER :: AllocStat,ilaenv,info
    lwork=ilaenv(1,"dgeqrf"," ",ndim,ndim,ndim,-1)
    lwork=ndim*lwork
    ALLOCATE(prop_buf(ndim,ndim),lyapunov(ndim),loclyap(ndim),ensemble(ndim,ndim),tau(ndim),prop(ndim,ndim), &
    & work2(ndim),work(lwork),STAT=AllocStat) 
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    lyapunov=0.0d0
    loclyap=0.0d0
    CALL init_one(prop)
    CALL random_number(ensemble)
    CALL DGEQRF(ndim,ndim,ensemble,ndim,tau,work,lwork, info) ! qr decomposition
  END SUBROUTINE init_lyap

  !> Multiplies prop_int from the left with the prop matrix defined in this
  !> module and saves theresult to prop_int
  SUBROUTINE multiply_prop(prop_mul)
    REAL(KIND=8), dimension(ndim,ndim),intent(in) :: prop_mul
    REAL(KIND=8), dimension(ndim,ndim) :: prop_buf
    prop_buf=prop    
    CALL DGEMM ('n', 'n', ndim, ndim, ndim, 1.0d0, prop_mul, ndim,prop_buf, ndim,0.0d0, prop, ndim)
    
  END SUBROUTINE

  !> Performs the bennettin step in integration. Multiplies the aggregated
  !> propagators in prop with ensemble and performs QR decomposition (Gram-Schmidt
  !> orthogonalization gives Q and upper triangular matrix R). Computes also the
  !> Lyapunov exponents via the diagonal of R. WATCH OUT: prop is changed during
  !> the subroutine and NOT restored
  SUBROUTINE bennettin_step
    
    INTEGER :: info,k

    ! multiply the Propagator from the right side with the non transposed q matrix from the qr decomposition
    CALL DORM2R("r","n",ndim,ndim,ndim,ensemble,ndim,tau,prop,ndim,work2,info)
    ! prop contains now prop*ensemble => copy to ensemble 
    ensemble=prop

    ! from here on ensemble contains the new information prop*ensemble
    CALL DGEQRF(ndim,ndim,ensemble,ndim,tau,work,lwork, info) ! qr decomposition
    
    DO k=1,ndim

      loclyap(k)=log(abs(ensemble(k,k)))/tw
    END DO

    !
    ! Add here save for 
    !

    !lyapunov=lyapunov+loclyap
    END
    
END MODULE lyap_vectors
     
