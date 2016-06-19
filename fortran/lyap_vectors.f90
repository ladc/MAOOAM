
! lyap_vectors.f90
!
!> Module for computation of Lyapunov exponents and vectors
!
!> @copyright                                                               
!> 2016 Sebastian Schubert.
!> See LICENSE.txt for license information.
!
!---------------------------------------------------------------------------!
!                                                                           
!>  @remark                                                                 
!>  This module contains the necessary tools to perform the Benettin
!>  steps to compute the lyapunov_BLV exponents. (Ginelli for CLV will be added later)
!>
!>  References :
!>  Benettin, G., Galgani, L., Giorgilli, A., & Strelcyn, J. M. (1980). Lyapunov
!>  characteristic exponents for smooth dynamical systems; a method for computing
!>  all of them. Part 2: Numerical application. \a Meccanica \a, 15, 21-30.
!                                                                           
!---------------------------------------------------------------------------


MODULE lyap_vectors

  !-----------------------------------------------------!
  !                                                     !
  ! Preamble and variables declaration                  !
  !                                                     !
  !-----------------------------------------------------!

  USE params, only: ndim,dt,rescaling_time, compute_BLV, conv_BLV,compute_FLV, conv_FLV,compute_CLV, conv_CLV,length_lyap,offset 
  USE util, only: init_one
  IMPLICIT NONE

  PRIVATE
  
  PUBLIC :: benettin_step,loclyap_BLV,lyapunov_BLV,ensemble,init_lyap,multiply_prop,get_prop
 
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: loclyap_BLV    !< Buffer containing the local Lyapunov exponenti of BLV
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: lyapunov_BLV   !< Buffer containing the averaged Lyapunov exponent of BLV
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: loclyap_FLV    !< Buffer containing the local Lyapunov exponent of FLV
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: lyapunov_FLV   !< Buffer containing the averaged Lyapunov exponent Ff FLV
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: loclyap_CLV    !< Buffer containing the local Lyapunov exponent of CLV
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: lyapunov_CLV   !< Buffer containing the averaged Lyapunov exponent of CLV
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: ensemble !< Buffer containing the QR decompsoition of the ensemble
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: prop     !< Buffer holding the propagator matrix
  
  INTEGER :: lwork
  REAL(kind=8), DIMENSION(:), ALLOCATABLE :: work       !< Temporary buffer for QR decomposition
  REAL(kind=8), DIMENSION(:), ALLOCATABLE :: work2      !< Temporary buffer for QR decomposition
  REAL(kind=8), DIMENSION(:), ALLOCATABLE :: tau        !< Temporary buffer for QR decomposition
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: prop_buf !< Buffer holding the local propagator matrix

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


  !> Initialize Lyapunov computation (possibly also vectors in later version)
  !> and initializes also a random orthogonal matrix for the matrix ensemble. 
  SUBROUTINE init_lyap
    INTEGER :: AllocStat,ilaenv,info
    lwork=ilaenv(1,"dgeqrf"," ",ndim,ndim,ndim,-1)
    lwork=ndim*lwork
    ALLOCATE(prop_buf(ndim,ndim),lyapunov_BLV(ndim),loclyap_BLV(ndim),ensemble(ndim,ndim),tau(ndim),prop(ndim,ndim), &
    & work2(ndim),work(lwork),STAT=AllocStat) 
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    lyapunov_BLV=0.0d0
    loclyap_BLV=0.0d0
    CALL init_one(prop)
    IF (compute_BLV .OR. compute_FLV .OR. compute_CLV) THEN
      CALL random_number(ensemble)
      ensemble=2*(ensemble-0.5)
      CALL DGEQRF(ndim,ndim,ensemble,ndim,tau,work,lwork, info) ! qr decomposition
    END IF
  END SUBROUTINE init_lyap

  !> Multiplies prop_mul from the left with the prop matrix defined in this
  !> module and saves the result to prop_mul
  !> @param prop_mul local propagator to multiply with the global one
  SUBROUTINE multiply_prop(prop_mul)
    REAL(KIND=8), DIMENSION(ndim,ndim),INTENT(IN) :: prop_mul
    prop_buf=prop    
    CALL DGEMM ('n', 'n', ndim, ndim, ndim, 1.0d0, prop_mul, ndim,prop_buf, ndim,0.0d0, prop, ndim)
  END SUBROUTINE multiply_prop

  !> Performs the benettin step in integration. Multiplies the aggregated
  !> propagators in prop with ensemble and performs QR decomposition (Gram-Schmidt
  !> orthogonalization gives Q and upper triangular matrix R). Computes also the
  !> Lyapunov exponents via the diagonal of R. WATCH OUT: prop is changed during
  !> the subroutine and restored to a unit matrix
  SUBROUTINE benettin_step
    INTEGER :: info,k

    ! Multiply the Propagator prop from the right side with the non transposed q matrix
    ! from the qr decomposition which is stored in ensemble.
    CALL DORM2R("r","n",ndim,ndim,ndim,ensemble,ndim,tau,prop,ndim,work2,info)
    ! prop contains prop*ensemble but QR decomposed(tau is needed for that as
    ! well !) => copy to ensemble 
    ensemble=prop

    ! From here on ensemble contains the new information prop*ensemble
    CALL DGEQRF(ndim,ndim,ensemble,ndim,tau,work,lwork, info) ! qr decomposition
    
    DO k=1,ndim
      loclyap_BLV(k)=log(abs(ensemble(k,k)))/rescaling_time
    END DO

    !
    ! Add here save for 
    !

    ! Initialise prop again with unit matrix
    CALL init_one(prop) 
    
   END SUBROUTINE benettin_step

   !> Routine that returns the current global propagator
   SUBROUTINE get_prop(prop_ret)
     REAL(KIND=8), DIMENSION(ndim,ndim),INTENT(OUT) :: prop_ret
     prop_ret=prop
   END SUBROUTINE get_prop

END MODULE lyap_vectors
     
