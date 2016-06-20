
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

  USE params, only: ndim,dt,rescaling_time, compute_BLV,compute_BLV_LE, conv_BLV,compute_FLV,compute_FLV_LE, conv_FLV,compute_CLV_LE,compute_CLV, length_lyap,offset 
  USE util, only: init_one
  IMPLICIT NONE

  PRIVATE
  
  PUBLIC :: benettin_step,ginelli,ensemble,init_lyap,multiply_prop,get_prop,compute_vectors,compute_exponents
  PUBLIC :: loclyap_BLV,lyapunov_BLV,loclyap_FLV,lyapunov_FLV,loclyap_CLV,lyapunov_CLV, init_ensemble
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: loclyap_BLV    !< Buffer containing the local Lyapunov exponenti of BLV
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: lyapunov_BLV   !< Buffer containing the averaged Lyapunov exponent of BLV
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: loclyap_FLV    !< Buffer containing the local Lyapunov exponent of FLV
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: lyapunov_FLV   !< Buffer containing the averaged Lyapunov exponent Ff FLV
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: loclyap_CLV    !< Buffer containing the local Lyapunov exponent of CLV
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: lyapunov_CLV   !< Buffer containing the averaged Lyapunov exponent of CLV
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: R          !< Upper triangular propagator in packed storage
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: BLV      !< Buffer containing the Backward Lyapunov Vectors
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: FLV      !< Buffer containing the Forward Lyapunov Vectors
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: CLV      !< Buffer containing the Covariant Lyapunov Vectors
   
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: ensemble !< Buffer containing the QR decompsoition of the ensemble
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: prop     !< Buffer holding the propagator matrix
  
  INTEGER :: lwork
  REAL(kind=8), DIMENSION(:), ALLOCATABLE :: work       !< Temporary buffer for QR decomposition
  REAL(kind=8), DIMENSION(:), ALLOCATABLE :: work2      !< Temporary buffer for QR decomposition
  REAL(kind=8), DIMENSION(:), ALLOCATABLE :: tau        !< Temporary buffer for QR decomposition
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: prop_buf !< Buffer holding the local propagator matrix
 
  INTEGER, DIMENSION(:), ALLOCATABLE :: IPIV                      !< Necessary for dgetrs
  
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: one  
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
    ALLOCATE(one(ndim,ndim))
    CALL init_one(one)
    lwork=ilaenv(1,"dgeqrf"," ",ndim,ndim,ndim,-1)
    lwork=ndim*lwork
    ALLOCATE(prop_buf(ndim,ndim),ensemble(ndim,ndim),tau(ndim),prop(ndim,ndim), &
    & work2(ndim),work(lwork),STAT=AllocStat) 
    IF (compute_BLV .OR. compute_BLV_LE) ALLOCATE(BLV(ndim,ndim))
    IF (compute_BLV_LE) THEN
      ALLOCATE(lyapunov_BLV(ndim),loclyap_BLV(ndim));loclyap_BLV=0.D0;lyapunov_BLV=0.D0
    END IF
    IF (compute_FLV .OR. compute_FLV_LE) ALLOCATE(FLV(ndim,ndim),IPIV(ndim))
    IF (compute_FLV_LE) THEN 
      ALLOCATE(lyapunov_FLV(ndim),loclyap_FLV(ndim));loclyap_FLV=0.D0;lyapunov_FLV=0.D0
    END IF
    IF (compute_CLV .OR. compute_CLV_LE) ALLOCATE(CLV(ndim,ndim),R(ndim*(ndim+1)/2))
    IF (compute_CLV_LE) THEN
      ALLOCATE(lyapunov_CLV(ndim),loclyap_CLV(ndim));loclyap_CLV=0.D0;lyapunov_CLV=0.D0
    END IF
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
    
    CALL init_ensemble 
    IF (compute_CLV .OR. compute_CLV_LE) THEN
      CALL random_number(CLV)
      CALL normMAT(CLV,loclyap_CLV)
      loclyap_CLV=-log(abs(loclyap_CLV))/rescaling_time
    END IF
  END SUBROUTINE init_lyap
 
  !> Routine to initialise ensmble. Will be called externally (therefore
  !> separate routine) 
  SUBROUTINE init_ensemble
    INTEGER :: info 
    CALL init_one(prop)
    CALL random_number(ensemble)
    ensemble=2*(ensemble-0.5)
    CALL DGEQRF(ndim,ndim,ensemble,ndim,tau,work,lwork, info) ! qr decomposition
  END SUBROUTINE init_ensemble

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
  SUBROUTINE benettin_step(forward,step)
    INTEGER :: info,k,step
    LOGICAL :: forward 
    
    IF (.NOT. forward) THEN
      CALL read_lyapvec(step,prop,23)
      CALL DGETRF(ndim,ndim,prop,ndim,IPIV,info)
      CALL DGETRS('n',ndim,ndim,prop,ndim,IPIV,one,ndim,info)
    ENDIF
    ! Multiply the Propagator prop from the right side with the non transposed q matrix
    ! from the qr decomposition which is stored in ensemble.
    CALL DORM2R("r","n",ndim,ndim,ndim,ensemble,ndim,tau,prop,ndim,work2,info)
    ! prop contains prop*ensemble but QR decomposed(tau is needed for that as
    ! well !) => copy to ensemble 
    ensemble=prop
    
    ! From here on ensemble contains the new information prop*ensemble
    CALL DGEQRF(ndim,ndim,ensemble,ndim,tau,work,lwork, info) ! qr decomposition

    IF (forward) THEN
      IF (compute_BLV_LE) THEN
        DO k=1,ndim
          loclyap_BLV(k)=log(abs(ensemble(k,k)))/rescaling_time
        END DO
      END IF
    ELSE
      IF (compute_FLV_LE) THEN
        DO k=1,ndim
          loclyap_FLV(k)=-log(abs(ensemble(k,k)))/rescaling_time
        END DO
      END IF
    END IF
    
    ! Note propagator is not restored!
    
   END SUBROUTINE benettin_step
   
   !> This routine performs the backward ginelli step
   SUBROUTINE ginelli(step)
   INTEGER :: step,info
      READ(unit=33,rec=step) R
      CALL DTPTRS('u','n','n',ndim,ndim,R,CLV,ndim,info)
      CALL normMAT(CLV,loclyap_CLV)
      loclyap_CLV=-log(abs(loclyap_CLV))/rescaling_time
   END SUBROUTINE ginelli

   !> Routine that returns the current global propagator
   SUBROUTINE get_prop(prop_ret)
     REAL(KIND=8), DIMENSION(ndim,ndim),INTENT(OUT) :: prop_ret
     prop_ret=prop
   END SUBROUTINE get_prop
   
   !> Routine that saves the BLV, FLV and CLV if in right time period according
   !> to namelist parameters in int_params.nml
   SUBROUTINE compute_vectors(t,step,forward)
   INTEGER :: step
   REAL(KIND=8) :: t
   LOGICAL :: forward
   LOGICAL :: past_conv_BLV,before_conv_FLV
   INTEGER :: info

   past_conv_BLV=(t-offset.gt.conv_BLV)
   before_conv_FLV=(t-offset.lt.length_lyap-conv_FLV)
   IF (past_conv_BLV) THEN
     IF (forward) THEN

       IF (compute_BLV .AND. before_conv_FLV) THEN
         BLV=ensemble ! make copy of QR decomposed ensemble     
         CALL DORGQR(ndim,ndim,ndim,BLV,ndim,tau,work,lwork,info) !retrieve Q (BLV) matrix 
         CALL write_lyapvec(step+1,BLV,12) !write Q (BLV) matrix
       END IF

       IF (compute_FLV .OR. compute_FLV_LE) call write_lyapvec(step,prop,23)

       IF (compute_CLV .OR. compute_CLV_LE) THEN
         CALL packTRI(ensemble,R)
         WRITE(unit=33,rec=step) R
       END IF
       CALL init_one(prop)

     ELSE
       IF (compute_FLV .AND. before_conv_FLV) THEN
         FLV=ensemble ! make copy of QR decomposed ensemble     
         CALL DORGQR(ndim,ndim,ndim,FLV,ndim,tau,work,lwork,info) !retrieve Q (BLV) matrix 
         CALL write_lyapvec(step,FLV,22)
       END IF

       IF (compute_CLV .AND. before_conv_FLV) THEN
         CALL read_lyapvec(step,BLV,12)
         CALL DGEMM ('n', 'n', ndim, ndim,ndim, 1.0d0, BLV, ndim,CLV, ndim,0.D0,FLV,ndim) !FLV is overwritten
         CALL write_lyapvec(step,FLV,32) ! THIS is the CLV not the FLV
       END IF

     END IF
   END IF   
   END SUBROUTINE compute_vectors
   
   SUBROUTINE compute_exponents(t,step,forward)
   INTEGER :: step
   REAL(KIND=8) :: t
   LOGICAL :: forward
   LOGICAL :: past_conv_BLV,before_conv_FLV
   INTEGER :: info

   past_conv_BLV=(t-offset.gt.conv_BLV)
   before_conv_FLV=(t-offset.lt.length_lyap-conv_FLV)
   IF (past_conv_BLV) THEN
     IF (forward) THEN
     
         IF (compute_BLV_LE .AND. before_conv_FLV) WRITE(11,rec=step) loclyap_BLV
   
     ELSE
       IF (compute_FLV .AND. before_conv_FLV) WRITE(21,rec=step) loclyap_FLV

       IF (compute_CLV .AND. before_conv_FLV) WRITE(31,rec=step) loclyap_CLV

     END IF
   END IF      
   END SUBROUTINE compute_exponents
      
   !> Routine to read lyapunov vectors
   SUBROUTINE read_lyapvec(i,vectors,unitI)
   INTEGER :: i,unitI
   REAL(KIND=8), DIMENSION(ndim,ndim), INTENT(OUT) :: vectors
     READ(unit=unitI,rec=i) vectors
   END SUBROUTINE read_lyapvec
  
   !> Routine to write lyapunov vectors
   SUBROUTINE write_lyapvec(i,vectors,unitI)
   INTEGER :: i,unitI
   REAL(KIND=8), DIMENSION(ndim,ndim), INTENT(IN) :: vectors
     write(unit=unitI,rec=i) vectors
   END SUBROUTINE write_lyapvec

   !> Routine that normalizes uppertriangular packed storage matrix (LAPACK
   !> standard)
   SUBROUTINE normMAT(M,norms) 
   INTEGER :: i
   REAL(KIND=8), dimension(ndim,ndim), INTENT(INOUT) :: M
   REAL(KIND=8), dimension(ndim) :: norms
   DO i=1,ndim
     norms(i)=sqrt(sum(M(:,i)**2.0d0)) 
     M(:,i)=M(:,i)/norms(i)
   END DO
   END SUBROUTINE normMAT

   !> Routine that normalizes uppertriangular packed storage matrix (LAPACK
   !> standard)
   SUBROUTINE normTRI(packedR) 
   INTEGER :: k,j,i
   REAL(KIND=8), dimension(ndim*(ndim+1)/2), INTENT(INOUT) :: packedR
   k = 0
   DO j=1,ndim
     packedR(k+1:k+j)=packedR(k+1:k+j)/sqrt(sum(packedR(k+1:k+j)**2.0d0))
     k = k+j
   END DO
   END SUBROUTINE normTRI

   !> Routine that transforms uppertriangular part into packed storage (LAPACK
   !> standard)
   SUBROUTINE packTRI(M,packedR) 
   INTEGER :: k,j,i
   REAL(KIND=8), dimension(ndim,ndim), INTENT(IN) :: M
   REAL(KIND=8), dimension(ndim*(ndim+1)/2), INTENT(OUT) :: packedR
   k = 0
   DO j=1,ndim
     DO i=1,j
       k = k+1
       packedR(k)=M(i,j)
     END DO
   END DO
   END SUBROUTINE packTRI
  
   !> Routine that transforms uppertriangular part into normal storage (LAPACK
   !> standard)
   SUBROUTINE unpackTRI(packedR,M) 
   INTEGER :: k,j,i
   REAL(KIND=8), dimension(ndim,ndim), INTENT(OUT) :: M
   REAL(KIND=8), dimension(ndim*(ndim+1)/2), INTENT(IN) :: packedR
   M=0.D0
   k = 0
   DO j=1,ndim
     DO i=1,j
       k = k+1
       M(i,j)=packedR(k)
     END DO
   END DO
   END SUBROUTINE unpackTRI
END MODULE lyap_vectors
     
