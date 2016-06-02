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
  USE IFPORT, only: rand
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
    REAL :: seed
    lwork=ilaenv(1,"dgeqrf"," ",ndim,ndim,ndim,-1)
    lwork=ndim*lwork
    ALLOCATE(prop_buf(ndim,ndim),lyapunov(ndim),loclyap(ndim),ensemble(ndim,ndim),tau(ndim),prop(ndim,ndim), &
    & work2(ndim),work(lwork),STAT=AllocStat) 
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    lyapunov=0.0d0
    loclyap=0.0d0
    CALL init_one(prop)
    CALL CPU_TIME(seed)
    CALL init_random_seed()
    CALL random_number(ensemble)
    CALL DGEQRF(ndim,ndim,ensemble,ndim,tau,work,lwork, info) ! qr decomposition
  END SUBROUTINE init_lyap

  !> Multiplies prop_int from the left with the prop matrix defined in this
  !> module and saves theresult to prop_int
  SUBROUTINE multiply_prop(prop_mul)
    REAL(KIND=8), dimension(ndim,ndim),intent(in) :: prop_mul
    prop_buf=prop    
    CALL DGEMM ('n', 'n', ndim, ndim, ndim, 1.0d0, prop_mul, ndim,prop_buf, ndim,0.0d0, prop, ndim)
    
  END SUBROUTINE multiply_prop

  !> Performs the bennettin step in integration. Multiplies the aggregated
  !> propagators in prop with ensemble and performs QR decomposition (Gram-Schmidt
  !> orthogonalization gives Q and upper triangular matrix R). Computes also the
  !> Lyapunov exponents via the diagonal of R. WATCH OUT: prop is changed during
  !> the subroutine and restored to a unit matrix
  SUBROUTINE bennettin_step
    
    INTEGER :: info,k

    ! > multiply the Propagator prop from the right side with the non transposed q matrix from the qr decomposition
    ! > which is stored in ensemble.
    CALL DORM2R("r","n",ndim,ndim,ndim,ensemble,ndim,tau,prop,ndim,work2,info)
    ! > prop contains prop*ensemble but QR decomposed(tau is needed for that as
    ! > well !) => copy to ensemble 
    ensemble=prop

    ! from here on ensemble contains the new information prop*ensemble
    CALL DGEQRF(ndim,ndim,ensemble,ndim,tau,work,lwork, info) ! qr decomposition
    
    DO k=1,ndim
      loclyap(k)=log(abs(ensemble(k,k)))/tw
    END DO

    !
    ! Add here save for 
    !

    !> Initialise prop again with unit matrix
    CALL init_one(prop) 
    
   END SUBROUTINE bennettin_step
   SUBROUTINE init_random_seed()
     USE iso_fortran_env, only: int64
     USE IFPORT, only: getpid
     IMPLICIT NONE
     INTEGER, ALLOCATABLE :: seed(:)
     INTEGER :: i, n, un, istat, dt(8), pid
     INTEGER(int64) :: t
         
     CALL random_seed(size = n)
     ALLOCATE(seed(n))
     ! First try IF the OS provides a random number generator
     OPEN(newunit=un, file="/dev/urandom", access="stream", &
          form="unformatted", action="read", status="old", iostat=istat)
     IF (istat == 0) THEN
        READ(un) seed
        CLOSE(un)
     ELSE
        ! Fallback to XOR:ing the current time and pid. The PID is
        ! useful in case one launches multiple instances of the same
        ! program in parallel.
        CALL system_clock(t)
        IF (t == 0) THEN
           CALL date_and_time(values=dt)
           t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
                + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
                + dt(3) * 24_int64 * 60 * 60 * 1000 &
                + dt(5) * 60 * 60 * 1000 &
                + dt(6) * 60 * 1000 + dt(7) * 1000 &
                + dt(8)
        END IF
        pid = getpid()
        t = ieor(t, int(pid, kind(t)))
        DO i = 1, n
           seed(i) = lcg(t)
        END DO
     END IF
     CALL random_seed(put=seed)
   contains
     ! This simple PRNG might not be good enough for real work, but is
     ! sufficient for seeding a better PRNG.
     FUNCTION lcg(s)
       integer :: lcg
       integer(int64) :: s
       IF (s == 0) THEN
          s = 104729
       ELSE
          s = mod(s, 4294967296_int64)
       END IF
       s = mod(s * 279470273_int64, 4294967291_int64)
       lcg = int(mod(s, int(huge(0), int64)), kind(0))
     END FUNCTION lcg
   END SUBROUTINE init_random_seed

END MODULE lyap_vectors
     
