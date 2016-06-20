
!  maooam_lyapvectors.f90
!
!> Fortran 90 implementation of the modular arbitrary-order ocean-atmosphere
!> model MAOOAM computing the Lyapunov spectrum.
!
!> @copyright                                                               
!> 2016 Lesley De Cruz, Sebastian Schubert & Jonathan Demaeyer.
!> See LICENSE.txt for license information.                                  
!
!---------------------------------------------------------------------------!

PROGRAM maooam_lyapvectors
  USE params, only: ndim, dt, tw, t_trans, t_run, writeout, rescaling_time, compute_BLV_LE, compute_BLV, conv_BLV,compute_FLV, compute_FLV_LE, conv_FLV,compute_CLV, compute_CLV_LE,length_lyap,offset 
  USE aotensor_def, only: init_aotensor
  USE maooam_tl_ad, only: init_tltensor
  USE IC_def, only: load_IC, IC
  USE integrator, only: init_integrator,step
  USE tl_ad_integrator, only: init_tl_ad_integrator,prop_step
  USE lyap_vectors, only:  lyapunov_BLV,loclyap_BLV,lyapunov_FLV,loclyap_FLV,lyapunov_CLV,loclyap_CLV,init_lyap,multiply_prop,benettin_step,compute_vectors,compute_exponents,init_ensemble,ginelli
  USE stat
  USE lyap_stat
  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: X             !< State variable in the model
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: Xnew          !< Updated state variable
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: prop_buf    !< Buffer for Integrator propagator
  REAL(KIND=8) :: t=0.D0                                   !< Time variable
  REAL(KIND=8) :: resc=1.D-9                               !< Variable rescaling factor for the divergence method
  REAL(KIND=8) :: t_up
  INTEGER :: IndexBen,WRSTAT
  CHARACTER(LEN=19) :: FMTX

  PRINT*, 'Model MAOOAM v1.0'
  PRINT*, '      - with computation of the Lyapunov spectrum'
  PRINT*, 'Loading information...'

  CALL init_aotensor    ! Compute the tensors
  CALL init_tltensor   
  CALL load_IC          ! Load the initial condition

  CALL init_integrator  ! Initialize the integrator
  CALL init_tl_ad_integrator  ! Initialize tangent linear integrator
  CALL init_lyap        ! Initialize Lyapunov computation
  write(FMTX,'(A10,i3,A6)') '(F10.2,4x,',ndim,'E15.5)'
  t_up=dt/t_trans*100.D0

  IF (writeout) THEN
     OPEN(10,file='evol_field.dat')
     IF (compute_BLV .OR. compute_CLV_LE) OPEN(12,file='BLV_vec.dat',status='replace',form='UNFORMATTED',access='DIRECT',recl=8*ndim**2)
     IF (compute_BLV_LE) OPEN(11,file='BLV_exp.dat',status='replace',form='UNFORMATTED',access='DIRECT',recl=8*ndim)
     IF (compute_FLV) OPEN(22,file='FLV_vec.dat',status='replace',form='UNFORMATTED',access='DIRECT',recl=8*ndim**2)
     IF (compute_FLV_LE) OPEN(21,file='FLV_exp.dat',status='replace',form='UNFORMATTED',access='DIRECT',recl=8*ndim)
     IF (compute_FLV .OR. compute_FLV_LE) OPEN(23,file='propagator.dat',status='replace',form='UNFORMATTED',access='DIRECT',recl=8*ndim**2)
     IF (compute_CLV) OPEN(32,file='CLV_vec.dat',status='replace',form='UNFORMATTED',access='DIRECT',recl=8*ndim**2)
     IF (compute_CLV_LE) OPEN(31,file='CLV_exp.dat',status='replace',form='UNFORMATTED',access='DIRECT',recl=8*ndim)
     IF (compute_CLV .OR. compute_CLV_LE) OPEN(33,file='R.dat',status='replace',form='UNFORMATTED',access='DIRECT',recl=8*ndim*(ndim+1)/2)
  END IF

  ALLOCATE(X(0:ndim),Xnew(0:ndim),prop_buf(ndim,ndim))
  X=IC
  PRINT*, 'Starting the transient time evolution... t_trans = ',t_trans

  DO WHILE (t<t_trans)
     CALL step(X,t,dt,Xnew)
     X=Xnew
     IF (mod(t/t_trans*100.D0,0.1)<t_up) WRITE(*,'(" Progress ",F6.1," %",A,$)') t/t_trans*100.D0,char(13)
  END DO

  PRINT*, 'Starting the forward time evolution... t_run = ',t_run,'; offset = ',offset

  CALL init_stat
  CALL lyap_init_stat
  t=0.D0
  t_up=dt/t_run*100.D0
  
  !
  ! Possibility for Offset for Lyapunov Test (TODO: implement read background
  ! trajectory if offset is bigger than zero, new module
  ! read_(rk4)_tl_ad_integrator.f90
  ! that replaces prop_step function with read command for X)
  !
  
  t=offset
    
  !
  ! Start forward part of run of run
  !
  IndexBen=0 ! Index for lyapunov vector calculations
  DO WHILE (t<offset+length_lyap)

     CALL prop_step(X,prop_buf,t,dt,Xnew,.false.) ! Obtains propagator prop_buf at X
     CALL multiply_prop(prop_buf) ! Multiplies prop_buf with prop
     X=Xnew
     IF (mod(t,rescaling_time)<dt) THEN
        IndexBen=IndexBen+1
        CALL benettin_step(.true.,IndexBen) ! Performs QR step with prop
        CALL compute_exponents(t,IndexBen,.true.)
        !CALL lyap_acc(loclyap_BLV)
        CALL acc(X)
        CALL compute_vectors(t,IndexBen,.true.)
     END IF
     IF (mod(t,tw)<dt) THEN
        !! Uncomment if you want the trajectory (may generate a huge file!)
        IF (writeout) WRITE(10,FMTX) t,X(1:ndim) 
        CONTINUE
     END IF
   
     IF (mod(t/t_run*100.D0,0.1)<t_up) WRITE(*,'(" Progress ",F6.1," %",A,$)') t/t_run*100.D0,char(13)
  END DO
  PRINT*, 'Forward evolution finished.'

  PRINT*, 'Starting the backward evolution ...'
  IF (compute_FLV .OR. compute_FLV_LE) THEN
    CALL init_ensemble
  END IF

  DO WHILE (t>offset .AND. IndexBen>0)

    IF (compute_FLV .OR. compute_FLV_LE) CALL benettin_step(.false.,IndexBen) ! Performs QR step with prop
    IF (compute_CLV .OR. compute_CLV_LE) CALL ginelli(IndexBen)               ! Performs Ginelli step with prop
    CALL compute_exponents(t,IndexBen,.false.)
    CALL compute_vectors(t,IndexBen,.false.)
    IndexBen=IndexBen-1
    t=t-dt
    IF (mod(t/t_run*100.D0,0.1)<t_up) WRITE(*,'(" Progress ",F6.1," %",A,$)') t/t_run*100.D0,char(13)
  END DO

!  IF (writeout) THEN
!     OPEN(10,file='mean_lyapunov.dat')
!     lyapunov_BLV=lyap_mean()
!     WRITE(10,*) 'mean',lyapunov(1:ndim)
!     lyapunov_BLV=lyap_var()
!     WRITE(10,*) 'var',lyapunov(1:ndim)
!  END IF

  IF (writeout) THEN
     OPEN(10,file='mean_field.dat')
     X=mean()
     WRITE(10,*) 'mean',X(1:ndim)
     X=var()
     WRITE(10,*) 'var',X(1:ndim)
  END IF

END PROGRAM maooam_lyapvectors
