
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
  USE params, only: ndim, dt, tw, t_trans, t_run, writeout, rescaling_time, sampling_time, compute_BLV_LE,&
 & compute_BLV, conv_BLV,compute_FLV, compute_FLV_LE, conv_FLV,compute_CLV, compute_CLV_LE,length_lyap,offset 
  USE aotensor_def, only: init_aotensor
  USE tl_ad_tensor, only: init_tltensor
  USE IC_def_ext, only: load_IC, IC, write_IC
  USE integrator, only: init_integrator,step
  USE tl_ad_integrator, only: init_tl_ad_integrator,prop_step
  USE lyap_vectors, only:  lyapunov_BLV,loclyap_BLV,lyapunov_FLV,loclyap_FLV,lyapunov_CLV, &
 & loclyap_CLV,init_lyap,multiply_prop,benettin_step,compute_vectors,compute_exponents,init_ensemble,ginelli,get_lyap_state
  USE stat
  USE lyap_stat
  USE util, only:str
  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: X             !< State variable in the model
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: Xnew          !< Updated state variable
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: prop_buf    !< Buffer for Integrator propagator
  REAL(KIND=8) :: t=0.D0                                   !< Time variable
  REAL(KIND=8) :: resc=1.D-9                               !< Variable rescaling factor for the divergence method
  REAL(KIND=8) :: t_up
  INTEGER :: IndexBen,IndexSample,WRSTAT
  LOGICAL :: write_sample
  CHARACTER(LEN=20) :: FMTX

  PRINT*, 'Model MAOOAM v1.0'
  PRINT*, '      - with computation of the Lyapunov spectrum'
  PRINT*, 'Loading information...'

  CALL init_aotensor    ! Compute the tensors
  CALL init_tltensor   
  CALL load_IC          ! Load the initial condition

  CALL init_integrator  ! Initialize the integrator
  CALL init_tl_ad_integrator  ! Initialize tangent linear integrator
  CALL init_lyap        ! Initialize Lyapunov computation and open files for output of LVs and Les
  write(FMTX,'(A10,i3,A7)') '(F10.2,4x,',ndim,'E30.15)'

  IF (writeout) OPEN(10,file='evol_field.dat')

  ALLOCATE(X(0:ndim),Xnew(0:ndim),prop_buf(ndim,ndim))
  X=IC
  
  !
  ! If offset is bigger than zero than the existing
  ! "IC_after_transient_state.nml" is used as initial condition
  !

  IF (offset>0) THEN
    t_trans=offset  ! transient phase in offset mode is just going
    PRINT*, '*** Skipping transient time evolution. ***' 
    PRINT*, '*** Instead loading IC.nml immediately and integrating until offset = ',offset,' ***'
  ELSE
    PRINT*, 'Starting transient time evolution... t_trans = ',t_trans
  ENDIF
  
  t = 0.0D0
  IF (t_trans.ne.0) t_up = dt/t_trans*100.D0
  DO WHILE (t<t_trans)
     CALL step(X,t,dt,Xnew)
     X=Xnew
     IF (mod(t/t_trans*100.D0,0.1)<t_up) WRITE(*,'(" Progress ",F6.1," %",A,$)') t/t_trans*100.D0,char(13)
  END DO

  PRINT*, 'Starting the forward time evolution with tangent linear model ... t_run = ',t_run,'; offset = ',offset

  CALL init_stat
  CALL lyap_init_stat
  t=0.D0
  t_up=dt/t_run*100.D0
  t=offset
  
  ! Possibility for Offset for Lyapunov Test (TODO: implement read background
  ! trajectory if offset is bigger than zero, new module
  ! read_(rk4)_tl_ad_integrator.f90
  ! that replaces prop_step function with read command for X
  ! COMMENT from Sebastian: I try a slightly different implementation that is
  ! more straightforward since most in between steps have to recomputed anyway
  ! we can just repeat the whole computation.
  
    
  IF (offset .eq. 0) THEN
    CALL write_IC('after_transient_state',X) 
  ELSE
    CALL write_IC('offset_'//trim(str(int(offset))),X) 
  END IF
  ! 
  ! Start forward part of run of run
  !
    
  IndexBen=Int(floor(offset/rescaling_time)) ! Index for lyapunov vector calculations
  IndexSample=Int(floor(offset/sampling_time)) ! Index for lyapunov vector calculations
  DO WHILE (t .LE. length_lyap)
     CALL prop_step(X,prop_buf,t,dt,Xnew,.false.) ! Obtains propagator prop_buf at X
     CALL multiply_prop(prop_buf) ! Multiplies prop_buf with prop
     X=Xnew

     IF (mod(t,rescaling_time)<dt .AND. t.GT.offset) THEN
        write_sample = mod(t,sampling_time)<dt
        IndexBen=IndexBen+1
        CALL benettin_step(.true.,IndexBen) ! Performs QR step with prop
        IF (write_sample) THEN
          IndexSample=IndexSample+1
          ! compute_exponents in fact only writes, so we don't need to call it every rescaling_time.
          CALL compute_exponents(t,IndexSample,.true.)
        END IF
        CALL acc(X)
        ! Because compute_vectors also accumulates prodR and computes BLVs, we need to call it every rescaling_time even though it's
        ! not always writing out the result (only at write_sample=true)
        CALL compute_vectors(t,IndexBen,IndexSample,.true.,write_sample)
     END IF
     IF (mod(t,tw)<dt) THEN
        !! Uncomment if you want the trajectory (may generate a huge file!)
        IF (writeout) WRITE(10,FMTX) t,X(1:ndim) 
     END IF
   
     IF (mod(t/t_run*100.D0,0.1)<t_up) WRITE(*,'(" Progress ",F6.1," %",A,$)') t/t_run*100.D0,char(13)
  END DO
  PRINT*, 'Forward evolution finished.'
  
  IF (offset .eq. 0) THEN
    CALL write_IC('final_state',X) 
  ELSE
    CALL write_IC('final_state_offset_'//trim(str(int(offset))),X) 
  END IF
!
! Start Backward Integration
!
  IF (compute_CLV .OR. compute_CLV_LE .OR. compute_FLV .OR. compute_FLV_LE) THEN
    PRINT*, 'Starting the backward evolution ...'
    IF (compute_FLV .OR. compute_FLV_LE) THEN
      CALL init_ensemble
    END IF
    DO WHILE (t>offset .AND. IndexBen>0 .AND. IndexSample>0)

      IF (compute_FLV .OR. compute_FLV_LE) CALL benettin_step(.false.,IndexBen) ! Performs QR step with prop
      IF (compute_CLV .OR. compute_CLV_LE) CALL ginelli(IndexSample)               ! Performs Ginelli step with prop
      write_sample = mod(t,sampling_time)<dt
      CALL compute_vectors(t,IndexBen,IndexSample,.false.,write_sample)
      IF (write_sample) THEN
        CALL compute_exponents(t,IndexSample,.false.)
        IndexSample=IndexSample-1
      END IF
      IndexBen=IndexBen-1
      t=t-rescaling_time
      IF (mod(t/t_run*100.D0,0.1)<t_up) WRITE(*,'(" Progress ",F6.1," %",A,$)') t/t_run*100.D0,char(13)
    END DO
  END IF
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
