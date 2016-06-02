
!  maooam_lyap.f90
!
!> Fortran 90 implementation of the modular arbitrary-order ocean-atmosphere
!> model MAOOAM allowing for the simultaneous integration of both the nonlinear
!> and tangent linear model. This version of the model also computes the
!> Lyapunov spectrum.
!
!> @copyright                                                               
!> 2015 Lesley De Cruz & Jonathan Demaeyer.
!> See LICENSE.txt for license information.                                  
!
!---------------------------------------------------------------------------!

PROGRAM maooam_lyap
  USE params, only: ndim, dt, tw, t_trans, t_run, writeout
  USE aotensor_def, only: init_aotensor
  USE maooam_tl_ad, only: init_tltensor
  USE IC_def, only: load_IC, IC
  USE ICdelta_def, only: load_ICdelta, ICdelta
  USE integrator, only: init_integrator,step_rk4
  USE rk4_tl_ad_integrator, only: init_tl_integrator,tl_prop_step
  USE lyap_vectors, only: lyapunov,loclyap,init_lyap,multiply_prop, & 
                               & init_one,prop,bennettin_step
  USE stat
  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: X,Xp       !< State variable in the model
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: Xnew,Xnewp    !< Updated state variable
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: prop_buf!< Buffer for Integrator propagator
  REAL(KIND=8) :: t=0.D0                             !< Time variable
  REAL(KIND=8) :: resc=1.D-9                             !< rescale variable
  INTEGER :: IndexBen,WRSTAT
  CHARACTER(LEN=19) :: FMTX

  PRINT*, 'Model MAOOAM v1.0'
  PRINT*, '      - with computation of Lyapunov spectrum'
  PRINT*, 'Loading information...'

  CALL init_aotensor    ! Compute the tensors
  CALL init_tltensor   
  CALL load_IC          ! Load the initial condition

  CALL init_integrator  ! Initialize the integrator
  CALL init_tl_integrator  ! Initialize tangent linear integrator
  CALL init_lyap ! Initialize Lyapunov computation
  write(FMTX,'(A10,i3,A6)') '(F10.2,4x,',ndim,'E15.5)'
  write(*,*) FMTX

  IF (writeout) THEN
     OPEN(10,file='evol_field.dat')
     OPEN(11,file='lyapunov_exponents.dat',status='replace',form='UNFORMATTED',access='DIRECT',recl=8*ndim)
  END IF

  ALLOCATE(X(0:ndim),Xnew(0:ndim),prop_buf(ndim,ndim))
  ALLOCATE(Xp(0:ndim),Xnewp(0:ndim))
  X=IC
  PRINT*, 'Starting the transient time evolution... t_trans = ',t_trans
 
  DO WHILE (t<t_trans)
     CALL step_rk4(X,t,dt,Xnew)
     X=Xnew
  END DO

  PRINT*, 'Starting the time evolution... t_run = ',t_run

  CALL init_stat
  Xp(0)=1.0d0
  Xp(1:ndim)=X(1:ndim)+resc/sqrt(dble(ndim))
  t=0.D0
  IndexBen=0
  DO WHILE (t<t_run)

     CALL tl_prop_step(X,prop_buf,t,dt,Xnew,.false.) !obtains propagator prop_buf at X
     CALL multiply_prop(prop_buf) ! Multiplies prop_buf with prop
     X=Xnew
     CALL step_rk4(Xp,t,dt,Xnewp)
     Xp=Xnewp;t=t-dt
     IF (mod(t,tw)<dt) THEN
        IndexBen=IndexBen+1
        CALL  bennettin_step !Performs QR step with prop
        IF (writeout) WRITE(10,FMTX) t,X(1:ndim)
        IF (writeout) WRITE(11,rec=IndexBen,iostat=WRSTAT) loclyap
        lyapunov=lyapunov+loclyap
        CALL acc(X)
        write(234,*) log(sqrt(sum((X(1:ndim)-Xp(1:ndim))**2.0d0))/resc)/tw
        Xp(1:ndim)=X(1:ndim)+(Xp(1:ndim)-X(1:ndim))/sqrt(sum((X(1:ndim)-Xp(1:ndim))**2.0d0))*resc

     END IF
  END DO
  lyapunov=lyapunov/dble(IndexBen)

  PRINT*, 'Evolution finished.'

  IF (writeout) CLOSE(10)
  IF (writeout) CLOSE(11)

  IF (writeout) OPEN(10,file='mean_lyapunov.dat')

  IF (writeout) WRITE(10,*) lyapunov(1:ndim)
  IF (writeout) CLOSE(10)

IF (writeout) OPEN(10,file='mean_field.dat')

  X=mean()
  IF (writeout) WRITE(10,*) X(1:ndim)
  IF (writeout) CLOSE(10)

END PROGRAM maooam_lyap
