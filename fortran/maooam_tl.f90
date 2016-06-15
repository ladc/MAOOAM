
!  maooam_tl.f90
!
!> Fortran 90 implementation of the modular arbitrary-order ocean-atmosphere
!> model MAOOAM allowing for the simultaneous integration of both the nonlinear
!> and tangent linear model.
!
!> @copyright                                                               
!> 2016 Lesley De Cruz & Jonathan Demaeyer.
!> See LICENSE.txt for license information.                                  
!
!---------------------------------------------------------------------------!

PROGRAM maooam_tl
  USE params, only: ndim, dt, tw, t_trans, t_run, writeout
  USE aotensor_def, only: init_aotensor
  USE maooam_tl_ad, only: init_tltensor
  USE IC_def, only: load_IC, IC
  USE ICdelta_def, only: load_ICdelta, ICdelta
  USE integrator, only: init_integrator,step
  USE rk4_tl_ad_integrator, only: init_tl_ad_integrator,evolve_tl_step
  USE stat
  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: X       !< State variable in the model
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: Xnew    !< Updated state variable
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: deltaX       !< State variable in the model
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: deltaXnew    !< Updated state variable
  REAL(KIND=8) :: t=0.D0                             !< Time variable

  PRINT*, 'Model MAOOAM v1.0'
  PRINT*, '       - Simultaneous integration of the original and TL models'
  PRINT*, 'Loading information...'

  CALL init_aotensor    ! Compute the tensors
  CALL init_tltensor   

  CALL load_IC          ! Load the initial condition
  CALL load_ICdelta

  CALL init_integrator  ! Initialize the integrator
  CALL init_tl_ad_integrator  ! Initialize the integrator

  IF (writeout) THEN
     OPEN(10,file='evol_field.dat')
     OPEN(11,file='evol_tl_field.dat')
  END IF

  ALLOCATE(X(0:ndim),Xnew(0:ndim),deltaX(0:ndim),deltaXnew(0:ndim))

  X=IC
  deltaX=ICdelta
  PRINT*, 'Starting the transient time evolution...'

  DO WHILE (t<t_trans)
     CALL step(X,t,dt,Xnew)
     X=Xnew
  END DO

  PRINT*, 'Starting the time evolution...'

  CALL init_stat
  
  t=0.D0

  DO WHILE (t<t_run)
     CALL evolve_tl_step(X,deltaX,t,dt,Xnew,deltaXnew)
     X=Xnew
     deltaX=deltaXnew
     IF (mod(t,tw)<dt) THEN
        IF (writeout) WRITE(10,*) t,X(1:ndim)
        IF (writeout) WRITE(11,*) t,deltaX(1:ndim)
        CALL acc(X)
     END IF
  END DO

  PRINT*, 'Evolution finished.'

  IF (writeout) CLOSE(10)
  IF (writeout) CLOSE(11)

  IF (writeout) OPEN(10,file='mean_field.dat')

  X=mean()
  IF (writeout) WRITE(10,*) X(1:ndim)
  IF (writeout) CLOSE(10)

END PROGRAM maooam_tl
