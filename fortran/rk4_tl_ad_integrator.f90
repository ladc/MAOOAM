
! tl_ad_integrator.f90
!
!> Tangent Linear (TL) and Adjoint (AD) model versions of MAOOAM.
!> Integrators module.
!
!> @copyright                                                               
!> 2016 Lesley De Cruz & Jonathan Demaeyer.
!> See LICENSE.txt for license information.                                  
!
!---------------------------------------------------------------------------!
!                                                                           
!>  @remark                                                                 
!>  This module actually contains the Heun algorithm routines.
!>  The user can modify it according to its preferred integration scheme.
!>  For higher-order schemes, additional buffers will probably have to be defined.
!                                                                           
!---------------------------------------------------------------------------



MODULE rk4_tl_ad_integrator

  USE params, only: ndim
  USE tensor, only:sparse_mul3
  USE aotensor_def, only: aotensor

  USE maooam_tl_ad, only: ad,tl
  IMPLICIT NONE

  PRIVATE

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: ad_buf_y1 !< Buffer to hold the intermediate position of the adjoint model
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: ad_buf_y11 !< Buffer to hold the intermediate position of the adjoint model
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: ad_buf_kA !< Buffer to hold tendencies in the RK4 scheme for the adjoint model
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: ad_buf_kB !< Buffer to hold tendencies in the RK4 scheme for the adjoint model
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: ad_buf_kAA !< Buffer to hold tendencies in the RK4 scheme for the adjoint model
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: ad_buf_kBA !< Buffer to hold tendencies in the RK4 scheme for the adjoint model


  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: tl_buf_y1 !< Buffer to hold the intermediate position of the tangent linear model
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: tl_buf_y11 !< Buffer to hold the intermediate position of the tangent linear model
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: tl_buf_kA !< Buffer to hold tendencies in the RK4 scheme for the tangent linear model
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: tl_buf_kB !< Buffer to hold tendencies in the RK4 scheme for the tangent linear model
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: tl_buf_kAA !< Buffer to hold tendencies in the RK4 scheme for the tangent linear model
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: tl_buf_kBA !< Buffer to hold tendencies in the RK4 scheme for the tangent linear model

    
  PUBLIC :: init_ad_integrator, ad_step, init_tl_integrator, tl_step

CONTAINS

  !> Routine computing the tendencies of the nonlinear model
  !> @param t Time at which the tendencies have to be computed. Actually not needed for autonomous systems.
  !> @param y Point at which the tendencies have to be computed.
  !> @param res vector to store the result.
  !> @remark Note that it is NOT safe to pass `y` as a result buffer, 
  !> as this operation does multiple passes.
  SUBROUTINE tendencies(t,y,res)
    REAL(KIND=8), INTENT(IN) :: t
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(IN) :: y
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(OUT) :: res
    CALL sparse_mul3(aotensor, y, y, res)
  END SUBROUTINE tendencies


  !-----------------------------------------------------!
  !                                                     !
  ! Adjoint model integrator                            !
  !                                                     !
  !-----------------------------------------------------!


  
  !> Routine to initialise the adjoint model integration buffers.
  SUBROUTINE init_ad_integrator
    INTEGER :: AllocStat
    ALLOCATE(ad_buf_y1(0:ndim),ad_buf_kA(0:ndim),ad_buf_kB(0:ndim), ad_buf_kAA(0:ndim),ad_buf_kBB(0:ndim) ,STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
  END SUBROUTINE init_ad_integrator

  !> Routine to perform an integration step (RK4 algorithm) of the adjoint model. The incremented time is returned.
  !> @param y Initial point.
  !> @param ystar Adjoint model at the point ystar.
  !> @param t Actual integration time
  !> @param dt Integration timestep.
  !> @param res Final point after the step.
  SUBROUTINE ad_step(y,ystar,t,dt,res)
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(IN) :: y,ystar
    REAL(KIND=8), INTENT(INOUT) :: t
    REAL(KIND=8), INTENT(IN) :: dt
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(OUT) :: res

    CALL ad(t,ystar,y,ad_buf_kA)
    ad_buf_y1 = y+0.5*dt*ad_buf_kA
    CALL ad(t+0.5*dt,ystar,ad_buf_y1,ad_buf_kB)
    ad_buf_y1 = y+0.5*dt*ad_buf_kB
    ad_buf_kA = ad_buf_kA+2*ad_buf_kB
    CALL ad(t+0.5*dt,ystar,ad_buf_y1,ad_buf_kB)
    ad_buf_y1 = y+0.5*dt*ad_buf_kB
    ad_buf_kA = ad_buf_kA+2*ad_buf_kB
    CALL ad(t+dt,ystar,ad_buf_y1,ad_buf_kB)
    ad_buf_kA = ad_buf_kA+ad_buf_kB
    res=y+ad_buf_kA*dt/6
    t=t+dt
  END SUBROUTINE ad_step

  !> Routine to perform a simultaneous integration step (RK4 algorithm) of the nonlinear and tangent linear model togheter. The incremented time is returned.
  !> @param y Model variable at time t
  !> @param deltay Perturbation at time t
  !> @param t Actual integration time
  !> @param dt Integration timestep.
  !> @param ynew Model variable at time t+dt
  !> @param deltaynew Perturbation at time t+dt
  SUBROUTINE evolve_ad_step(y,deltay,t,dt,ynew,deltaynew)
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(IN) :: y,deltay
    REAL(KIND=8), INTENT(INOUT) :: t
    REAL(KIND=8), INTENT(IN) :: dt
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(OUT) :: ynew,deltaynew

    CALL tendencies(t,y,ad_buf_kA)
    CALL ad(t,y,deltay,ad_buf_kAA)

    ad_buf_y1 = y + 0.5*dt*ad_buf_kA
    ad_buf_y11 = deltay + 0.5*dt*ad_buf_kAA

    CALL tendencies(t+0.5*dt,ad_buf_y1,ad_buf_kB)
    CALL ad(t+0.5*dt,ad_buf_y1,ad_buf_y11,ad_buf_kBB)
    
    ad_buf_y1 = y + 0.5*dt*ad_buf_kB
    ad_buf_y11 = deltay + 0.5*dt*ad_buf_kBB

    ad_buf_kA = ad_buf_kA + 2*ad_buf_kB
    ad_buf_kAA = ad_buf_kAA + 2*ad_buf_kBB
    
    CALL tendencies(t+0.5*dt,ad_buf_y1,ad_buf_kB)
    CALL ad(t+0.5*dt,ad_buf_y1,ad_buf_y11,ad_buf_kBB)
    
    ad_buf_y1 = y + dt*ad_buf_kB
    ad_buf_y11 = deltay + dt*ad_buf_kBB
    
    ad_buf_kA = ad_buf_kA + 2*ad_buf_kB
    ad_buf_kAA = ad_buf_kAA + 2*ad_buf_kBB
    
    CALL tendencies(t+dt,ad_buf_y1,ad_buf_kB)
    CALL ad(t+dt,ad_buf_y1,ad_buf_y11,ad_buf_kBB)

    ad_buf_kA = ad_buf_kA + ad_buf_kB
    ad_buf_kAA = ad_buf_kAA + ad_buf_kBB
    
    t=t+dt
    ynew=y+ad_buf_kA*dt/6
    deltaynew=deltay+ad_buf_kAA*dt/6
  END SUBROUTINE evolve_ad_step


  !-----------------------------------------------------!
  !                                                     !
  ! Tangent linear model integrator                     !
  !                                                     !
  !-----------------------------------------------------!

  !> Routine to initialise the tangent linear model integration buffers.
  SUBROUTINE init_tl_integrator
    INTEGER :: AllocStat
    ALLOCATE(tl_buf_y1(0:ndim),tl_buf_kA(0:ndim),tl_buf_kB(0:ndim), tl_buf_kAA(0:ndim),tl_buf_kBB(0:ndim) ,STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
  END SUBROUTINE init_tl_integrator

  !> Routine to perform an integration step (RK4 algorithm) of the tangent linear model. The incremented time is returned.
  !> @param y Initial point.
  !> @param ystar Adjoint model at the point ystar.
  !> @param t Actual integration time
  !> @param dt Integration timestep.
  !> @param res Final point after the step.
  SUBROUTINE tl_step(y,ystar,t,dt,res)
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(IN) :: y,ystar
    REAL(KIND=8), INTENT(INOUT) :: t
    REAL(KIND=8), INTENT(IN) :: dt
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(OUT) :: res

    CALL tl(t,ystar,y,tl_buf_kA)
    tl_buf_y1 = y+0.5*dt*tl_buf_kA
    CALL tl(t+0.5*dt,ystar,tl_buf_y1,tl_buf_kB)
    tl_buf_y1 = y+0.5*dt*tl_buf_kB
    tl_buf_kA = tl_buf_kA+2*tl_buf_kB
    CALL tl(t+0.5*dt,ystar,tl_buf_y1,tl_buf_kB)
    tl_buf_y1 = y+0.5*dt*tl_buf_kB
    tl_buf_kA = tl_buf_kA+2*tl_buf_kB
    CALL tl(t+dt,ystar,tl_buf_y1,tl_buf_kB)
    tl_buf_kA = tl_buf_kA+tl_buf_kB
    res=y+tl_buf_kA*dt/6
    t=t+dt
  END SUBROUTINE tl_step

  !> Routine to perform a simultaneous integration step (RK4 algorithm) of the nonlinear and tangent linear model togheter. The incremented time is returned.
  !> @param y Model variable at time t
  !> @param deltay Perturbation at time t
  !> @param t Actual integration time
  !> @param dt Integration timestep.
  !> @param ynew Model variable at time t+dt
  !> @param deltaynew Perturbation at time t+dt
  SUBROUTINE evolve_tl_step(y,deltay,t,dt,ynew,deltaynew)
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(IN) :: y,deltay
    REAL(KIND=8), INTENT(INOUT) :: t
    REAL(KIND=8), INTENT(IN) :: dt
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(OUT) :: ynew,deltaynew

    CALL tendencies(t,y,tl_buf_kA)
    CALL tl(t,y,deltay,tl_buf_kAA)

    tl_buf_y1 = y + 0.5*dt*tl_buf_kA
    tl_buf_y11 = deltay + 0.5*dt*tl_buf_kAA

    CALL tendencies(t+0.5*dt,tl_buf_y1,tl_buf_kB)
    CALL tl(t+0.5*dt,tl_buf_y1,tl_buf_y11,tl_buf_kBB)
    
    tl_buf_y1 = y + 0.5*dt*tl_buf_kB
    tl_buf_y11 = deltay + 0.5*dt*tl_buf_kBB

    tl_buf_kA = tl_buf_kA + 2*tl_buf_kB
    tl_buf_kAA = tl_buf_kAA + 2*tl_buf_kBB
    
    CALL tendencies(t+0.5*dt,tl_buf_y1,tl_buf_kB)
    CALL tl(t+0.5*dt,tl_buf_y1,tl_buf_y11,tl_buf_kBB)
    
    tl_buf_y1 = y + dt*tl_buf_kB
    tl_buf_y11 = deltay + dt*tl_buf_kBB
    
    tl_buf_kA = tl_buf_kA + 2*tl_buf_kB
    tl_buf_kAA = tl_buf_kAA + 2*tl_buf_kBB
    
    CALL tendencies(t+dt,tl_buf_y1,tl_buf_kB)
    CALL tl(t+dt,tl_buf_y1,tl_buf_y11,tl_buf_kBB)

    tl_buf_kA = tl_buf_kA + tl_buf_kB
    tl_buf_kAA = tl_buf_kAA + tl_buf_kBB
    
    t=t+dt
    ynew=y+tl_buf_kA*dt/6
    deltaynew=deltay+tl_buf_kAA*dt/6
  END SUBROUTINE evolve_tl_step
    
END MODULE rk4_tl_ad_integrator
