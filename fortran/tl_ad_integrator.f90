
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



MODULE tl_ad_integrator

  USE params, only: ndim
  USE maooam_tl_ad, only: ad,tl
  IMPLICIT NONE

  PRIVATE

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: ad_buf_y1 !< Buffer to hold the intermediate position (Heun algorithm) of the adjoint model
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: ad_buf_f0 !< Buffer to hold tendencies at the initial position of the adjoint model
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: ad_buf_f1 !< Buffer to hold tendencies at the intermediate position of the adjoint model

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: tl_buf_y1 !< Buffer to hold the intermediate position (Heun algorithm) of the tangent linear model
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: tl_buf_f0 !< Buffer to hold tendencies at the initial position of the tangent linear model
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: tl_buf_f1 !< Buffer to hold tendencies at the intermediate position of the tangent linear model

    
  PUBLIC :: init_ad_integrator, ad_step, init_tl_integrator, tl_step, tl_prop_step

CONTAINS

  !-----------------------------------------------------!
  !                                                     !
  ! Adjoint model integrator                            !
  !                                                     !
  !-----------------------------------------------------!
  
  !> Routine to initialise the adjoint model integration buffers.
  SUBROUTINE init_ad_integrator
    INTEGER :: AllocStat
    ALLOCATE(ad_buf_y1(0:ndim),ad_buf_f0(0:ndim),ad_buf_f1(0:ndim) ,STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
  END SUBROUTINE init_ad_integrator

  !> Routine to perform an integration step (Heun algorithm) of the adjoint model. The incremented time is returned.
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
    
    CALL ad(t,ystar,y,ad_buf_f0)
    ad_buf_y1 = y+dt*ad_buf_f0
    CALL ad(t,ystar,ad_buf_y1,ad_buf_f1)
    res=y+0.5*(ad_buf_f0+ad_buf_f1)*dt
    t=t+dt
  END SUBROUTINE ad_step

  !-----------------------------------------------------!
  !                                                     !
  ! Tangent linear model integrator                     !
  !                                                     !
  !-----------------------------------------------------!

  !> Routine to initialise the tangent linear model integration buffers.
  SUBROUTINE init_tl_integrator
    INTEGER :: AllocStat
    ALLOCATE(tl_buf_y1(0:ndim),tl_buf_f0(0:ndim),tl_buf_f1(0:ndim) ,STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
  END SUBROUTINE init_tl_integrator

  !> Routine to perform an integration step (Heun algorithm) of the tangent linear model. The incremented time is returned.
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

    CALL tl(t,ystar,y,tl_buf_f0)
    tl_buf_y1 = y+dt*tl_buf_f0
    CALL tl(t,ystar,tl_buf_y1,tl_buf_f1)
    res=y+0.5*(tl_buf_f0+tl_buf_f1)*dt
    t=t+dt
  END SUBROUTINE tl_step

  !> Routine to perform a simultaneously an integration step (heun algorithm) of the nonlinear and computes the RK4 tangent linear propagator. The boolean variable adjoint allows for an adjoint forward integration. The incremented time is returned.
  !> @param y Model variable at time t
  !> @param prop heun Propagator at time t
  !> @param t Actual integration time
  !> @param dt Integration timestep.
  !> @param ynew Model variable at time t+dt
  SUBROUTINE tl_prop_step(y,propagator,t,dt,ynew,adjoint)
    REAL(KIND=8), INTENT(INOUT) :: t
    REAL(KIND=8), INTENT(IN) :: dt
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(IN) :: y
    LOGICAL, INTENT(IN) :: adjoint
    REAL(KIND=8), DIMENSION(ndim,ndim), INTENT(OUT) :: propagator
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(OUT) :: ynew
 
    CALL tendencies(t,y,tl_buf_kA)
    tl_buf_j1=jacobian_mat(y)

    tl_buf_y1 = y + dt*tl_buf_kA

    CALL tendencies(t+dt,tl_buf_y1,tl_buf_kB)
    tl_buf_j2=jacobian_mat(tl_buf_y1)
    
    tl_buf_j1h=tl_buf_j1
    tl_buf_j2h=tl_buf_j2
    call dgemm ('n', 'n', ndim, ndim, ndim, dt, tl_buf_j2, ndim,tl_buf_j1h, ndim,1.0d0, tl_buf_j2h, ndim)
     
    ynew=y  + dt/2.0d0*(tl_buf_kA + tl_buf_kB)
    IF (adjoint) THEN
            propagator=one - dt/2.0d0*(tl_buf_j1h + tl_buf_j2h)
    ELSE
            propagator=one + dt/2.0d0*(tl_buf_j1h + tl_buf_j2h)
    END IF        
    t=t+dt
   
  END SUBROUTINE tl_prop_step


  
END MODULE tl_ad_integrator
