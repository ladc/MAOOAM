
! rk4_tl_ad_integrator.f90
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
!>  This module actually contains the RK4 algorithm routines.
!>  The user can modify it according to its preferred integration scheme.
!>  For higher-order schemes, additional bufers will probably have to be defined.
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
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: ad_buf_kBB !< Buffer to hold tendencies in the RK4 scheme for the adjoint model


  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: tl_buf_y1 !< Buffer to hold the intermediate position of the tangent linear model
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: tl_buf_y11 !< Buffer to hold the intermediate position of the tangent linear model
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: tl_buf_kA !< Buffer to hold tendencies in the RK4 scheme for the tangent linear model
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: tl_buf_kB !< Buffer to hold tendencies in the RK4 scheme for the tangent linear model
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: tl_buf_kC !< Buffer to hold tendencies in the RK4 scheme for the tangent linear model
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: tl_buf_kD !< Buffer to hold tendencies in the RK4 scheme for the tangent linear model
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: tl_buf_j1 !< Buffer to hold jacobians in the RK4 scheme for the tangent linear model
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: tl_buf_j2 !< Buffer to hold jacobians in the RK4 scheme for the tangent linear model
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: tl_buf_j3 !< Buffer to hold jacobians in the RK4 scheme for the tangent linear model
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: tl_buf_j4 !< Buffer to hold jacobians in the RK4 scheme for the tangent linear model
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: tl_buf_j1h !< Buffer to hold jacobians in the RK4 scheme for the tangent linear model
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: tl_buf_j2h !< Buffer to hold jacobians in the RK4 scheme for the tangent linear model
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: tl_buf_j3h !< Buffer to hold jacobians in the RK4 scheme for the tangent linear model
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: tl_buf_j4h !< Buffer to hold jacobians in the RK4 scheme for the tangent linear model
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: tl_buf_kAA !< Buffer to hold tendencies in the RK4 scheme for the tangent linear model
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: tl_buf_kBB !< Buffer to hold tendencies in the RK4 scheme for the tangent linear model

    
  PUBLIC :: init_ad_integrator, ad_step, evolve_ad_step, init_tl_integrator, tl_step, evolve_tl_step

CONTAINS

  !> Routine computing the tendencies of the nonlinear model
  !> @param t Time at which the tendencies have to be computed. Actually not needed for autonomous systems.
  !> @param y Point at which the tendencies have to be computed.
  !> @param res vector to store the result.
  !> @remark Note that it is NOT safe to pass `y` as a result bufer, 
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


  
  !> Routine to initialise the adjoint model integration bufers.
  SUBROUTINE init_ad_integrator
    INTEGER :: AllocStat,ii
    ALLOCATE(ad_buf_y1(0:ndim),ad_buf_kA(0:ndim),ad_buf_kB(0:ndim), ad_buf_kAA(0:ndim),ad_buf_kBB(0:ndim),tl_buf_j1h(ndim,ndim),tl_buf_j2h(ndim,ndim),tl_buf_j3h(ndim,ndim),tl_buf_j4h(ndim,ndim),tl_buf_j1(ndim,ndim),tl_buf_j2(ndim,ndim),tl_buf_j3(ndim,ndim),tl_buf_j4(ndim,ndim),one(ndim,ndim),STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
    
    one=0.0d0
    do ii=1,ndim
      one(ii,ii)=1.0d0
    enddo

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

  !> Routine to initialise the tangent linear model integration bufers.
  SUBROUTINE init_tl_integrator
    INTEGER :: AllocStat
    ALLOCATE(tl_buf_y1(0:ndim),tl_buf_kA(0:ndim),tl_buf_kB(0:ndim),tl_buf_kC(0:ndim),tl_buf_kD(0:ndim), tl_buf_kAA(0:ndim),tl_buf_kBB(0:ndim) ,STAT=AllocStat)
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

  !> Routine to perform a simultaneously an integration step (RK4 algorithm) of the nonlinear and computes the RK4 tangent linear propagator. The incremented time is returned.
  !> @param y Model variable at time t
  !> @param prop RK4 Propagator at time t
  !> @param t Actual integration time
  !> @param dt Integration timestep.
  !> @param ynew Model variable at time t+dt
  SUBROUTINE tl_prop_step(y,ensemble,t,dt,ynew)
    REAL(KIND=8), INTENT(INOUT) :: t
    REAL(KIND=8), INTENT(IN) :: dt
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(IN) :: y

    REAL(KIND=8), DIMENSION(ndim,ndim), INTENT(INOUT) :: ensemble
    REAL(KIND=8), DIMENSION(ndim,ndim) :: prop
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(OUT) :: ynew
 
    CALL tendencies(t,y,tl_buf_kA)
    tl_buf_j1=jacobian_matrix(y)

    tl_buf_y1 = y + 0.5*dt*tl_buf_kA

    CALL tendencies(t+0.5*dt,tl_buf_y1,tl_buf_kB)
    tl_buf_j2=jacobian_matrix(tl_buf_y1)
    
    tl_buf_y1 = y + 0.5*dt*tl_buf_kB

    CALL tendencies(t+0.5*dt,tl_buf_y2,tl_buf_kC)
    tl_buf_j3=jacobian_matrix(tl_buf_y1)

    tl_buf_y1 = y + dt*tl_buf_kC

    CALL tendencies(t+dt,tl_buf_y2,tl_buf_kD)
    tl_buf_j4=jacobian_matrix(tl_buf_y1)
    
    tl_buf_j1h=tl_buf_j1
    tl_buf_j2h=tl_buf_j2
    tl_buf_j3h=tl_buf_j3
    tl_buf_j4h=tl_buf_j4
    call dgemm ('n', 'n', ndim, ndim, ndim, dt/2., tl_buf_j2, ndim,tl_buf_j1h, ndim,1.0d0, tl_buf_j2h, ndim)
    call dgemm ('n', 'n', ndim, ndim, ndim, dt/2., tl_buf_j3, ndim,tl_buf_j2h, ndim,1.0d0, tl_buf_j3h, ndim)
    call dgemm ('n', 'n', ndim, ndim, ndim, dt , tl_buf_j4, ndim,tl_buf_j3h, ndim,1.0d0, tl_buf_j4h, ndim)
     
    ynew=y  + dt6*(tl_buf_kA + 2.*tl_buf_kB + 2.*tl_buf_kC + tl_buf_kD)
    prop=one+ dt6*(tl_buf_j4 + 2.*tl_buf_j2 + 2.*tl_buf_j3 + tl_buf_j4)
    tl_buf_j1h=ensemble
    call dgemm ('n', 'n', ndim, ndim, ndim, 1.0d0 , prop, ndim,tl_buf_j1h, ndim,0.0d0, ensemble, ndim)

    t=t+dt
   
  END SUBROUTINE tl_prop_step

END MODULE rk4_tl_ad_integrator
