
! icdelta_def.f90
!
!>  Module to load the perturbation initial condition.
!
!> @copyright                                                               
!> 2015 Lesley De Cruz & Jonathan Demaeyer.
!> See LICENSE.txt for license information.                                  
!
!---------------------------------------------------------------------------!

MODULE icdelta_def

  USE params, only: natm,noc,ndim
  USE util, only: str,rstr
  USE inprod_analytic, only:awavenum,owavenum
  IMPLICIT NONE

  PRIVATE

  LOGICAL :: exists !< Boolean to test for file existence.
  
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE, PUBLIC :: ICdelta !< Initial condition vector

  PUBLIC ::load_ICdelta

CONTAINS

  !> Subroutine to load the initial condition if ICdelta.nml exists.
  !> If it does not, then write ICdelta.nml with random initial condition.
  SUBROUTINE load_ICdelta
    INTEGER :: i,AllocStat
    REAL(KIND=8) :: x
    CHARACTER(len=20) :: fm
    CHARACTER(len=20) :: fm2

    NAMELIST /IClist/ ICdelta

    CALL srand(1254)

    fm(1:6)='(F3.1)'
    fm2(1:7)='(F10.7)'
   
    IF (ndim == 0) STOP "*** Number of dimensions is 0! ***"
    ALLOCATE(ICdelta(0:ndim), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    INQUIRE(FILE='./ICdelta.nml',EXIST=exists)

    IF (exists) THEN
       OPEN(8, file="ICdelta.nml", status='OLD', recl=80, delim='APOSTROPHE')
       READ(8,nml=IClist)
       CLOSE(8)
    ELSE
       OPEN(8, file="ICdelta.nml", status='NEW')
       WRITE(8,'(a)') "!------------------------------------------------------------------------------!"
       WRITE(8,'(a)') "! Namelist file :                                                              !"
       WRITE(8,'(a)') "! Initial condition.                                                           !"
       WRITE(8,'(a)') "!------------------------------------------------------------------------------!"
       WRITE(8,*) ""
       WRITE(8,'(a)') "&ICLIST"
       WRITE(8,*) " ! psi variables"
       DO i=1,natm
          x=dble(2*rand()-1)
          WRITE(8,*) " ICDELTA("//TRIM(str(i))//") = "//rstr(x,fm2)//"   ! typ= "&
               &//awavenum(i)%typ//", Nx= "//TRIM(rstr(awavenum(i)&
               &%Nx,fm))//", Ny= "//TRIM(rstr(awavenum(i)%Ny,fm))
          ICdelta(i)=x
       END DO
       WRITE(8,*) " ! theta variables"
       DO i=1,natm
          x=dble(2*rand()-1)
          WRITE(8,*) " ICDELTA("//TRIM(str(i+natm))//") = "//rstr(x,fm2)//"   ! typ= "&
               &//awavenum(i)%typ//", Nx= "//TRIM(rstr(awavenum(i)&
               &%Nx,fm))//", Ny= "//TRIM(rstr(awavenum(i)%Ny,fm))
          ICdelta(i+natm)=x
       END DO

       WRITE(8,*) " ! A variables"
       DO i=1,noc
          x=dble(2*rand()-1)
          WRITE(8,*) " ICDELTA("//TRIM(str(i+2*natm))//") = "//rstr(x,fm2)//"   ! Nx&
               &= "//TRIM(rstr(owavenum(i)%Nx,fm))//", Ny= "&
               &//TRIM(rstr(owavenum(i)%Ny,fm))
          ICdelta(i+2*natm)=x
       END DO
       WRITE(8,*) " ! T variables"
       DO i=1,noc
          x=dble(2*rand()-1)
          WRITE(8,*) " ICDELTA("//TRIM(str(i+noc+2*natm))//") = "//rstr(x,fm2)//"   &
               &! Nx= "//TRIM(rstr(owavenum(i)%Nx,fm))//", Ny= "&
               &//TRIM(rstr(owavenum(i)%Ny,fm))
          ICdelta(i+noc+2*natm)=x
       END DO

       WRITE(8,'(a)') "&END"
       WRITE(8,*) ""
       CLOSE(8)
       WRITE(6,*) "*** ICdelta.nml namelist written. Starting with 0 as initial condition !***"
    ENDIF
    ICdelta(0)=1.0D0
  END SUBROUTINE load_ICdelta
END MODULE icdelta_def
