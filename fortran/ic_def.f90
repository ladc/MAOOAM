
! ic_def.f90
!
!>  Module to load the initial condition.
!
!> @copyright                                                               
!> 2015 Lesley De Cruz & Jonathan Demaeyer.
!> See LICENSE.txt for license information.                                  
!
!---------------------------------------------------------------------------!

MODULE ic_def

  USE params, only: natm,noc,ndim
  USE util, only: str,rstr
  USE inprod_analytic, only:awavenum,owavenum
  IMPLICIT NONE

  PRIVATE

  LOGICAL :: exists !< Boolean to test for file existence.
  
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE, PUBLIC :: IC !< Initial condition vector

  PUBLIC ::load_IC

CONTAINS

  !> Subroutine to load the initial condition if IC.nml exists.
  !> If it does not, then write IC.nml with 0 as initial condition.
  SUBROUTINE load_IC
    INTEGER :: i,AllocStat
    CHARACTER(len=20) :: fm
    REAL(KIND=8) :: size_of_random_noise
    CHARACTER(LEN=4) :: init_type 
    NAMELIST /IClist/ IC
    NAMELIST /RAND/ init_type,size_of_random_noise 



    fm(1:6)='(F3.1)'
   
    IF (ndim == 0) STOP "*** Number of dimensions is 0! ***"
    ALLOCATE(IC(0:ndim), STAT=AllocStat)
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"

    INQUIRE(FILE='./IC.nml',EXIST=exists)

    IF (exists) THEN
       OPEN(8, file="IC.nml", status='OLD', recl=80, delim='APOSTROPHE')
       READ(8,nml=IClist)
       READ(8,nml=RAND)
       SELECT CASE (init_type)
         CASE ('rand')
           CALL init_random_seed()
           CALL random_number(IC)
           IC=IC*size_of_random_noise*10.D0
           IC(0)=1.0d0
           WRITE(6,*) "*** IC.nml namelist written. Starting with random initial condition !***"
         CASE ('zero')
           IC=0
           IC(0)=1.0d0
           WRITE(6,*) "*** IC.nml namelist written. Starting with initial condition in IC.nml !***"
         CASE ('read')
           !nothing has to be done IC has already the right values
           WRITE(6,*) "*** IC.nml namelist written. Starting with initial condition in IC.nml !***"
       END SELECT
       CLOSE(8)
    ELSE
       IC=0
       IC(0)=1.0D0
       init_type="zero"
       size_of_random_noise=0.D0
       WRITE(6,*) "*** IC.nml namelist written. Starting with 0 as initial condition !***"

    END IF
    OPEN(8, file="IC.nml", status='REPLACE')
    WRITE(8,'(a)') "!------------------------------------------------------------------------------!"
    WRITE(8,'(a)') "! Namelist file :                                                              !"
    WRITE(8,'(a)') "! Initial condition.                                                           !"
    WRITE(8,'(a)') "!------------------------------------------------------------------------------!"
    WRITE(8,*) ""
    WRITE(8,'(a)') "&ICLIST"
    WRITE(8,*) " ! psi variables"
    DO i=1,natm
       WRITE(8,*) " IC("//TRIM(str(i))//") = ",IC(i+natm),"   ! typ= "&
            &//awavenum(i)%typ//", Nx= "//TRIM(rstr(awavenum(i)&
            &%Nx,fm))//", Ny= "//TRIM(rstr(awavenum(i)%Ny,fm))
    END DO
    WRITE(8,*) " ! theta variables"
    DO i=1,natm
       WRITE(8,*) " IC("//TRIM(str(i+natm))//") = ",IC(i+natm),"   ! typ= "&
            &//awavenum(i)%typ//", Nx= "//TRIM(rstr(awavenum(i)&
            &%Nx,fm))//", Ny= "//TRIM(rstr(awavenum(i)%Ny,fm))
    END DO

    WRITE(8,*) " ! A variables"
    DO i=1,noc
       WRITE(8,*) " IC("//TRIM(str(i+2*natm))//") = ",IC(i+2*natm),"   ! Nx&
            &= "//TRIM(rstr(owavenum(i)%Nx,fm))//", Ny= "&
            &//TRIM(rstr(owavenum(i)%Ny,fm))
    END DO
    WRITE(8,*) " ! T variables"
    DO i=1,noc
       WRITE(8,*) " IC("//TRIM(str(i+noc+2*natm))//") = ",IC(i+2*natm+noc),"   &
            &! Nx= "//TRIM(rstr(owavenum(i)%Nx,fm))//", Ny= "&
            &//TRIM(rstr(owavenum(i)%Ny,fm))
    END DO

    WRITE(8,'(a)') "&END"
    WRITE(8,*) ""
    WRITE(8,'(a)') "!------------------------------------------------------------------------------!"
    WRITE(8,'(a)') "! Namelist file :                                                              !"
    WRITE(8,'(a)') "! Initialisation type.                                                         !"
    WRITE(8,'(a)') "!------------------------------------------------------------------------------!"
    WRITE(8,'(a)') "! type = 'read': use IC; 'rand': random state; 'zero': zero condition "
    WRITE(8,*) ""
    WRITE(8,'(a)') "&RAND"
    WRITE(8,'(a)') "  init_type= '"//init_type//"'" 
    WRITE(8,'(a,d)') "  size_of_random_noise = ",size_of_random_noise
    WRITE(8,'(a)') "&END"
    WRITE(8,*) ""
    CLOSE(8)
    
  END SUBROUTINE load_IC
  
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
END MODULE ic_def
