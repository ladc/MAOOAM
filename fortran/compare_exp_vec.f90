!
! This program is supposed to compare LEs and LVs of different runs.
! First compare.sh has to be executed. 
! This program reads the same namelists as the main program
! It creates files give the correlation of the vectors and the distance of the
! lyapunov exponents.
! 
!


PROGRAM compare_exp_vec
USE params, only: ndim, dt, tw, t_trans, t_run, writeout, rescaling_time, compute_BLV_LE,&
& compute_BLV, conv_BLV,compute_FLV, compute_FLV_LE, conv_FLV,compute_CLV, compute_CLV_LE,length_lyap,offset, &
& init_params
USE util, only: str
implicit none
CHARACTER(LEN=:) , ALLOCATABLE :: affix
!INTEGER(KIND=8), DIMENSION(:), ALLOCATABLE :: BLV_exp_bytesize, CLV_exp_bytesize, FLV_exp_bytesize, BLV_vec_bytesize, CLV_vec_bytesize, FLV_vec_bytesize
!INTEGER(KIND=8), DIMENSION(:), ALLOCATABLE :: BLV_exp_recordsize, CLV_exp_recordsize, FLV_exp_recordsize, BLV_vec_recordsize, CLV_vec_recordsize, FLV_vec_recordsize
INTEGER, DIMENSION(:), ALLOCATABLE :: BLV_exp_nrec, CLV_exp_nrec, FLV_exp_nrec, BLV_vec_nrec, CLV_vec_nrec, FLV_vec_nrec
REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: corrBLV, corrFLV, corrCLV,diffBLE, diffFLE, diffCLE
REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: BLE0,FLE0,CLE0 
REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: BLE1,FLE1,CLE1 
REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: BLV0,FLV0,CLV0 
REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: BLV1,FLV1,CLV1 

INTEGER :: k
INTEGER :: num_runs=2
INTEGER :: record
INTEGER :: nrecords
INTEGER, parameter :: diagout = 7
! Diagnosis outputfile
OPEN(UNIT=diagout, FILE='diag_compare.log')
! This introduces some flexibility for nameing the files
affix="_"
! number of runs
num_runs = 2


! Allocation part

ALLOCATE(BLV_exp_nrec(0:num_runs-1), CLV_exp_nrec(0:num_runs-1), FLV_exp_nrec(0:num_runs-1), BLV_vec_nrec(0:num_runs-1), CLV_vec_nrec(0:num_runs-1), FLV_vec_nrec(0:num_runs-1))

! Read parameters & namelists

CALL init_params

! number of records

nrecords = floor(dble(length_lyap+offset)/dble(rescaling_time))

! 

DO k=0,num_runs-1
  IF (compute_BLV) THEN
    IF (k.eq.0) ALLOCATE(BLV0(ndim,ndim),BLV1(ndim,ndim,num_runs-1),corrBLV(ndim))
    OPEN(UNIT=100*k+12,FILE='BLV_vec'//affix//trim(str(k))//'.dat', STATUS='OLD',ACCESS='DIRECT',RECL=8*ndim**2)
    BLV_vec_nrec(k)=nrecords
    WRITE(diagout,'(i2,a,i2,a,i15)') k,' BLV_vec_nrec(',k,'): ', BLV_vec_nrec(k)
  END IF
  IF (compute_BLV_LE) THEN
    IF (k.eq.0) ALLOCATE(BLE0(ndim),BLE1(ndim,num_runs-1),diffBLE(ndim))
    OPEN(UNIT=100*k+11,FILE='BLV_exp'//affix//trim(str(k))//'.dat', STATUS='OLD',ACCESS='DIRECT',RECL=8*ndim)
    BLV_vec_nrec(k)=nrecords
    BLV_exp_nrec(k)=nrecords
    WRITE(diagout,'(i2,a,i2,a,i15)') k,' BLV_exp_nrec(',k,'): ', BLV_exp_nrec(k)
  END IF
  IF (compute_FLV) THEN
    IF (k.eq.0) ALLOCATE(FLV0(ndim,ndim),FLV1(ndim,ndim,num_runs-1),corrFLV(ndim))
    OPEN(UNIT=100*k+22,FILE='FLV_vec'//affix//trim(str(k))//'.dat', STATUS='OLD',ACCESS='DIRECT',RECL=8*ndim**2)
    BLV_vec_nrec(k)=nrecords
    FLV_vec_nrec(k)=nrecords
    WRITE(diagout,'(i2,a,i2,a,i15)') k,' FLV_vec_nrec(',k,'): ', FLV_vec_nrec(k)
  END IF
  IF (compute_FLV_LE) THEN
    IF (k.eq.0) ALLOCATE(FLE0(ndim),FLE1(ndim,num_runs-1),diffFLE(ndim))
    OPEN(UNIT=100*k+21,FILE='FLV_exp'//affix//trim(str(k))//'.dat', STATUS='OLD',ACCESS='DIRECT',RECL=8*ndim)
    BLV_vec_nrec(k)=nrecords
    FLV_exp_nrec(k)=nrecords
    WRITE(diagout,'(i2,a,i2,a,i15)') k,' FLV_exp_nrec(',k,'): ', FLV_exp_nrec(k)
  END IF
  IF (compute_CLV) THEN
    IF (k.eq.0) ALLOCATE(CLV0(ndim,ndim),CLV1(ndim,ndim,num_runs-1),corrCLV(ndim))
    OPEN(UNIT=100*k+32,FILE='CLV_vec'//affix//trim(str(k))//'.dat', STATUS='OLD',ACCESS='DIRECT',RECL=8*ndim**2)
    BLV_vec_nrec(k)=nrecords
    CLV_vec_nrec(k)=nrecords
    WRITE(diagout,'(i2,a,i2,a,i15)') k,' CLV_vec_nrec(',k,'): ', CLV_vec_nrec(k)
  END IF
  IF (compute_CLV_LE) THEN
    IF (k.eq.0) ALLOCATE(CLE0(ndim),CLE1(ndim,num_runs-1),diffCLE(ndim))
    OPEN(UNIT=100*k+31,FILE='CLV_exp'//affix//trim(str(k))//'.dat', STATUS='OLD',ACCESS='DIRECT',RECL=8*ndim)
    BLV_vec_nrec(k)=nrecords
    CLV_exp_nrec(k)=nrecords
    WRITE(diagout,'(i2,a,i2,a,i15)') k,' CLV_exp_nrec(',k,'): ', CLV_exp_nrec(k)
  END IF
END DO

IF (COMPUTE_BLV_LE) THEN

  WRITE(*,*) '*** Start to cycle through BLE files ***'
  OPEN(UNIT=51, FILE='compare_diff_BLE.dat', ACCESS='DIRECT', RECL=ndim*8)
  DO record=1,BLV_exp_nrec(0) 
    READ(12,rec=record) BLE0(:)
    DO k=1,num_runs-1
      READ(100*k+12,rec=record) BLE1(:,k)
    END DO
    diffBLE(:)=0.0
    DO k=1,num_runs-1
      diffBLE(:) =diffBLE(:) + ABS(BLE0-BLE1(:,k))/dble(num_runs-1)
    END DO
    WRITE(51,rec=record) diffBLE(:)
  END DO

END IF

IF (COMPUTE_FLV_LE) THEN

  WRITE(*,*) '*** Start to cycle through FLE files ***'
  OPEN(UNIT=61, FILE='compare_diff_FLE.dat', ACCESS='DIRECT', RECL=ndim*8)
  DO record=1,FLV_exp_nrec(0) 
    READ(22,rec=record) FLE0(:)
    DO k=1,num_runs-1
      READ(100*k+22,rec=record) FLE1(:,k)
    END DO
    diffFLE(:)=0.0
    DO k=1,num_runs-1
      diffFLE(:) =diffFLE(:) + ABS(FLE0-FLE1(:,k))/dble(num_runs-1)
    END DO
    WRITE(61,rec=record) diffFLE(:)
  END DO

END IF

IF (COMPUTE_CLV_LE) THEN

  WRITE(*,*) '*** Start to cycle through CLE files ***'
  OPEN(UNIT=71, FILE='compare_diff_CLE.dat', ACCESS='DIRECT', RECL=ndim*8)
  DO record=1,CLV_exp_nrec(0) 
    READ(32,rec=record) CLE0(:)
    DO k=1,num_runs-1
      READ(100*k+32,rec=record) CLE1(:,k)
    END DO
    diffCLE(:)=0.0
    DO k=1,num_runs-1
      diffCLE(:) =diffCLE(:) + ABS(CLE0-CLE1(:,k))/dble(num_runs-1)
    END DO
    WRITE(71,rec=record) diffCLE(:)
  END DO

END IF

IF (COMPUTE_BLV) THEN

  WRITE(*,*) '*** Start to cycle through BLV files ***'
  OPEN(UNIT=51, FILE='compare_correlation_BLV.dat', ACCESS='DIRECT', RECL=ndim*8)
  DO record=1,BLV_vec_nrec(0) 
    READ(12,rec=record) BLV0(:,:)
    DO k=1,num_runs-1
      READ(100*k+12,rec=record) BLV1(:,:,k)
    END DO
    corrBLV(:)=0.0
    DO k=1,num_runs-1
      corrBLV(:) =corrBLV(:) + ABS(SUM(BLV0*BLV1(:,:,k),1))/dble(num_runs-1)
    END DO
    WRITE(51,rec=record) corrBLV(:)
  END DO

END IF

IF (COMPUTE_FLV) THEN

  WRITE(*,*) '*** Start to cycle through FLV files ***'
  OPEN(UNIT=61, FILE='compare_correlation_FLV.dat', ACCESS='DIRECT', RECL=ndim*8)
  DO record=1,FLV_vec_nrec(0) 
    READ(22,rec=record) FLV0(:,:)
    DO k=1,num_runs-1
      READ(100*k+22,rec=record) FLV1(:,:,k)
    END DO
    corrFLV(:)=0.0
    DO k=1,num_runs-1
      corrFLV(:) =corrFLV(:) + ABS(SUM(FLV0*FLV1(:,:,k),1))/dble(num_runs-1)
    END DO
    WRITE(61,rec=record) corrFLV(:)
  END DO

END IF

IF (COMPUTE_CLV) THEN

  WRITE(*,*) '*** Start to cycle through CLV files ***'
  OPEN(UNIT=71, FILE='compare_correlation_CLV.dat', ACCESS='DIRECT', RECL=ndim*8)
  DO record=1,CLV_vec_nrec(0) 
    READ(32,rec=record) CLV0(:,:)
    DO k=1,num_runs-1
      READ(100*k+32,rec=record) CLV1(:,:,k)
    END DO
    corrCLV(:)=0.0
    DO k=1,num_runs-1
      corrCLV(:) =corrCLV(:) + ABS(SUM(CLV0*CLV1(:,:,k),1))/dble(num_runs-1)
    END DO
    WRITE(71,rec=record) corrCLV(:)
  END DO

END IF
end program


