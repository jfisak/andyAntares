SUBROUTINE read_input(n_pack, iseed)

  USE types
USE constants

  IMPLICIT NONE

  INTEGER    :: n_pack, iseed, idx, npar
  CHARACTER(LEN=180)  :: LINE, ACTPAR
  CHARACTER(LEN=180)  :: cur_calcmode
  CHARACTER(LEN=filename_length) :: linefile

  INTEGER                               :: calc_brtm_int

  OPEN (UNIT=1, FILE='input.dat', STATUS='OLD')

! Default values
  n_pack = 10000
  n_nubin = 100
  nx_cell = 100 ! number of cell in x direction
  ny_cell = 100 ! number of cell in y direction
  nz_cell = 100 ! number of cell in z direction
  model_type = 1
  iseed = -1  
  eldensfile = 0
  inputpopfile = ''
  inputmodelFile = ''
  inputComposition = ''
  abs_surface = 0
  saved_grid = 0
  enable_diffusion = 0
  sobolev_approximation = 1
  calc_brtm = .FALSE.

! 001 n_pack
! 002 n_nubin
! 003 nx_cell
! 004 ny_cell
! 005 nz_cell
! 006 model_type
! 007 xmax
! 008 ymax
! 009 zmax
! 010 iseed
! 011 inputflux
! 012 inputmodel
! 013 Nvirtpoint
! 014 dyngrid
! 015 nlte
! 016 velApprox
! 017 n_pack_save
! 018 abs_surface
! 019 inputmodelFile
! 020 inputcomposition
! 021 eldensfile
! 022 inputpopfile
! 023 saved grid
! 024 diffusion
! 025 BRTM
! 026 Sobolev approximation
! 027 Note
! 028 Calculation mode
! 029 
  DO
    READ (1, '(A)', END=99) LINE

    ! 023 saved grid
    ! 001 n_pack
    CALL SARGV(LINE,1,ACTPAR)
    IF (ACTPAR .EQ. 'n_pack') THEN
    CALL SARGC (LINE, NPAR)
    IF (NPAR .LT. 2) GOTO 90
    CALL SARGV(LINE,2,ACTPAR)
    READ (ACTPAR, '(I20)', ERR=91) n_pack

    ! 002 n_nubin
    ELSE IF (ACTPAR .EQ. 'n_nubin') THEN
    CALL SARGC (LINE, NPAR)
    IF (NPAR .LT. 2) GOTO 90
    CALL SARGV(LINE,2,ACTPAR)
    READ (ACTPAR, '(I20)', ERR=91) n_nubin

    ! 003 nx_cell
    ELSE IF (ACTPAR .EQ. 'nx_cell') THEN
    CALL SARGC (LINE, NPAR)
    IF (NPAR .LT. 2) GOTO 90
    CALL SARGV(LINE,2,ACTPAR)
    READ (ACTPAR, '(I20)', ERR=91) nx_cell

    ! 004 ny_cell
    ELSE IF (ACTPAR .EQ. 'ny_cell') THEN
    CALL SARGC (LINE, NPAR)
    IF (NPAR .LT. 2) GOTO 90
    CALL SARGV(LINE,2,ACTPAR)
    READ (ACTPAR, '(I20)', ERR=91) ny_cell

    ! 005 nz_cell
    ELSE IF (ACTPAR .EQ. 'nz_cell') THEN
    CALL SARGC (LINE, NPAR)
    IF (NPAR .LT. 2) GOTO 90
    CALL SARGV(LINE,2,ACTPAR)
    READ (ACTPAR, '(I20)', ERR=91) nz_cell

    ! 006 model_type
    ELSE IF (ACTPAR .EQ. 'model_type') THEN
    CALL SARGC (LINE, NPAR)
    IF (NPAR .LT. 2) GOTO 90
    CALL SARGV(LINE,2,ACTPAR)
    READ (ACTPAR, '(I20)', ERR=91) model_type

    ! 007 xmax
    ELSE IF (ACTPAR .EQ. 'xmax') THEN
    CALL SARGC (LINE, NPAR)
    IF (NPAR .LT. 2) GOTO 90
    CALL SARGV(LINE,2,ACTPAR)
    READ (ACTPAR, '(F20.0)', ERR=91) xmax

    ! 008 ymax
    ELSE IF (ACTPAR .EQ. 'ymax') THEN
    CALL SARGC (LINE, NPAR)
    IF (NPAR .LT. 2) GOTO 90
    CALL SARGV(LINE,2,ACTPAR)
    READ (ACTPAR, '(F20.0)', ERR=91) ymax

    ! 009 zmax
    ELSE IF (ACTPAR .EQ. 'zmax') THEN
    CALL SARGC (LINE, NPAR)
    IF (NPAR .LT. 2) GOTO 90
    CALL SARGV(LINE,2,ACTPAR)
    READ (ACTPAR, '(F20.0)', ERR=91) zmax

    ! 010 iseed
    ELSE IF (ACTPAR .EQ. 'random_seed') THEN
    CALL SARGC (LINE, NPAR)
    IF (NPAR .LT. 2) GOTO 90
    CALL SARGV(LINE,2,ACTPAR)
    READ (ACTPAR, '(I20)', ERR=94) iseed

    ! 011 inputflux
    ELSE IF (ACTPAR .EQ. 'inputflux') THEN
    CALL SARGC (LINE, NPAR)
    IF (NPAR .LT. 2) GOTO 90
    CALL SARGV(LINE,2,ACTPAR)
    READ (ACTPAR, '(I20)', ERR=94) inputflux 

    ! 012 inputmodel
    ELSE IF (ACTPAR .EQ. 'inputmodel') THEN
    CALL SARGC (LINE, NPAR)
    IF (NPAR .LT. 2) GOTO 90
    CALL SARGV(LINE,2,ACTPAR)
    READ (ACTPAR, '(I20)', ERR=94) inputmodel

    ! 013 Nvirtpoint
    ELSE IF (ACTPAR .EQ. 'Nvirtpoint') THEN
    CALL SARGC (LINE, NPAR)
    IF (NPAR .LT. 2) GOTO 90
    CALL SARGV(LINE,2,ACTPAR)
    READ (ACTPAR, '(I20)', ERR=94) Nvirtpoint

    ! 014 dyngrid
    ELSE IF (ACTPAR .EQ. 'dyngrid') THEN
    CALL SARGC (LINE, NPAR)
    IF (NPAR .LT. 2) GOTO 90
    CALL SARGV(LINE,2,ACTPAR)
    READ (ACTPAR, '(I20)', ERR=94) dyngrid

    ! 015 nlte
    ELSE IF (ACTPAR .EQ. 'nlte') THEN
    CALL SARGC (LINE, NPAR)
    IF (NPAR .LT. 2) GOTO 90
    CALL SARGV(LINE,2,ACTPAR)
    READ (ACTPAR, '(I20)', ERR=94) nlte

    ! 016 velApprox
    ELSE IF (ACTPAR .EQ. 'velApprox') THEN
    CALL SARGC (LINE, NPAR)
    IF (NPAR .LT. 2) GOTO 90
    CALL SARGV(LINE,2,ACTPAR)
    READ (ACTPAR, '(I20)', ERR=94) velApprox

    ! 017 n_pack_save
    ELSE IF (ACTPAR .EQ. 'n_pack_save') THEN
    CALL SARGC (LINE, NPAR)
    IF (NPAR .LT. 2) GOTO 90
    CALL SARGV(LINE,2,ACTPAR)
    READ (ACTPAR, '(I20)', ERR=94) n_pack_save

    ! 018 abs_surface
    ELSE IF (ACTPAR .EQ. 'abs_surface') THEN
    CALL SARGC (LINE, NPAR)
    IF (NPAR .LT. 2) GOTO 90
    CALL SARGV(LINE,2,ACTPAR)
    READ (ACTPAR, '(I20)', ERR=94) abs_surface

    ! 19 inputmodelFile
    ELSE IF (ACTPAR .EQ. 'inputmodelFile') THEN
    CALL SARGC (LINE, NPAR)
    IF (NPAR .LT. 2) GOTO 90
    CALL SARGV(LINE,2,ACTPAR)
    READ (ACTPAR, '(A60)', ERR=94) inputmodelFile

    ! 020 inputComposition
    ELSE IF (ACTPAR .EQ. 'inputComposition') THEN
    CALL SARGC (LINE, NPAR)
    IF (NPAR .LT. 2) GOTO 90
    CALL SARGV(LINE,2,ACTPAR)
    READ (ACTPAR, '(A)', ERR=94) inputcomposition

    ! 21 eldensfile
    ELSE IF (ACTPAR .EQ. 'eldensfile') THEN
    CALL SARGC (LINE, NPAR)
    IF (NPAR .LT. 2) GOTO 90
    CALL SARGV(LINE,2,ACTPAR)
    READ (ACTPAR, '(I20)', ERR=94) eldensfile

    ! 22 inputpopfile
    ELSE IF (ACTPAR .EQ. 'inputpopfile') THEN
    CALL SARGC (LINE, NPAR)
    IF (NPAR .LT. 2) GOTO 90
    CALL SARGV(LINE,2,ACTPAR)
    READ (ACTPAR, '(A20)', ERR=94) inputpopfile

    ! 23 saved grid
    ELSE IF (ACTPAR .EQ. 'saved_grid') THEN
    CALL SARGC (LINE, NPAR)
    IF (NPAR .LT. 2) GOTO 90
    CALL SARGV(LINE,2,ACTPAR)
    READ (ACTPAR, '(I2)', ERR=94) saved_grid

    ! 24 diffusion approximation
    ELSE IF (ACTPAR .EQ. 'diffusive') THEN
    CALL SARGC (LINE, NPAR)
    IF (NPAR .LT. 2) GOTO 90
    CALL SARGV(LINE,2,ACTPAR)
    READ (ACTPAR, '(I2)', ERR=94) enable_diffusion

    ! 025 BRTM
    ELSE IF (ACTPAR .EQ. 'calc_brtm_int') THEN
    CALL SARGC (LINE, NPAR)
    IF (NPAR .LT. 2) GOTO 90
    CALL SARGV(LINE,2,ACTPAR)
    READ (ACTPAR, '(I2)', ERR=94) calc_brtm_int
    IF(calc_brtm_int == 1) calc_brtm = .true.

    ! 26 Sobolev approximation
    ELSE IF (ACTPAR .EQ. 'sobolev') THEN
     CALL SARGC (LINE, NPAR)
     IF (NPAR .LT. 2) GOTO 90
     CALL SARGV(LINE,2,ACTPAR)
    READ (ACTPAR, '(I2)', ERR=94) sobolev_approximation
     ! 27 Note
     ! ELSE IF (ACTPAR .EQ. 'note') THEN
     ! CALL SARGC (LINE, NPAR)
     ! IF (NPAR .LT. 2) GOTO 90
     ! CALL SARGV(LINE,2,ACTPAR)
     ! DO
     ! write(*,*) 'read_input: note = ', put_a_note
     ! end all ifs
    ! 026 mode of calculation
    ELSE IF (ACTPAR .EQ. 'calcmode') THEN
     CALL SARGC (LINE, NPAR)
     IF (NPAR .LT. 2) GOTO 90
     CALL SARGV(LINE,2,ACTPAR)
     READ (ACTPAR, '(A)', ERR=94) cur_calcmode
     IF(cur_calcmode .EQ. 'oneline') oneline = .true.
     IF(cur_calcmode .EQ. 'multiline') oneline = .false.
     IF(cur_calcmode .EQ. 'brtm') only_brtm = .true.
    ! 027 a filename including data of the line
    ELSE IF (ACTPAR .EQ. 'linefile') THEN
     CALL SARGC (LINE, NPAR)
     IF (NPAR .LT. 2) GOTO 90
     CALL SARGV(LINE,2,ACTPAR)
     READ (ACTPAR, '(A)', ERR=94) linefile
     IF(oneline) THEN
      singleline_file = linefile
     END IF
    END IF
  END DO

  n_packets = n_pack

99 CONTINUE
  CLOSE (1)
  IF (iseed .LE. 0) THEN 
     WRITE (99,'(A)') 'Random-seed value is random '
  ELSE
     WRITE (99,'(A,I6)') 'random_seed = ', iseed
  ENDIF

RETURN

!! Error branches
90 WRITE (99,*) '*** Not enough parameters'
   WRITE (99,*) '*** The error occured in the following line:'
   WRITE (99,*) LINE(:IDX(LINE))
   GOTO 100

91 WRITE (99,*) '*** Error when decoding floating point number'
   WRITE (99,*) '*** The error occured in the following line:'
   WRITE (99,*) LINE(:IDX(LINE))
   GOTO 100

! 92 WRITE (99,*) '*** Error: vel_dopp not given in input file'
!    GOTO 100

! 93 WRITE (99,*) '*** Error: L0 not given in input file'
!    GOTO 100

94 WRITE (99,*) '*** Error when decoding integer number'
   WRITE (99,*) '*** The error occured in the following line:'
   WRITE (99,*) LINE(:IDX(LINE))
   GOTO 100

100 STOP '**** FATAL ERROR in READ_INPUT'

END

