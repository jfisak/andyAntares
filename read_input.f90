SUBROUTINE read_input(n_pack, iseed)

  USE types

  IMPLICIT NONE

  INTEGER    :: n_pack, iseed, idx, npar
  CHARACTER  :: LINE*80, ACTPAR*20

  OPEN (UNIT=1, FILE='input.dat', STATUS='OLD')

! Default values
  n_pack = 10000
  n_nubin = 100
  nx_cell = 100 ! number of cell in x direction
  ny_cell = 100 ! number of cell in y direction
  nz_cell = 100 ! number of cell in z direction
  model_type = 1
  xmax = 50.    ! coordinates of outer bourder of the wind in units of stellar radius 
  ymax = 50.  
  zmax = 50. 
  iseed = -1  

  DO
    READ (1, '(A)', END=99) LINE

    CALL SARGV(LINE,1,ACTPAR)
    IF (ACTPAR .EQ. 'n_pack') THEN
    CALL SARGC (LINE, NPAR)
    IF (NPAR .LT. 2) GOTO 90
    CALL SARGV(LINE,2,ACTPAR)
    READ (ACTPAR, '(I20)', ERR=91) n_pack

    ELSE IF (ACTPAR .EQ. 'n_nubin') THEN
    CALL SARGC (LINE, NPAR)
    IF (NPAR .LT. 2) GOTO 90
    CALL SARGV(LINE,2,ACTPAR)
    READ (ACTPAR, '(I20)', ERR=91) n_nubin

    ELSE IF (ACTPAR .EQ. 'nx_cell') THEN
    CALL SARGC (LINE, NPAR)
    IF (NPAR .LT. 2) GOTO 90
    CALL SARGV(LINE,2,ACTPAR)
    READ (ACTPAR, '(I20)', ERR=91) nx_cell

    ELSE IF (ACTPAR .EQ. 'ny_cell') THEN
    CALL SARGC (LINE, NPAR)
    IF (NPAR .LT. 2) GOTO 90
    CALL SARGV(LINE,2,ACTPAR)
    READ (ACTPAR, '(I20)', ERR=91) ny_cell

    ELSE IF (ACTPAR .EQ. 'nz_cell') THEN
    CALL SARGC (LINE, NPAR)
    IF (NPAR .LT. 2) GOTO 90
    CALL SARGV(LINE,2,ACTPAR)
    READ (ACTPAR, '(I20)', ERR=91) nz_cell

    ELSE IF (ACTPAR .EQ. 'model_type') THEN
    CALL SARGC (LINE, NPAR)
    IF (NPAR .LT. 2) GOTO 90
    CALL SARGV(LINE,2,ACTPAR)
    READ (ACTPAR, '(I20)', ERR=91) model_type

    ELSE IF (ACTPAR .EQ. 'xmax') THEN
    CALL SARGC (LINE, NPAR)
    IF (NPAR .LT. 2) GOTO 90
    CALL SARGV(LINE,2,ACTPAR)
    READ (ACTPAR, '(F20.0)', ERR=91) xmax

    ELSE IF (ACTPAR .EQ. 'ymax') THEN
    CALL SARGC (LINE, NPAR)
    IF (NPAR .LT. 2) GOTO 90
    CALL SARGV(LINE,2,ACTPAR)
    READ (ACTPAR, '(F20.0)', ERR=91) ymax

    ELSE IF (ACTPAR .EQ. 'zmax') THEN
    CALL SARGC (LINE, NPAR)
    IF (NPAR .LT. 2) GOTO 90
    CALL SARGV(LINE,2,ACTPAR)
    READ (ACTPAR, '(F20.0)', ERR=91) zmax

    ELSE IF (ACTPAR .EQ. 'random_seed') THEN
    CALL SARGC (LINE, NPAR)
    IF (NPAR .LT. 2) GOTO 90
    CALL SARGV(LINE,2,ACTPAR)
    READ (ACTPAR, '(I20)', ERR=94) iseed

    ELSE IF (ACTPAR .EQ. 'inputflux') THEN
    CALL SARGC (LINE, NPAR)
    IF (NPAR .LT. 2) GOTO 90
    CALL SARGV(LINE,2,ACTPAR)
    READ (ACTPAR, '(I20)', ERR=94) inputflux 

    ELSE IF (ACTPAR .EQ. 'inputmodel') THEN
    CALL SARGC (LINE, NPAR)
    IF (NPAR .LT. 2) GOTO 90
    CALL SARGV(LINE,2,ACTPAR)
    READ (ACTPAR, '(I20)', ERR=94) inputmodel

    ELSE IF (ACTPAR .EQ. 'Nvirtpart') THEN
    CALL SARGC (LINE, NPAR)
    IF (NPAR .LT. 2) GOTO 90
    CALL SARGV(LINE,2,ACTPAR)
    READ (ACTPAR, '(I20)', ERR=94) Nvirtpart

    ELSE IF (ACTPAR .EQ. 'dyngrid') THEN
    CALL SARGC (LINE, NPAR)
    IF (NPAR .LT. 2) GOTO 90
    CALL SARGV(LINE,2,ACTPAR)
    READ (ACTPAR, '(I20)', ERR=94) dyngrid

    ELSE IF (ACTPAR .EQ. 'nlte') THEN
    CALL SARGC (LINE, NPAR)
    IF (NPAR .LT. 2) GOTO 90
    CALL SARGV(LINE,2,ACTPAR)
    READ (ACTPAR, '(I20)', ERR=94) nlte

    ELSE IF (ACTPAR .EQ. 'velApprox') THEN
    CALL SARGC (LINE, NPAR)
    IF (NPAR .LT. 2) GOTO 90
    CALL SARGV(LINE,2,ACTPAR)
    READ (ACTPAR, '(I20)', ERR=94) velApprox

    ELSE IF (ACTPAR .EQ. 'n_pack_save') THEN
    CALL SARGC (LINE, NPAR)
    IF (NPAR .LT. 2) GOTO 90
    CALL SARGV(LINE,2,ACTPAR)
    READ (ACTPAR, '(I20)', ERR=94) n_pack_save
    ENDIF
  END DO

99 CONTINUE
  CLOSE (1)

  WRITE (99,'(A,I10)')   'npackages        = ', n_pack
  WRITE (99,'(A,I10)')   'n_nubin          = ', n_nubin
  WRITE (99,'(A,I10)')   'nx_cell          = ', nx_cell
  WRITE (99,'(A,I10)')   'ny_cell          = ', ny_cell  
  WRITE (99,'(A,I10)')   'nz_cell          = ', nz_cell
  WRITE (99,'(A,I10)')   'model_type       = ', model_type
  WRITE (99,'(A,F10.3)') 'grid_xmax (Rsun) = ', xmax
  WRITE (99,'(A,F10.3)') 'grid_ymax (Rsun) = ', ymax
  WRITE (99,'(A,F10.3)') 'grid_zmax (Rsun) = ', zmax
  WRITE (99,'(A,I10)') 'inputflux = ', inputflux
  WRITE (99,'(A,I10)') 'inputmodel = ', inputmodel
  WRITE (99,'(A,I10)') 'dyngrid = ', dyngrid
  WRITE (99,'(A,I10)') 'nlte = ', nlte
  IF (iseed .LE. 0) THEN 
     WRITE (99,'(A)') 'Random-seed value is random '
  ELSE
     WRITE (99,'(A,I6)') 'random_seed = ', iseed
  ENDIF
  WRITE (99,'(3/)')

  ! Convert quantities to cgs 
!  xmax = xmax * r_sun
!  ymax = ymax * r_sun  
!  zmax = zmax * r_sun
!!
!
!  WRITE (99,'(A,I10)')   'npackages    = ', n_pack
!  WRITE (99,'(A,I10)')   'n_nubin      = ', n_nubin
!  WRITE (99,'(A,I10)')   'nx_cell      = ', nx_cell
!  WRITE (99,'(A,I10)')   'ny_cell      = ', ny_cell  
!  WRITE (99,'(A,I10)')   'nz_cell      = ', nz_cell
!  WRITE (99,'(A,I10)')   'model_type   = ', model_type
!  WRITE (99,'(A,G10.3)') 'grid_xmax    = ', xmax
!  WRITE (99,'(A,G10.3)') 'grid_ymax    = ', ymax
!  WRITE (99,'(A,G10.3)') 'grid_zmax    = ', zmax
!  IF (iseed .le. 0) then 
!     WRITE (99,'(A)') 'Random-seed value is random '
!  ELSE
!     WRITE (99,'(A,I6)') 'random_seed = ', iseed
!  ENDIF
!  WRITE (99,'(3/)')

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

92 WRITE (99,*) '*** Error: vel_dopp not given in input file'
   GOTO 100

93 WRITE (99,*) '*** Error: L0 not given in input file'
   GOTO 100

94 WRITE (99,*) '*** Error when decoding integer number'
   WRITE (99,*) '*** The error occured in the following line:'
   WRITE (99,*) LINE(:IDX(LINE))
   GOTO 100

100 STOP '**** FATAL ERROR in READ_INPUT'

END

