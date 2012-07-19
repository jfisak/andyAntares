SUBROUTINE read_input(n_pack, n_bin, iseed)

  USE types

  IMPLICIT NONE

  INTEGER             :: n_pack, n_bin, iseed, idx, npar
  CHARACTER  :: LINE*80, ACTPAR*20

  OPEN (UNIT=1, FILE='input.dat', STATUS='OLD')

! Default values
  n_pack = 10000
  n_bin = 100
  nx_cell = 100 ! number of cell in x direction
  ny_cell = 100 ! number of cell in y direction
  nz_cell = 100 ! number of cell in z direction
  n_modelgrid = 100
  xmax = 50.    ! coordinates of outer bourder of the wind in units of stellar radius 
  ymax = 50.  
  zmax = 50. 
  R_inf = 20.
  R_star = 2. 
  T_eff = 20000.
  M_dot = 1.D-6
  V_inf = 3000.D0
  iseed = -1  

  DO
    READ (1, '(A)', END=99) LINE

    CALL SARGV(LINE,1,ACTPAR)
    IF (ACTPAR .EQ. 'n_pack') THEN
    CALL SARGC (LINE, NPAR)
    IF (NPAR .LT. 2) GOTO 90
    CALL SARGV(LINE,2,ACTPAR)
    READ (ACTPAR, '(I20)', ERR=91) n_pack

    ELSE IF (ACTPAR .EQ. 'n_bin') THEN
    CALL SARGC (LINE, NPAR)
    IF (NPAR .LT. 2) GOTO 90
    CALL SARGV(LINE,2,ACTPAR)
    READ (ACTPAR, '(I20)', ERR=91) n_bin

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

    ELSE IF (ACTPAR .EQ. 'n_modelgrid') THEN
    CALL SARGC (LINE, NPAR)
    IF (NPAR .LT. 2) GOTO 90
    CALL SARGV(LINE,2,ACTPAR)
    READ (ACTPAR, '(I20)', ERR=91) n_modelgrid

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

    ELSE IF (ACTPAR .EQ. 'R_inf') THEN
    CALL SARGC (LINE, NPAR)
    IF (NPAR .LT. 2) GOTO 90
    CALL SARGV(LINE,2,ACTPAR)
    READ (ACTPAR, '(F20.0)', ERR=91) R_inf

    ELSE IF (ACTPAR .EQ. 'R_star') THEN
    CALL SARGC (LINE, NPAR)
    IF (NPAR .LT. 2) GOTO 90
    CALL SARGV(LINE,2,ACTPAR)
    READ (ACTPAR, '(F20.0)', ERR=91) R_star 

    ELSE IF (ACTPAR .EQ. 'T_eff') THEN
    CALL SARGC (LINE, NPAR)
    IF (NPAR .LT. 2) GOTO 90
    CALL SARGV(LINE,2,ACTPAR)
    READ (ACTPAR, '(F20.0)', ERR=91) T_eff 

    ELSE IF (ACTPAR .EQ. 'M_dot') THEN
    CALL SARGC (LINE, NPAR)
    IF (NPAR .LT. 2) GOTO 90
    CALL SARGV(LINE,2,ACTPAR)
    READ (ACTPAR, '(F20.0)', ERR=91) M_dot

    ELSE IF (ACTPAR .EQ. 'V_inf') THEN
    CALL SARGC (LINE, NPAR)
    IF (NPAR .LT. 2) GOTO 90
    CALL SARGV(LINE,2,ACTPAR)
    READ (ACTPAR, '(F20.0)', ERR=91) V_inf

    ELSE IF (ACTPAR .EQ. 'random_seed') THEN
    CALL SARGC (LINE, NPAR)
    IF (NPAR .LT. 2) GOTO 90
    CALL SARGV(LINE,2,ACTPAR)
    READ (ACTPAR, '(I20)', ERR=94) iseed
    ENDIF
  END DO

99 CONTINUE
  CLOSE (1)

  WRITE (*,'(A,I10)')   'npackages    = ', n_pack
  WRITE (*,'(A,I10)')   'nbin         = ', n_bin
  WRITE (*,'(A,I10)')   'nx_cell      = ', nx_cell
  WRITE (*,'(A,I10)')   'ny_cell      = ', ny_cell  
  WRITE (*,'(A,I10)')   'nz_cell      = ', nz_cell
  WRITE (*,'(A,I10)')   'n_modelgrid  = ', n_modelgrid
  WRITE (*,'(A,F10.3)') 'grid_xmax    = ', xmax
  WRITE (*,'(A,F10.3)') 'grid_ymax    = ', ymax
  WRITE (*,'(A,F10.3)') 'grid_zmax    = ', zmax
  WRITE (*,'(A,F10.3)') 'R_wind       = ', R_inf
  WRITE (*,'(A,F10.3)') 'R_star       = ', R_star ! R_star*r_sun
  WRITE (*,'(A,F10.3)') 'T_eff        = ', T_eff
  WRITE (*,'(A,F10.3)') 'M_dot        = ', M_dot 
  WRITE (*,'(A,F10.3)') 'V_inf        = ', V_inf
  IF (iseed .le. 0) then 
     WRITE (*,'(A)') 'Random-seed value is random '
  ELSE
     WRITE (*,'(A,I6)') 'random_seed = ', iseed
  ENDIF
  WRITE (*,'(3/)')

RETURN

!! Error branches
90 WRITE (*,*) '*** Not enough parameters'
   WRITE (*,*) '*** The error occured in the following line:'
   WRITE (*,*) LINE(:IDX(LINE))
   GOTO 100

91 WRITE (*,*) '*** Error when decoding floating point number'
   WRITE (*,*) '*** The error occured in the following line:'
   WRITE (*,*) LINE(:IDX(LINE))
   GOTO 100

92 WRITE (*,*) '*** Error: vel_dopp not given in input file'
   GOTO 100

93 WRITE (*,*) '*** Error: L0 not given in input file'
   GOTO 100

94 WRITE (*,*) '*** Error when decoding integer number'
   WRITE (*,*) '*** The error occured in the following line:'
   WRITE (*,*) LINE(:IDX(LINE))
   GOTO 100

100 STOP '**** FATAL ERROR in READ_INPUT'

END

