! this sbr will find next line which a r-package can interact with
! now it is possible compute only lines with the delta function profile
SUBROUTINE next_line_red(approximation, pack_index, init_line, next1line, n_eqf_lines, tooRed)
USE types
IMPLICIT NONE

! input variables
INTEGER                                 :: init_line
INTEGER                                 :: approximation, pack_index
! loop variables
INTEGER                                 :: I
! line which we are now computing

! output variables
INTEGER                                 :: next1line, n_eqf_lines
LOGICAL                                 :: tooRed

INTEGER, PARAMETER                      :: test_pi = 0

IF(init_line == ntransitions .OR. init_line == ntransitions + 1) THEN
 next1line = ntransitions + 1
 tooRed = .TRUE.
 RETURN
END IF

IF(ntransitions == 0) THEN
 next1line = ntransitions + 1
 tooRed = .TRUE.
 RETURN
END IF

SELECT CASE(approximation)
! the Sobolev approximation: line profiles are delta functions
CASE(1)
! write(*,*) 'next_line: last_line = ', package(pack_index)%last_line
! write(*,*) 'next_line: n_eqf_lines = ', n_eqf_lines
IF(package(pack_index)%freq_cmf > linelist(ntransitions)%freq)THEN
 tooRed = .FALSE.
 IF (init_line == no_line) THEN
  ! In case we have more lines
  DO I = 1, ntransitions
   ! write(*,*) 'next_line: I = ', I, ' freq_cmf: ', package(pack_index)%freq_cmf, 'linelist:', linelist(I)%freq
   IF (package(pack_index)%freq_cmf > linelist(I)%freq) THEN
    ! write(*,*) 'next_line: new last_line = ', I - 1
    ! IF(I - 1 == 0) package(pack_index)%last_line = no_line
    next1line = I
    EXIT
    ! write(*,*)  'next_line: package(pack_index)%last_line = I-1', I-1
   END IF
  END DO 
 ELSE IF(init_line <= ntransitions .AND. init_line > 0) THEN! the last line /= no_line
  n_eqf_lines = 1
  DO I = init_line + 1, ntransitions
   IF(package(pack_index)%freq_cmf > linelist(I)%freq) THEN
    next1line = I
    EXIT
   END IF
  END DO
  !IF(package(pack_index)%freq_cmf > linelist(I)%freq) THEN
  ! write(*,*) 'next1line: lastline = ', lastline
  ! next1line = init_line + n_eqf_lines
 ELSE
  next1line = ntransitions + 1
  tooRed = .TRUE.
  RETURN
 END IF
ELSE ! f_cmf < f_reddest_line
 next1line = ntransitions + 1
 tooRed = .TRUE.
 RETURN
END IF

! write(*,*) 'next_line: next1line = ', next1line

linelist(next1line)%counted = .TRUE.

n_eqf_lines = 1
DO I = next1line + 1, ntransitions
 IF(linelist(I)%freq == linelist(next1line)%freq) THEN
  ! write(*,*) 'next1line: after next line: ', I
  linelist(I)%counted = .TRUE.
  n_eqf_lines = n_eqf_lines + 1
  CYCLE
 END IF
 EXIT
END DO
! package(pack_index)%last_line = next1line - 1
! write(*,*) 'next_line: next1line = ', next1line
 ! write(*,*) 'next_line: next1line = ', next1line, ' n_eqf_lines = ', n_eqf_lines
 ! write(*,*) 'next_line: f_cmf / f_line = ', package(pack_index)%freq_cmf / linelist(next1line)%freq
 ! write(*,*) 'next_line:(n-1) f_cmf / f_line = ', package(pack_index)%freq_cmf / linelist(next1line - 1)%freq
 ! write(*,*) 'next_line: f_cmf / f_line(n1 - 1) = ', package(pack_index)%freq_cmf / linelist(next1line - 1)%freq
 ! STOP 'next_line: testing'
! number of lines with the same frequency
CASE DEFAULT
 write(*,*) 'no default case'
END SELECT

! final checks
IF(next1line /= ntransitions + 1) THEN
 IF(package(pack_index)%freq_cmf < linelist(next1line)%freq) THEN
  write(*,*) 'next_line: init_line = ', init_line, ' next1line = ', next1line
  write(*,*) 'next1line: f_cmf / f_line = ', package(pack_index)%freq_cmf / linelist(next1line)%freq
  STOP 'next_line: f_cmf < f_line'
 END IF
END IF

IF(package(pack_index)%freq_cmf < linelist(ntransitions)%freq) THEN
 next1line = ntransitions + 1
 tooRed = .TRUE.
END IF
IF(package(pack_index)%freq_cmf <= linelist(next1line)%freq) THEN
 write(*,*) 'pack_index = ', pack_index
 write(*,*) ' f_cmf / f_line = ', package(pack_index)%freq_cmf / linelist(next1line)%freq
 STOP 'next_line'
 next1line = ntransitions + 1
 tooRed = .TRUE.
END IF

! IF(pack_index == test_pi) THEN
!  write(*,*) 'next_line: init_line = ', init_line, ' next1line = ', next1line
! END IF
! write(*,*) 'next_line: f_cmf / f_line = ', package(pack_index)%freq_cmf / linelist(next1line)%freq

IF(next1line > SIZE(linelist) + 1) THEN
 write(*,*) 'next_line: dim(linelist) = ', SIZE(linelist)
 write(*,*) 'next_line: next1line = ', next1line, ' n_eqf_lines = ', n_eqf_lines
 write(*,*) 'next_line: pack_index = ', pack_index
 STOP
END IF
END SUBROUTINE
