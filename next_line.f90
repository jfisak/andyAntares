! this sbr will find next line which a r-package can interact with
! now it is possible compute only lines with the delta function profile
SUBROUTINE next_line(approximation, pack_index, next1line, n_eqf_lines)
USE types
IMPLICIT NONE

! input variables
INTEGER                                 :: approximation, pack_index
! loop variables
INTEGER                                 :: I
! line which we are now computing

! output variables
INTEGER                                 :: next1line, n_eqf_lines
INTEGER                                 :: lastline

SELECT CASE(approximation)
! the Sobolev approximation
CASE(1)
! write(*,*) 'next_line: last_line = ', package(pack_index)%last_line
! print*, 'next_line: n_eqf_lines = ', n_eqf_lines
IF (package(pack_index)%last_line .EQ. no_line) THEN
 ! In case we have more lines
 DO I = 1, ntransitions
  ! write(*,*) 'next_line: I = ', I, ' freq_cmf: ', package(pack_index)%freq_cmf, 'linelist:', linelist(I)%freq
  IF (package(pack_index)%freq_cmf > linelist(I)%freq) THEN
   ! write(*,*) 'next_line: new last_line = ', I - 1
   package(pack_index)%last_line = I - 1
   IF(I - 1 == 0) package(pack_index)%last_line = no_line
   next1line = I
   EXIT
   ! write(*,*)  'next_line: package(pack_index)%last_line = I-1', I-1
  END IF
 END DO 
 ! In case that package frequency can interact only with one more line from the line list,
 ! then index of the last line with which package interacted is ntransitions.
 ! We put (ntransitions - 1) only to be consistence with calculation of next_line, with which
 ! package may interact,should be general for any line interaction
 IF (package(pack_index)%last_line .EQ. no_line) THEN
  package(pack_index)%last_line = ntransitions - 1
  next1line = package(pack_index)%last_line + 1
 END IF
ELSE ! the last line /= no_line
 lastline = package(pack_index)%last_line
 ! write(*,*) 'next1line: lastline = ', lastline
 IF(lastline == 0) lastline = 1
 DO I = lastline, ntransitions
  ! write(*,*) 'next_line: I = ', I, ' freq / f_line = ', package(pack_index)%freq_cmf / linelist(I)%freq
  IF(package(pack_index)%freq_cmf > linelist(I)%freq) THEN
   next1line = I
   ! write(*,*) 'next_line: next1line = ', I
   EXIT
  END IF
  IF(I == ntransitions) package(pack_index)%last_line = ntransitions
 END DO
END IF

! write(*,*) 'next_line: after next1line = ', next1line

IF(package(pack_index)%last_line == ntransitions) THEN
 ! write(*,*) 'next_line: last line == ntransitions'
 next1line = ntransitions
END IF

n_eqf_lines = 1
DO I = next1line + 1, ntransitions
 IF(linelist(I)%freq == linelist(next1line)%freq) THEN
  n_eqf_lines = n_eqf_lines + 1
  CYCLE
 END IF
 EXIT
END DO
 ! write(*,*) 'next_line: next1line = ', next1line, ' n_eqf_lines = ', n_eqf_lines
 ! write(*,*) 'next_line: f_cmf / f_line = ', package(pack_index)%freq_cmf / linelist(next1line)%freq
 ! write(*,*) 'next_line: f_cmf / f_line(n1 - 1) = ', package(pack_index)%freq_cmf / linelist(next1line - 1)%freq
! number of lines with the same frequency
CASE DEFAULT
END SELECT

! final checks
! IF(package(pack_index)%freq_cmf < linelist(next1line)%freq) THEN
!  STOP 'next_line: f_cmf < f_line'
! END IF
IF(next1line > SIZE(linelist)) THEN
 write(*,*) 'next_line: dim(linelist) = ', SIZE(linelist)
 write(*,*) 'next_line: next1line = ', next1line, ' n_eqf_lines = ', n_eqf_lines
 write(*,*) 'next_line: pack_index = ', pack_index
 STOP
END IF
END SUBROUTINE
