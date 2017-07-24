! this sbr will find next line which a r-package can interact with
! now it is possible compute only lines with the delta function profile
SUBROUTINE next_line(approximation, pack_index, next1line, n_eqf_lines)
USE types
IMPLICIT NONE

! input variables
INTEGER                                 :: approximation, pack_index
! loop variables
INTEGER                                 :: I
! the next first line
INTEGER                                 :: actLine
! line which we are now computing

! output variables
INTEGER                                 :: next1line, n_eqf_lines


SELECT CASE(approximation)
CASE(1)
IF (package(pack_index)%last_line .EQ. no_line) THEN
 ! In case we have more lines
 DO I = 1, ntransitions
    !print*, 'photon: ', pack_index, 'freq_cmf: ', package(pack_index)%freq_cmf, 'linelist:', linelist(I)%freq
  IF (package(pack_index)%freq_cmf .GT. linelist(I)%freq) THEN
   package(pack_index)%last_line = I - 1
   !print*, 'package(pack_index)%last_line = I-1', I-1
  END IF
 END DO 
 ! In case that package frequency can interact only with one more line from the line list,
 ! then index of the last line with which package interacted is ntransitions.
 ! We put (ntransitions - 1) only to be consistence with calculation of next_line, with which
 ! package may interact,should be general for any line interaction
 IF (package(pack_index)%last_line .EQ. no_line) package(pack_index)%last_line = ntransitions - 1
END IF


next1line = package(pack_index)%last_line + 1
! number of lines with the same frequency
! if there are not other lines this variable will be equal to one 
n_eqf_lines = 1
I = 1
DO
 actLine = next1line + I
 IF(actLine < ntransitions) THEN
  ! testing a frequency of next line: if it has the same frequency we will add this line to the "list"
  ! if the frequency differs, we will exit the loop because the rest of frequencies are totally not equal
  IF(linelist(actLine)%freq == linelist(actLine + 1)%freq) THEN
   n_eqf_lines = n_eqf_lines + 1
   I = I + 1
  ELSE
   EXIT
  END IF
 ELSE
  EXIT
 END IF
END DO

CASE DEFAULT
END SELECT
END SUBROUTINE
