SUBROUTINE next_line_bluered(approximation, pack_index, b_dist, init_line, next_line, n_lines)

USE types
IMPLICIT NONE

INTEGER                                 :: approximation
INTEGER                                 :: pack_index
INTEGER                                 :: init_line
INTEGER                                 :: next_line, n_lines
DOUBLE PRECISION                        :: b_dist, f_cmf

INTEGER                                 :: I
INTEGER                                 :: cur_index
DOUBLE PRECISION                        :: cur_fline


LOGICAL                                 :: redshift, rob
LOGICAL                                 :: pack_redshift

DOUBLE PRECISION                        :: f_line, f_nextline

n_lines = 1
redshift = rob(pack_index, b_dist)

! write(*,*) 'next_line_bluered: redshift = ', redshift
! write(*,*) 'next_line_bluered: init_line = ', init_line

f_cmf = package(pack_index)%freq_cmf
pack_redshift = package(pack_index)%redshift

SELECT CASE(approximation)
 CASE(1)
 !!!!!!!!!!!!!!!!!
 ! if we have no initial line, we will have to go through the whole linelist
 ! write(*,*) 'next_line_bluered: init_line = ', init_line
 if(init_line == no_line) THEN
  IF(redshift) THEN
   IF(f_cmf > linelist(ntransitions)%freq) THEN
    DO I = 1, ntransitions
     cur_fline = linelist(I)%freq
     ! write(*,*) 'next_line_bluered: f_lu = ', cur_fline/f_cmf
     IF(f_cmf >  cur_fline) THEN
      next_line = I
      package(pack_index)%redshift = .true.
      EXIT
     END IF ! f_cmf > cur_fline
    END DO
   ELSE
    next_line = ntransitions + 1
   END IF ! cmf > f_lastline
  !!!!!!!!!!!!!!!!!
  ELSE ! blueshift
   IF(f_cmf < linelist(1)%freq) THEN
    DO I = 1, ntransitions
     cur_index = ntransitions - I + 1
     cur_fline = linelist(cur_index)%freq
     IF(f_cmf < cur_fline) THEN
      next_line = cur_index
      package(pack_index)%redshift = .false.
      EXIT
     END IF ! f_cmf > cur_fline
    END DO
   ELSE
    next_line = ntransitions + 1
   END IF
  END IF
 else ! init line is not no line
  IF(redshift) THEN
   if(pack_redshift) then
    next_line = init_line + 1
   else
    next_line = init_line
   end if
   package(pack_index)%redshift = .true.
  ELSE
   if(pack_redshift) then
    next_line = init_line
   else
    next_line = init_line - 1
   end if 
   package(pack_index)%redshift = .false.
  END IF ! is redshift
 end if ! init line
 if (next_line < 1 .or. next_line > ntransitions) then
  next_line = ntransitions + 1
 end if
 ! write(*,*) 'next_line_bluered: next_line = ', next_line
 ! number of lines with the same frequency
 if (next_line < ntransitions + 1) THEN
  f_nextline = linelist(next_line)%freq
  if(redshift) then
   do I = next_line + 1, ntransitions
    f_line = linelist(I)%freq
    if(f_nextline == f_line) then
     n_lines = n_lines + 1
    else ! line frequency differs
     EXIT
    end if
   end do
  else ! blueshift
   do I = next_line + 1, ntransitions
    cur_index = ntransitions - I + 1
    f_line = linelist(cur_index)%freq
    if(f_nextline == f_line) then
     n_lines = n_lines + 1
    else ! line frequency differs
     EXIT
    end if
   end do
  end if
 else
  next_line = ntransitions
  n_lines = 0
 end if ! next_line < ntransitions + 1
CASE DEFAULT
 write(*,*) 'next_line_bluered: the choice approximation = ', approximation, ' is not supported'
 STOP
END SELECT

! write(*,*) 'next_line_bluered: next_line = ', next_line
 
END SUBROUTINE next_line_bluered
