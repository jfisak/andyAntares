! finds a next line in the blue or the red direction from the last line
! if there is more than one line with the same frequency, the n_lines > 1
!
! INPUT: approximation(INT) -- assumed approximation
!        pack_index(INT) -- index of a packet
!        b_dist(DBLE) -- a distance from a boundary
!        init_line(INT) -- initial line
! OUTPUT: next_line(INT) -- index of a next line
!         n_lines(INT) -- number of the lines with the same frequency as next_line
!
SUBROUTINE next_line_bluered(approximation, pack_index, b_dist, init_line, next_line, n_lines)

USE types
USE constants
IMPLICIT NONE

INTEGER                                 :: approximation
INTEGER                                 :: pack_index
INTEGER                                 :: init_line
INTEGER                                 :: next_line, n_lines
DOUBLE PRECISION                        :: b_dist, f_cmf

INTEGER                                 :: ind_I
INTEGER                                 :: cur_index
DOUBLE PRECISION                        :: cur_fline


LOGICAL                                 :: redshift, rob
LOGICAL                                 :: pack_redshift

DOUBLE PRECISION                        :: f_line, f_nextline
INTEGER, PARAMETER                      :: ind_init_line = 1

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
    DO ind_I = 1, ntransitions
     cur_fline = linelist(ind_I)%freq
     ! write(*,*) 'next_line_bluered: f_lu = ', cur_fline/f_cmf
     IF(f_cmf >  cur_fline) THEN
      next_line = ind_I
      package(pack_index)%redshift = .true.
      EXIT
     END IF ! f_cmf > cur_fline
    END DO
   ELSE
    next_line = ntransitions + 1
   END IF ! cmf > f_lastline
  !!!!!!!!!!!!!!!!!
  ELSE ! blueshift
   IF(f_cmf < linelist(ind_init_line)%freq) THEN
    DO ind_I = 1, ntransitions
     cur_index = ntransitions - ind_I + 1
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
   do ind_I = next_line + 1, ntransitions
    f_line = linelist(ind_I)%freq
    if(f_nextline == f_line) then
     n_lines = n_lines + 1
    else ! line frequency differs
     EXIT
    end if
   end do
  else ! blueshift
   do ind_I = next_line + 1, ntransitions
    cur_index = ntransitions - ind_I + 1
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

if(debug == 5) then
 write(*,*) 'next_line_bluered: f_cmf/f_next = ', package(pack_index)%freq_cmf/linelist(next_line)%freq
 if(next_line > 1) write(*,*) 'next_line_bluered: f_cmf/f_- = ', package(pack_index)%freq_cmf/linelist(next_line-1)%freq
 if(next_line < ntransitions) write(*,*) 'next_line_bluered: f_cmf/f_+ = ', package(pack_index)%freq_cmf/linelist(next_line+1)%freq
end if
! write(*,*) 'next_line_bluered: next_line = ', next_line
 
END SUBROUTINE next_line_bluered
