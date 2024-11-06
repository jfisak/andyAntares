SUBROUTINE r_choose_line(pack_index, actirrates, n_next_lines, n_chosenline)


USE types
USE constants
USE rates_r
IMPLICIT NONE


TYPE(rrates)                    :: actirrates
INTEGER                         :: n_next_lines, pack_index
INTEGER                         :: n_chosenline
DOUBLE PRECISION                :: ran2


INTEGER                         :: ind_I, act_line
DOUBLE PRECISION                :: tot_lop, summ, ran_numb

IF(n_next_lines == 1) THEN
 n_chosenline = actirrates%nline(1)
 RETURN
END IF

tot_lop = 0.D0
DO ind_I = 1, n_next_lines
 tot_lop = tot_lop + actirrates%Lline(ind_I)
END DO

ran_numb = ran2(idum) * tot_lop
summ = 0.D0
! write(*,*) 'r_choose_line: n_next_lines = ', n_next_lines
DO ind_I = 1, n_next_lines
 act_line = actirrates%nline(ind_I)
 ! write(*,*) 'r_choose_line: ind_I = ', I, ' ran_numb = ', ran_numb, ' summ = ', summ
 ! write(*,*) 'r_choose_line: ind_I = ', I, ' Lline = ', actirrates%Lline(ind_I)
 IF(ran_numb > summ .AND. ran_numb < summ + actirrates%Lline(ind_I)) THEN
  ! write(*,*) 'event_dist: last_line = ', act_line
  package(pack_index)%l_ele = linelist(act_line)%indexe
  package(pack_index)%l_ion = linelist(act_line)%indexi
  package(pack_index)%l_lev = linelist(act_line)%upper
  package(pack_index)%last_line = act_line
  ! write(*,*) 'event_dist: #1 chosen line = ', act_line
  EXIT
 END IF
 summ = summ + actirrates%Lline(ind_I)
END DO

IF(package(pack_index)%last_line <= 0) THEN
 write(*,*) 'r_choose_line: number of lines = ', actirrates%nline(:)
 write(*,*) 'r_choose_line: Lline = ', actirrates%Lline(:)
 write(*,*) 'r_choose_line: the chosen line = ', package(pack_index)%last_line, ' < 0'
 STOP 'r_choose_line'
END IF

END SUBROUTINE r_choose_line
