SUBROUTINE r_choose_event(tau_rand, kappa_cont, tau_line, l_dist, cell_dist)

USE TYPES
USE constants
IMPLICIT NONE

LOGICAL                                 :: raninit
DOUBLE PRECISION                        :: ran_numb

raninit = .TRUE.
DO WHILE(raninit .EQV. .TRUE.)
 ran_numb = ran2(idum)  ! PUT IT IN SUBROUTINE - write is as do loop
 IF (ran_numb > 0.D0) THEN
  tau_rand = -LOG(ran_numb)
  ! write(*,*) 'event_dist: tau_rand = ', tau_rand
  raninit = .FALSE.
 END IF
END DO

! Initialize optical depth and distance
tau = 0.D0
dist = 0.D0

DO WHILE (do_loop == 1)
 ! calculation of tau_cont
 IF(lineint) THEN
  tau_cont = kappa_cont * l_dist
 ELSE
  tau_cont = kappa_cont * (cell_dist - dist)
  tline = 0.D0
 END IF
 tcont = tau_cont

 ! choosing an event
 IF(tau_rand + tau < tcont) THEN
  ! continuum event occurs
 ELSE IF(tau_rand >= tcont .AND. tau_rand < tline + tcont) THEN
  ! line event occurs
 ELSE
  ! choosing next line
 END IF
 

END DO



END SUBROUTINE r_choose_event
