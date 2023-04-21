SUBROUTINE event_dist2(pack_index, cell_dist, e_dist, event, actirrates)

 USE types
 USE rates_r

 IMPLICIT NONE    

 INTEGER                           :: I, pack_index, event, get_package_model_index
 INTEGER                           :: nextLine, current_mgi
 LOGICAL                           :: do_loop
!pointer to a field of continuum rates
 DOUBLE PRECISION                  :: e_dist, ran_numb, tau_rand, cell_dist, D, cont_dist
 DOUBLE PRECISION                  :: tau, l_dist, tau_line, tau_cont
 DOUBLE PRECISION                  :: electron_density, kappa_cont, dist
 DOUBLE PRECISION, PARAMETER       :: largeNumber = 1.D20
 ! number of lines with the same frequencies
 INTEGER                           :: n_next_lines
 INTEGER                           :: dummypackage
 INTEGER                           :: n_pack_d
 ! INTEGER                           :: OMP_GET_THREAD_NUM
 DOUBLE PRECISION                  :: freq_line
 TYPE(rrates)                      :: actirrates
 DOUBLE PRECISION                       :: ran2
 ! looking for next line
 INTEGER                           :: act_line
 DOUBLE PRECISION                  :: summ, tot_lop
 LOGICAL                                :: procout=.FALSE.
 LOGICAL                                :: inCell, tooRed
 LOGICAL                                :: raninit
 INTEGER                           :: lastLine
 INTEGER                           :: nloop

 ! temporary variables
 DOUBLE  PRECISION                      :: fr_line, loc_dist
 LOGICAL                                :: linCell
 INTEGER                                :: n_lines
      
 n_pack_d = SIZE(package)
 dummypackage = SIZE(package)
 ! write(*,*) 'event_dist: last_line = ', package(pack_index)%last_line

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
 do_loop = .TRUE.

 !Get the packet's current position on the model grid
 current_mgi = get_package_model_index(pack_index)

 electron_density = model_grid(current_mgi)%e_dens
 IF(current_mgi .EQ. n_modelgrid + 2) electron_density = 0.D0
 !print*, 'electron_density = ', electron_density

 ! calculates all continuum opacities
 CALL r_kappa_cont(pack_index, kappa_cont, actirrates)
 ! kappa_cont = 0.D0

 ! This is the opacity in co-moving frame. Must be transformed to the lab frame
 ! According to Mihalas and Mihalas Eq. 90.8 this is achieved by 
 CALL doppler_factor(pack_index, D)
 kappa_cont = D * kappa_cont


 ! initialization of n_next_lines to be equal to one
 n_next_lines = 1
 lastLine = package(pack_index)%last_line

nloop = 0
DO WHILE (do_loop) 

 nloop = nloop + 1
 IF(nloop == 1000) THEN
 !  write(*,*) 'event_dist: packet = ', pack_index
  STOP 'nloop == 100'
 END IF

 IF(kappa_cont /= 0) THEN
  cont_dist = (tau_rand - tau) / kappa_cont
 ELSE
  cont_dist = largeNumber
 END IF
 tau_cont = kappa_cont * cont_dist

 IF(nloop > 1) lastLine = nextLine + n_next_lines - 1
 ! write(*,*) 'event_dist: pack_index = ', pack_index, ' n_next_lines = ', n_next_lines
 ! write(*,*) 'event_dist: nloop = ', nloop, ' lastLine = ', lastLine, ' ntransitions = ', ntransitions
 
 IF (package(pack_index)%freq_cmf .GT. freq_line) THEN
  ! next_line: the first line the packet can interact with
  ! n_next_lines: the number of lines with the same frequency as the nextLine
  CALL next_line(1, pack_index, lastLine, nextLine, n_next_lines, tooRed)
  ! write(*,*) 'event_dist: nextLine = ', nextLine, ' n_next_lines = ', n_next_lines
  ALLOCATE(actirrates%Lline(n_next_lines), actirrates%nline(n_next_lines))
  IF(nextLine < ntransitions + 1) THEN
   freq_line = linelist(nextLine)%freq
   CALL resonance_distance(pack_index, nextLine, freq_line, cell_dist, l_dist, inCell, .TRUE.)
   CALL r_kappa_line(pack_index, current_mgi, nextLine, n_next_lines, actirrates, tau_line, l_dist)
  END IF
  
  ! write(*,*) 'event_dist2: nextLine = ', nextLine, ' tau_line = ', tau_line
  ! write(*,*) 'event_dist2: ntransitions = ', ntransitions, ' tooRed = ', tooRed
  if(l_dist < cont_dist .and. l_dist < cell_dist .and. .not. tooRed) then
   if(tau_cont > (tau_rand - tau).and. .not. tooRed) then
    if(tau_line + tau_cont > tau_rand - tau) then
     ! a line interraction occurs
     e_dist = dist + l_dist
     do_loop = .FALSE.
     ! nextLine: the line the packet will interact with
     CALL r_choose_line(pack_index, actirrates, n_next_lines, nextLine)
     package(pack_index)%last_line = nextLine
     event = rpkt_eventtype_lineinteraction
    else ! tau_line + tau_cont > tau
     ! choose next line
     ! write(*,*) 'event_dist2: '
     tau = tau + tau_line
    end if ! tau_line + tau_cont > tau
   end if ! tau_cont > tau
  else if (l_dist > cont_dist .and. cont_dist < cell_dist) THEN
   event = rpkt_eventtype_continuum
   e_dist = cont_dist
   do_loop = .FALSE.
  else ! l_dist < cont_dist & l_dist < cell_dist & not tooRed
   ! change of cell
   e_dist = dist + largeNumber
   event = rpkt_eventtype_changecell
   do_loop = .FALSE.
  end if
 ELSE
  IF (cont_dist < cell_dist) THEN
   event = rpkt_eventtype_continuum
   e_dist = cont_dist
   do_loop = .FALSE.
  ELSE ! cont_dist < cell_dist
   e_dist = dist + largeNumber
   event = rpkt_eventtype_changecell
   do_loop = .FALSE.
  END IF ! cont_dist < cell_dist
 END IF

 ! write(*,*) 'event_dist: inCell = ', inCell
 I = I + 1
 DEALLOCATE(actirrates%Lline, actirrates%nline)
 ! write(*,*) 'event_dist: eofloop, do_loop = ', do_loop
END DO
! STOP 'event_dist: testing'  
END SUBROUTINE event_dist2
