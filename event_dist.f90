SUBROUTINE event_dist(pack_index, cell_dist, e_dist, event, actirrates)

 USE types
USE constants
 USE rates_r

 IMPLICIT NONE    

 INTEGER                           :: I, pack_index, event, get_package_model_index
 INTEGER                           :: nextLine, current_mgi
 LOGICAL                           :: do_loop
!pointer to a field of continuum rates
 DOUBLE PRECISION                  :: e_dist, ran_numb, tau_rand, cell_dist, D
 DOUBLE PRECISION                  :: tau, l_dist, tau_line, tau_cont
 DOUBLE PRECISION                  :: electron_density, kappa_cont, dist
 DOUBLE PRECISION, PARAMETER       :: largeNumber = 1.D20
 ! number of lines with the same frequencies
 INTEGER                           :: n_next_lines
 INTEGER                           :: dummypackage
 INTEGER                           :: n_pack_d
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

 INTEGER                                :: cur_approx

IF(debug == 4 .or. debug == 5) procout = .TRUE.

      
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
 inCell=.FALSE.

 !Get the packet's current position on the model grid
 current_mgi = get_package_model_index(pack_index)
!  write(*,*) 'event_dist: pack_index = ', pack_index, ' current_mgi = ', current_mgi

 IF(current_mgi > n_modelgrid) THEN
  event = rpkt_eventtype_changecell
  e_dist = cell_dist + largeNumber
  RETURN
 END IF

 electron_density = model_grid(current_mgi)%e_dens

 ! calculates all continuum opacities
 ! write(*,*) 'event_dist: calling r_kappa_cont for a packet = ', pack_index
 CALL r_kappa_cont(pack_index, kappa_cont, actirrates)
 kappa_cont = 0.D0

 ! This is the opacity in co-moving frame. Must be transformed to the lab frame
 ! According to Mihalas and Mihalas Eq. 90.8 this is achieved by 
 CALL doppler_factor(pack_index, D)
 kappa_cont = D * kappa_cont
 ! write(*,*) 'event_dist: kappa_cont = ', kappa_cont

 ! initialization of n_next_lines to be equal to one
 lastLine = package(pack_index)%last_line

nloop = 0
DO WHILE (do_loop) 

 nloop = nloop + 1
 IF(nloop == 1000) THEN
  write(*,*) 'event_dist: nextLine = ', lastLine
  write(*,*) 'event_dist: packet = ', pack_index
  STOP 'nloop == 1000'
 END IF

 ! write(*,*) 'event_dist: pack_index = ', pack_index, ' n_next_lines = ', n_next_lines
 ! write(*,*) 'event_dist: nloop = ', nloop, ' lastLine = ', lastLine, ' ntransitions = ', ntransitions
 
 cur_approx = 1
 CALL next_line_bluered(cur_approx, pack_index, cell_dist, lastLine, nextLine, n_next_lines)
 ! write(*,*) 'event_dist: nextLine = ', nextLine, ' n_next_lines = ', n_next_lines
 ALLOCATE(actirrates%Lline(n_next_lines), actirrates%nline(n_next_lines))
 IF(nextLine < ntransitions + 1) THEN
  freq_line = linelist(nextLine)%freq
  CALL resonance_distance2(pack_index, nextLine, cell_dist, inCell, l_dist)
  IF(inCell) THEN
   CALL r_kappa_line(pack_index, current_mgi, nextLine, n_next_lines, l_dist, actirrates, tau_line)
  ELSE
   tau_line = 0.e0
  END IF
 ELSE
  tau_line = 0.e0
 END IF

  ! write(*,*) 'event_dist: pack_index = ', pack_index, ' l_dist = ', l_dist
  ! write(*,*) 'event_dist: l_dist/cell_dist = ', l_dist/cell_dist, ' nextLine = ', nextLine
  ! write(*,*) 'event_dist: tau_line = ', tau_line, ' tau_cont = ', tau_cont, ' tau_rand = ', tau_rand

 IF(inCell .AND. nextLine /= ntransitions + 1 .AND. .NOT. tooRed) THEN
 
 ! Calculate optical depth in the next line (Sobolev, dv/dr dependent)
 ! and continuum optical depth accumulated up to the line
 !print*, 'before moving package #', pack_index

 if(inCell) then
  tau_cont = kappa_cont * l_dist
 else
  tau_cont = kappa_cont * cell_dist
 end if

 
 if(procout) then
  write(*,*) 'event_dist: pack_index = ', pack_index, ' l_dist/cell_dist = ', l_dist/cell_dist, ' nextLine = ', nextLine
  write(*,*) 'event_dist: tau_line = ', tau_line, ' tau_cont = ', tau_cont, ' tau_rand = ', tau_rand
 end if
 
  IF(current_mgi .EQ. n_modelgrid + 2) tau_cont = 0.D0
  ! write(*,*) 'event_dist:', tau_line, tau_cont, tau_rand, tau
 
  ! Now do a step by step analysis of which event occurs and return the 
  ! distance and corresponding event
  IF ((tau_rand - tau) .GT. tau_cont) THEN
   ! if #02
   IF ((tau_rand - tau) .GT. (tau_cont + tau_line)) THEN
    dist = l_dist
    ! if #03
    IF (dist .GT. cell_dist) THEN
     ! In this case the package propagates to the next cell
     e_dist = cell_dist + largeNumber
     do_loop = .FALSE.
     event = rpkt_eventtype_changecell
     if(procout) write(*,*) 'event_dist: rpkt_eventtype_changecell'
    ! if #03
    ELSE
     ! choosing next line
     tau = tau + tau_cont + tau_line
     if(procout) write(*,*) 'event_dist: pack_index = ', pack_index, ' nextLine = ', nextLine
     if(package(pack_index)%redshift) then
      lastLine = nextLine + n_next_lines - 1
     else
      lastLine = nextLine - n_next_lines + 1
     end if
     if(procout) write(*,*) 'event_dist: choosing next line nextLine = ', nextLine
    ! if #03
    END IF
   ! if #02
   ELSE
    e_dist = l_dist
    do_loop = .FALSE.
    event = rpkt_eventtype_lineinteraction
    ! write(*,*) 1.D8 * light_speed / package(pack_index)%freq_cmf, &
    !  1.D8 * light_speed / package(pack_index)%freq_rf
    if(procout) write(*,*) 'event_dist: rpkt_eventtype_lineinteraction'
    ! choosing the line
    CALL r_choose_line(pack_index, actirrates, n_next_lines, nextLine)
    package(pack_index)%last_line = nextLine
   ! if #02
   END IF
  ELSE
   ! Continuum process will happen
   e_dist = dist + (tau_rand - tau) / kappa_cont
   do_loop = .FALSE. 
   event = rpkt_eventtype_continuum
   if(procout) write(*,*) 'event_dist: rpkt_eventtype_continuum 1'
 !   print*, 'cont.process happens',  tau_line, tau_cont
  END IF
 ELSE    
  ! The package cmf frequency is too red to interact to another
  ! line - No line interact anymore
  tau_cont = kappa_cont * (cell_dist - dist)
  ! IF(pack_index == 2222) write(*,*) 'event_dist: #2', tau_cont, tau_rand, tau
  IF ((tau_rand - tau) .GT. tau_cont) THEN
   e_dist = cell_dist + 1.D20
   do_loop = .FALSE.   
   event = rpkt_eventtype_changecell
   if(procout) write(*,*) 'event_dist: rpkt_eventtype_changecell'
  ELSE
   ! Continuum absorption happens
   e_dist = dist + (tau_rand - tau) / kappa_cont
   do_loop = .FALSE.                 
   event = rpkt_eventtype_continuum
   if(procout) write(*,*) 'event_dist: rpkt_eventtype_continuum 2'
  END IF
 END IF
 DEALLOCATE(actirrates%Lline, actirrates%nline)
 ! write(*,*) 'event_dist: eofloop, do_loop = ', do_loop
END DO
! STOP 'event_dist: testing'  
END SUBROUTINE event_dist
