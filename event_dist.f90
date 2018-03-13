SUBROUTINE event_dist(pack_index, cell_dist, e_dist, event, actirrates)

  USE types
  USE rates_r

  IMPLICIT NONE    

  INTEGER                           :: I, pack_index, event,do_loop, get_package_model_index
  INTEGER                           :: nextLine, indexe, indexi, lower_level, current_mgi
! pointer to a field of continuum rates
  DOUBLE PRECISION                  :: e_dist, ran_numb, tau_rand, cell_dist, D
  DOUBLE PRECISION                  :: tau, l_dist, tau_line, constant, pop_number, tau_cont
  DOUBLE PRECISION                  :: electron_density, kappa_cont, vec_length, dist
  DOUBLE PRECISION                  :: graund_level_pop, g_gl, g_ll, e_exc
  DOUBLE PRECISION                  :: f_ul
  DOUBLE PRECISION, PARAMETER       :: largeNumber = 1.D20
  ! number of lines with the same frequencies
  INTEGER                           :: n_next_lines
  DOUBLE PRECISION, DIMENSION(3)    :: vel_vec
  INTEGER                           :: approximation
  INTEGER                           :: my_rank, dummypackage
  INTEGER                           :: n_pack_d
  INTEGER                           :: OMP_GET_THREAD_NUM
  DOUBLE PRECISION                  :: freq_line
  TYPE(rrates)                      :: actirrates
  REAL(8)                           :: random
  ! looking for next line
  INTEGER                           :: act_line
  DOUBLE PRECISION                  :: summ, tot_lop
  DOUBLE PRECISION                  :: Blu! , exci_energy_l, exci_energy_u
  DOUBLE PRECISION                      :: low_pop, upp_pop 
  DOUBLE PRECISION                      :: stat_weight_l, stat_weight_u

  my_rank = OMP_GET_THREAD_NUM()
  !write(*,*)'event_dist: Thread rank: ', my_rank
       
  n_pack_d = SIZE(package)
  dummypackage = n_pack_d - n_dummy_packs + my_rank + 1


10  ran_numb = DBLE(random())  ! PUT IT IN SUBROUTINE - write is as do loop
    IF (ran_numb .EQ. 0.D0) GOTO 10    
    tau_rand = -LOG(ran_numb)

! allocating the field

  ! Initialize optical depth and distance
  tau = 0.D0
  dist = 0.D0
  do_loop = 1
  constant = (pi * e_charge**2)/( me_g * light_speed)

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
!  n_next_lines = 1
DO WHILE (do_loop .EQ. 1) 

 CALL next_line(1, pack_index, nextLine, n_next_lines)
 ! n_next_lines = 1
 ! write(*,*) 'event_dist: nextLine = ', nextLine, ' n_next_lines = ', n_next_lines
 freq_line = linelist(nextLine)%freq
 indexe = linelist(nextLine)%indexe
 indexi = linelist(nextLine)%indexi
 lower_level = linelist(nextLine)%lower
     
 IF (package(pack_index)%freq_cmf .GT. freq_line) THEN
  CALL resonance_distance(pack_index, freq_line, l_dist)

  ! Calculate optical depth in the next line (Sobolev, dv/dr dependent)
  ! and continuum optical depth accumulated up to the line
  package(dummypackage) = package(pack_index)
  !print*, 'before moving package #', pack_index
  CALL move_package(dummypackage, l_dist)
  CALL velo(dummypackage,vel_vec)

  !! Total mass density 
  !pop_number = M_dot / (4.D0 * pi * & 
  !     vec_length(package(dummypackage)%pos)**2 * vec_length(vel_vec))     
  !! Total number density of hydrogen
  !pop_number = pop_number/mp_g

  ! Current model grid cell
  ! Be sure to get the proper connection between element and elementindex, dito for ion
!  print*, 'ABC', current_mgi, element, package(pack_index)%pos, &
!  SQRT(package(pack_index)%pos(1)**2 + package(pack_index)%pos(2)**2 + package(pack_index)%pos(3)**2)/R_star

  !print*, 'event_dist: n_next_lines = ', n_next_lines
  tau_line = 0.D0
  DO I = 1, n_next_lines
   indexe = linelist(nextLine + I - 1)%indexe
   indexi = linelist(nextLine + I - 1)%indexi
   lower_level = linelist(nextLine + I - 1)%lower
   CALL populations(indexe, indexi, linelist(nextLine + I - 1)%lower, current_mgi, low_pop)
   CALL populations(indexe, indexi, linelist(nextLine + I - I)%upper, current_mgi, upp_pop)
   f_ul = linelist(nextLine + I - 1)%f_ul
   ! tau_line = tau_line + light_speed / freq_line * constant * &
   !            linelist(nextLine)%f_ul * pop_number * &
   ! vec_length(package(dummypackage)%pos) / vec_length(vel_vec)     
   ! the basic variables
   stat_weight_u = elements(indexe)%ions(indexi)%levels(linelist(nextLine + I - 1)%upper)%stat_waight
   stat_weight_l = elements(indexe)%ions(indexi)%levels(linelist(nextLine + I - I)%lower)%stat_waight
   ! exci_energy_u = elements(indexe)%ions(indexi)%levels(linelist(I)%upper)%exci_energy
   ! exci_energy_l = elements(indexe)%ions(indexi)%levels(linelist(I)%lower)%exci_energy
   Blu = light_speed**2.0 / (2.0 * h * linelist(I)%freq**3.0) * stat_weight_u / stat_weight_l &
     * linelist(I)%A_ul
   tau_line = tau_line +  low_pop * Blu * h * light_speed * (R_inf / V_inf) &
     / (4.0 * pi) * (1.D0 - (stat_weight_l * upp_pop) / (stat_weight_u * low_pop))
  END DO
  IF(current_mgi .EQ. n_modelgrid + 2) tau_line = 0.D0
  tau_cont = kappa_cont * l_dist
  IF(current_mgi .EQ. n_modelgrid + 2) tau_cont = 0.D0
  !print*, 'event_dist:', tau_line, tau_cont

  ! Now do a step by step analysis of which event occurs and return the 
  ! distance and corresponding event
  ! print*, 'event_dist: cont.process happens', tau_line, tau_cont
  IF ((tau_rand - tau) .GT. tau_cont) THEN
   IF ((tau_rand - tau) .GT. (tau_cont + tau_line)) THEN
    tau = tau + tau_cont + tau_line
    dist = dist + l_dist
    IF (dist .GT. cell_dist) THEN
     ! In this case the package propagates to the next cell
     e_dist = cell_dist + largeNumber
     do_loop = 0
     event = rpkt_eventtype_changecell
    END IF
   ELSE
    e_dist = dist + l_dist
    do_loop = 0
    event = rpkt_eventtype_lineinteraction
    ! choosing the next line
    IF(n_next_lines > 1) THEN
     tot_lop = 0.D0
     DO I = 1, n_next_lines
      tot_lop = tot_lop + linelist(nextLine + I - 1)%A_ul
     END DO
     ran_numb = DBLE(random()) * tot_lop
     summ = 0.D0
     DO I = 1, n_next_lines
      act_line = nextLine + I - 1
      IF(ran_numb > summ .AND. ran_numb < summ + linelist(act_line)%A_ul) THEN
       package(pack_index)%last_line = nextLine
      END IF
     END DO
    ELSE
     package(pack_index)%last_line = nextLine
    END IF 
   END IF
  ELSE
   ! Continuum process will happen
   e_dist = dist + (tau_rand - tau) / kappa_cont
   do_loop = 0 
   event = rpkt_eventtype_continuum
!   print*, 'cont.process happens',  tau_line, tau_cont
  END IF
 ELSE    
  ! The package cmf frequency is too red to interact to another
  ! line - No line interact anymore
  tau_cont = kappa_cont * (cell_dist - dist)
  IF ((tau_rand - tau) .GT. tau_cont) THEN
   e_dist = cell_dist + 1.D20
   do_loop = 0   
   event = rpkt_eventtype_changecell
  ELSE
   ! Continuum absorption happens
   e_dist = dist + (tau_rand - tau) / kappa_cont
   do_loop = 0                 
   event = rpkt_eventtype_continuum
  END IF
 END IF
END DO
  
END SUBROUTINE event_dist
