SUBROUTINE event_dist(pack_index, cell_dist, e_dist, event, Lcont)

  USE types

  IMPLICIT NONE    

  INTEGER                           :: I, pack_index, event,do_loop, get_package_model_index
  INTEGER                           :: nextLine, indexe, indexi, lower_level, current_mgi
! pointer to a field of continuum rates
  DOUBLE PRECISION                  :: e_dist, ran_numb, ran2, tau_rand, cell_dist, D
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
  DOUBLE PRECISION, POINTER               :: Lcont(:)

  my_rank = OMP_GET_THREAD_NUM()
  !write(*,*)'event_dist: Thread rank: ', my_rank
       
  n_pack_d = SIZE(package)
  dummypackage = n_pack_d - n_dummy_packs + my_rank + 1


10  ran_numb = ran2(idum)  ! PUT IT IN SUBROUTINE - write is as do loop
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

  ! Number density of hydrogen. Assuming that all hydrogen is at the ground state we can use this as population number
  ! pop_number = model_grid(get_package_model_index(pack_index))%rho/mp_g    
                         !cell(package(pack_index)%cell_numb)%model_index)%rho/mp_g  

  electron_density = model_grid(current_mgi)%e_dens
  IF(current_mgi .EQ. n_modelgrid + 2) electron_density = 0.D0
  !print*, 'electron_density = ', electron_density

  ! For now we neglect cont. opacities, but the routine was written generally to
  ! allow for adding  cont. opacity in the future
  ! kappa_cont = 0.D0                      ! no e scattering
  CALL r_kappa_cont(pack_index, kappa_cont)
  !kappa_cont = sigma_e * electron_density  ! add e scattering

  ! This is the opacity in co-moving frame. Must be transformed to the lab frame
  ! According to Mihalas and Mihalas Eq. 90.8 this is achieved by 
  CALL doppler_factor(pack_index, D)
  kappa_cont = D * kappa_cont

  ! initialization of n_next_lines to be equal to one
!  n_next_lines = 1
  DO WHILE (do_loop .EQ. 1) 

 CALL next_line(1, pack_index, nextLine, n_next_lines)
 freq_line = linelist(nextLine)%freq
 indexe = linelist(nextLine)%indexe
 indexi = linelist(nextLine)%indexi
 lower_level = linelist(nextLine)%lower
     
   !  print*, 'subroutine event_dist:'
!     IF(pack_index.EQ.1) print*, pack_index, package(pack_index)%last_line, package(pack_index)%freq_cmf, nextLine, freq_line

     IF (package(pack_index)%freq_cmf .GT. freq_line) THEN
        CALL resonance_distance(pack_index, freq_line, l_dist)
        ! near future:
        !CALL l_dist()
        !print*, 'AAAAA', ' l_dist = ', l_dist
        ! Calculate optical depth in the next line (Sobolev, dv/dr dependent)
        ! and continuum optical depth accumulated up to the line
        !CALL velo(pack_index,vel_vec)
        !tau_line = light_speed/freq_line * constant * osc_line * pop_number * &
        !     vec_length(package(pack_index)%pos)/vec_length(vel_vec)     
        !tau_cont = kappa_cont * l_dist

        ! Assine package(pack_index) to package(dummypackage) 
        ! Use dummypackage to move to the l_dist only to check whether or which kind of interaction will happen at l_dist 
        ! We move only dummypackage instead of package(pack_index) only to check whether or which kind of interaction will happen at l_dist 
        ! 
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
!        print*, 'ABC', current_mgi, element, package(pack_index)%pos, &
!        SQRT(package(pack_index)%pos(1)**2 + package(pack_index)%pos(2)**2 + package(pack_index)%pos(3)**2)/R_star

        !print*, 'event_dist: n_next_lines = ', n_next_lines
        tau_line = 0.D0
        DO I = 1, n_next_lines
         indexe = linelist(nextLine + I - 1)%indexe
         indexi = linelist(nextLine + I - 1)%indexi
         lower_level = linelist(nextLine + I - 1)%lower
         CALL populations(indexe, indexi, lower_level, current_mgi, pop_number)
         f_ul = linelist(nextLine + I - 1)%f_ul
         tau_line = tau_line + light_speed / freq_line * constant * linelist(nextLine)%f_ul * pop_number * &
          vec_length(package(dummypackage)%pos) / vec_length(vel_vec)     
        END DO
        IF(current_mgi .EQ. n_modelgrid + 2) tau_line = 0.D0
        tau_cont = kappa_cont * l_dist
        IF(current_mgi .EQ. n_modelgrid + 2) tau_cont = 0.D0
        !print*, 'event_dist:', tau_line, tau_cont

        ! Now do a step by step analysis of which event occurs and return the 
        ! distance and corresponding event
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
           END IF
        ELSE
           ! Continuum process will happen
           e_dist = dist + (tau_rand - tau) / kappa_cont
           do_loop = 0 
           event = rpkt_eventtype_continuum
!           print*, 'cont.process happens',  tau_line, tau_cont
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
!           print*, 'cont.process happens' tau_line, tau_cont
        END IF
     END IF
  END DO
  
END SUBROUTINE event_dist
