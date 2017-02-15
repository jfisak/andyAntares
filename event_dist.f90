SUBROUTINE event_dist(pack_index, cell_dist, e_dist, event)

  USE types

  IMPLICIT NONE    

  INTEGER                           :: I, pack_index, event,do_loop, get_package_model_index
  INTEGER                           :: next_line, indexe, indexi, lower_level, current_mgi
  DOUBLE PRECISION                  :: e_dist, ran_numb, ran2, tau_rand, cell_dist, D
  DOUBLE PRECISION                  :: tau, l_dist, tau_line, constant, pop_number, tau_cont
  DOUBLE PRECISION                  :: electron_density, kappa_cont, vec_length, dist
  DOUBLE PRECISION                  :: graund_level_pop, g_gl, g_ll, e_exc
  DOUBLE PRECISION, DIMENSION(30)   :: vel_vec

10  ran_numb = ran2(idum)  ! PUT IT IN SUBROUTINE - write is as do loop
    IF (ran_numb .EQ. 0.D0) GOTO 10    
    tau_rand = -LOG(ran_numb)

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
   kappa_cont = 0.D0                      ! no e scattering
  ! kappa_cont = sigma_e * electron_density * 10.D2 ! add e scattering

  ! This is the opacity in co-moving frame. Must be transformed to the lab frame
  ! According to Mihalas and Mihalas Eq. 90.8 this is achieved by 
  CALL doppler_factor(pack_index, D)
  kappa_cont = D * kappa_cont

  DO WHILE (do_loop .EQ. 1) 

     IF (package(pack_index)%last_line .EQ. no_line) THEN
         ! In case we have more lines
         DO I = 1, ntransitions
            ! If 
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

     next_line = package(pack_index)%last_line + 1
     freq_line = linelist(next_line)%freq
     indexe = linelist(next_line)%indexe
     indexi = linelist(next_line)%indexi
     lower_level = linelist(next_line)%lower
     
   !  print*, 'subroutine event_dist:'
     IF(pack_index.EQ.1) print*, pack_index, package(pack_index)%last_line, package(pack_index)%freq_cmf, next_line, freq_line

     IF (package(pack_index)%freq_cmf .GT. freq_line) THEN
        ! Calculate distance the photon needs to travel to come to
        ! resonance with the next line. This assumes homologous
        ! expansion i.e. velocity is proportional to r. Projected
        ! gradient of the projected velocity in a direction of the
        ! photon propagation is more complicated in case of no
        ! homologous expansion
        ! This is case when we assume that he have only hydrogen 
        l_dist = light_speed * (R_inf/V_inf) * ((package(pack_index)%freq_cmf - freq_line)/package(pack_index)%freq_rf)
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

        !This is not safe yet: possible mismatch of atomic number and element index!!!
        !print*, 'current_mgi = ', current_mgi, ' indexe = ', indexe, ' indexi = ', indexi
        graund_level_pop = model_grid(current_mgi)%grid_comp(indexe)%grid_ion(indexi)%gl_pop
!        print*, 'BBBBB'
        g_gl = elements(indexe)%ions(indexi)%levels(1)%stat_waight
        g_ll = elements(indexe)%ions(indexi)%levels(lower_level)%stat_waight
        e_exc = elements(indexe)%ions(indexi)%levels(lower_level)%exci_energy - elements(indexe)%ions(indexi)%levels(1)%exci_energy
        !print*, 'e_exc1 = ', elements(indexe)%ions(indexi)%levels(lower_level)%exci_energy, 'e_exc2 = ', &
        !elements(indexe)%ions(indexi)%levels(1)%exci_energy
        pop_number = graund_level_pop * g_ll / g_gl * exp(e_exc / BOLK / model_grid(current_mgi)%T )
        !print*, 'e_exc1 = ', elements(indexe)%ions(indexi)%levels(lower_level)%exci_energy, 'e_exc2 = ', &
        !        elements(indexe)%ions(indexi)%levels(1)%exci_energy
        !Needs proper treatment of empty cells
        !IF (vec_length(package(dummypackage)%pos) .GT. R_inf) pop_number = 0
        tau_line = light_speed / freq_line * constant * linelist(next_line)%f_ul * pop_number * &
                   vec_length(package(dummypackage)%pos) / vec_length(vel_vec)     
        IF(current_mgi .EQ. n_modelgrid + 2) tau_line = 0.D0
        tau_cont = kappa_cont * l_dist
        !print*, 'CCCCC', tau_line, tau_cont
 

        ! Now do a step by step analysis of which event occurs and return the 
        ! distance and corresponding event
        IF ((tau_rand - tau) .GT. tau_cont) THEN
           IF ((tau_rand - tau) .GT. (tau_cont + tau_line)) THEN
              tau = tau + tau_cont + tau_line
              dist = dist + l_dist
              IF (dist .GT. cell_dist) THEN
                 ! In this case the package propagates to the next cell
                 e_dist = cell_dist + 1.D20
                 do_loop = 0
                 event = rpkt_eventtype_changecell
              END IF
           ELSE
              e_dist = dist + l_dist
              do_loop = 0
              event = rpkt_eventtype_lineinteraction
           END IF
        ELSE
           ! Continuum process will happens
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
