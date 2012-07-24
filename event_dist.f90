SUBROUTINE event_dist(pack_index, cell_dist, e_dist, event)

  USE types

  IMPLICIT NONE    

  INTEGER                           :: pack_index, event,do_loop, get_package_model_index
  DOUBLE PRECISION                  :: e_dist, ran_numb, ran2, tau_rand, cell_dist, D
  DOUBLE PRECISION                  :: tau, l_dist, tau_line, constant, pop_number, tau_cont
  DOUBLE PRECISION                  :: electron_density, kappa_cont, vec_length, dist
  DOUBLE PRECISION, DIMENSION(30)   :: vel_vec

10  ran_numb = ran2(idum)  ! PUT IT IN SUBROUTINE - write is as do loop
    IF (ran_numb .EQ. 0.D0) GOTO 10    
    tau_rand = -LOG(ran_numb)

  tau = 0.D0
  dist = 0.D0
  do_loop = 1
  constant = (pi * e_charge**2)/( me_g * light_speed)
  ! Number density of hydrogen. Assuming that all hydrogen is at the graund state we can use this as population number
  pop_number = model_grid(get_package_model_index(pack_index))%rho/mp_g    
                         !cell(package(pack_index)%cell_numb)%model_index)%rho/mp_g  

  ! For now we neglect cont. opacities, but the routine was written generally to
  ! allow for adding  cont. opacity in the future
  electron_density = 0.D0 * pop_number !0.D0 !take this from modelgrid finally
  kappa_cont = sigma_e * electron_density
  ! This is the opacity in co-moving frame. Must be transformed to the lab frame
  ! According to Mihalas and Mihalas Eq. 90.8 this is achieved by 
  CALL doppler_factor(pack_index, D)
  kappa_cont = D * kappa_cont

  DO WHILE (do_loop .EQ. 1) 
     IF (package(pack_index)%freq_cmf .GT. freq_line) THEN
        ! Calculate distance the photon needs to travel to come to
        ! resonance with the next line. This assumes homologus
        ! expansion i.e. velocity is proportional to r. Projected
        ! gradient of the projected velocity in a direction of the
        ! photon propagation is more complicated in case of no
        ! homologus expannsion
        l_dist = light_speed * (R_inf/V_inf) * ((package(pack_index)%freq_cmf - freq_line)/freq_line)

        ! Calculate optical depth in the next line (Sobolev, dv/dr dependent)
        ! and continuum optical depth accumulated up to the line
        CALL velo(pack_index,vel_vec)
        tau_line = light_speed/freq_line * constant * osc_line * pop_number * &
             vec_length(package(pack_index)%pos)/vec_length(vel_vec)     
        tau_cont = kappa_cont * l_dist

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
           e_dist = dist + (tau_rand - tau)/kappa_cont
           do_loop = 0 
           event = rpkt_eventtype_continuum
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
           e_dist = dist + (tau_rand - tau)/kappa_cont
           do_loop = 0                 
           event = rpkt_eventtype_continuum
        END IF
     END IF
  END DO
  
END SUBROUTINE event_dist
