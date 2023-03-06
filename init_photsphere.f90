SUBROUTINE init_photsphere(n_pack)

  USE types
USE constants

  IMPLICIT NONE

  INTEGER                           :: I, J, n_pack, ind_cell_numb
  DOUBLE PRECISION                  :: L_star, sint, cost, sinp, cosp, freq, D
  !   DOUBLE PRECISION, PARAMETER       :: delta_t=1.D0
  DOUBLE PRECISION, DIMENSION(3)    :: direction, directionn
  DOUBLE PRECISION, DIMENSION(n_pack) :: frequencies

  DOUBLE PRECISION, DIMENSION(3)        :: cur_pos
  INTEGER                               :: cur_pgi, cur_mgi
  DOUBLE PRECISION                      :: cur_Teff
  LOGICAL, PARAMETER                    :: homogeneous=.true.

  DOUBLE PRECISION                      :: R_bound, ran2

  destroyed_pack = 0
  L_star = 4.D0*pi*(R_star)**2*sigma*T_eff**4
  write(99,*) 'init photsphere...'
  write(99,*) 'init_photsphere: R_star = ', R_star, ' T_eff = ', T_eff, ' L_star = ', L_star
  !print*, L_star, pi, R_star/r_sun,sigma, T_eff

  ! delete
  R_bound = R_star

  !    ind_x = nx_cell/2 + 1
  !    ind_y = ny_cell/2 + 1
  !    ind_z = nz_cell/2 + 1
  !    ind_cell_numb = (ind_x-1)*ny_cell*nz_cell + (ind_y-1)*nz_cell + ind_z
  !    write(99,*) ind_cell_numb
  !    write(99,*) R_star
! OPEN(16,FILE='photon_positions.dat')
  DO I = 1, n_pack
   ! Place photon on the photosphere's surface
   CALL random_unitvector1(direction, sint, cost, sinp, cosp)
   package(I)%pos = R_bound * direction
  
   IF(I > tot_saved_packets) THEN
    ! write(*,*) 'init_photsphere: R_star = ', R_star

    ! Then give it a random direction outward from the photosphere
    CALL random_unitvector2(directionn) !random_unitvector(direction) 
    direction(1)=directionn(3)*sint*cosp+directionn(1)*cost*cosp-directionn(2)*sinp
    direction(2)=directionn(3)*sint*sinp+directionn(1)*cost*sinp+directionn(2)*cosp
    direction(3)=directionn(3)*cost-directionn(1)*sint
    package(I)%dir = direction

    ! Now put the photon to the corresponding grid cell
    ! Determine the cell index where is the photon 
    ! This works only for regular grids!!!!
    ! write(*,*) 'init_photsphere: calling find_dyn_cell1, pack_index = ', I
    CALL find_dyn_cell1(package(I)%pos,ind_cell_numb)
    IF ((inputflux .EQ. 0) ) THEN
     IF(homogeneous) THEN
      CALL freq_from_planck(freq, T_eff)   ! here the frequency is sampled from a Planck law
     END IF
     IF(I > tot_saved_packets) package(I)%freq_rf = freq
    ELSE IF ((inputflux .EQ. 1 .OR. inputflux == 2) .AND. (I==1)) THEN
    CALL freq_from_file(n_pack,frequencies) ! frequency is sampled using an existing emergent flux
    DO J = tot_saved_packets + 1, n_pack
     package(J)%freq_rf = frequencies(J)
    END DO
   END IF
    ! write(*,*) 'init_photsphere: init cell numb = ', ind_cell_numb
    IF(ind_cell_numb > SIZE(dyn_cell)) THEN
     write(*,*) 'init_photsphere: wrong cell number'
     CALL abort()
    END IF
    package(I)%cell_numb = ind_cell_numb

    ! Flag the packet as an active r-pkt and allow all kind of cell crossings
    package(I)%active     = 1
    package(I)%typ        = type_rpkt
    package(I)%n_interactions = 0
    package(I)%next_cross = NONE

    ! Assign rf energy and frequency to the packet
    package(I)%e_rf = L_star/n_pack  

    ! Now convert the energy and frequency to their cmf values
    CALL doppler_factor(I, D)
    package(I)%freq_cmf = package(I)%freq_rf * D 
    package(I)%e_cmf    = package(I)%e_rf * D  

    ! Assine 1 to the last_line whith which package is in resonance
    package(I)%last_line = no_line
    package(I)%delta_s = 0.D0
   END IF
  END DO

OPEN(19,file="photonFdistr.dat")
 do I=1,n_pack
  write(19,*) -99, package(I)%freq_rf, package(I)%e_rf
 end do
CLOSE(19)
        

END SUBROUTINE init_photsphere

