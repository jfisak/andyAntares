SUBROUTINE init_photsphere(n_pack)

  USE types

  IMPLICIT NONE

  INTEGER                           :: I, n_pack, ind_cell_numb, ind_x, ind_y, ind_z 
  DOUBLE PRECISION                  :: L_star, sint, cost, sinp, cosp, length, freq, D
  !   DOUBLE PRECISION, PARAMETER       :: delta_t=1.D0
  DOUBLE PRECISION, DIMENSION(3)    :: direction, directionn

  L_star = 4.D0*pi*(R_star)**2*sigma*T_eff**4
  print*, L_star, pi, R_star/r_sun,sigma, T_eff

  !    ind_x = nx_cell/2 + 1
  !    ind_y = ny_cell/2 + 1
  !    ind_z = nz_cell/2 + 1
  !    ind_cell_numb = (ind_x-1)*ny_cell*nz_cell + (ind_y-1)*nz_cell + ind_z
  !    print*, ind_cell_numb
  !    print*, R_star
  DO I = 1, n_pack
     ! Place photon on the photosphere's surface
     CALL random_unitvector1(direction, sint, cost, sinp, cosp)
     package(I)%pos = R_star * direction

     ! Then give it a random direction outward from the photosphere
     CALL random_unitvector2(directionn) !random_unitvector(direction) 
     direction(1)=directionn(3)*sint*cosp+directionn(1)*cost*cosp-directionn(2)*sinp
     direction(2)=directionn(3)*sint*sinp+directionn(1)*cost*sinp+directionn(2)*cosp
     direction(3)=directionn(3)*cost-directionn(1)*sint
     package(I)%dir = direction

     ! Now put the photon to the corresponding grid cell
     ! Determine the cell index where is the photon 
     ! This works only for regular grids!!!!
     ind_x = FLOOR(package(I)%pos(1)/cell_width + DBLE(nx_cell)/2) + 1
     ind_y = FLOOR(package(I)%pos(2)/cell_width + DBLE(ny_cell)/2) + 1
     ind_z = FLOOR(package(I)%pos(3)/cell_width + DBLE(nz_cell)/2) + 1
     ind_cell_numb = (ind_x - 1) * ny_cell * nz_cell + (ind_y - 1) * nz_cell + ind_z
     IF ((ind_cell_numb .GT. nx_cell*ny_cell*nz_cell) .OR. (ind_cell_numb .LT. 1)) &
          STOP 'ERROR in cell_number'
     package(I)%cell_numb = ind_cell_numb

     ! Flag the packet as an active r-pkt and allow all kind of cell crossings
     package(I)%active     = 1
     package(I)%typ        = type_rpkt
     package(I)%last_cross = NONE

     ! Assign rf energy and frequency to the packet
     package(I)%e_rf = L_star/n_pack  
     IF (I .EQ. 1) print*, package(I)%e_rf
     CALL freq_from_planck(freq)   ! here the frequency is sampled from a Planck law
     package(I)%freq_rf = freq

     ! Now convert the energy and frequency to their cmf values
     CALL doppler_factor(I, D)
     package(I)%freq_cmf = package(I)%freq_rf * D 
     package(I)%e_cmf    = package(I)%e_rf * D  

     ! Assine 1 to the last_line whith which package is in resonance
     package(I)%last_line = no_line

     ! print*, package(I)%cell_numb,package(I)%dir !,  package(I)%pos, package(I)% e_rf
     ! e_cmf, freq_cmf, freq_rf, cell_numb, pack_numb, active
  
  END DO

  ! PRINT*, ind_x, ind_y, ind_z, ind_cell_numb

  ! I=10
  ! package(I)%pos = 95.
  ! package(I)%dir = -1.
  ! package(I)%e_rf = 0.

END SUBROUTINE init_photsphere

