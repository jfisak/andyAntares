! sets the packets up
! 
! INPUT: n_pack(INT): number of packets
! OUTPUT: NONE
!
SUBROUTINE init_photsphere(n_pack)

USE types
USE constants

IMPLICIT NONE

INTEGER                           :: ind_I, ind_J, n_pack, ind_cell_numb
DOUBLE PRECISION                  :: L_star, sint, cost, sinp, cosp, freq, doppler_D
DOUBLE PRECISION, DIMENSION(const_dimofspace)    :: direction, directionn
DOUBLE PRECISION, DIMENSION(n_pack) :: frequencies

INTEGER                               :: cur_mgi
DOUBLE PRECISION                      :: cur_Teff
LOGICAL                               :: homogeneous

DOUBLE PRECISION                      :: R_bound
INTEGER, PARAMETER                    :: ind_savephdistr = 12
LOGICAL                               :: is_diff

IF(T_eff <= 0) THEN
 L_star = 0.D0
 homogeneous = .false.
ELSE
 L_star = 4.D0*const_pi*(R_star)**2*const_stefbolz*T_eff**4
 homogeneous = .true.
END IF

destroyed_pack = 0
write(99,*) 'init photsphere...'
write(99,*) 'init_photsphere: R_star = ', R_star, ' T_eff = ', T_eff, ' L_star = ', L_star

! delete
R_bound = R_star

DO ind_I = 1, n_pack
 ! Place photon on the photosphere's surface
 CALL random_unitvector1(direction, sint, cost, sinp, cosp)
 package(ind_I)%pos = R_bound * direction

 IF(ind_I > tot_saved_packets) THEN
  ! write(*,*) 'init_photsphere: R_star = ', R_star

  ! Then give it a random direction outward from the photosphere
  CALL random_unitvector2(directionn) !random_unitvector(direction) 
  direction(ind_x)=directionn(ind_z)*sint*cosp+directionn(ind_x)*cost*cosp-directionn(ind_y)*sinp
  direction(ind_y)=directionn(ind_z)*sint*sinp+directionn(ind_x)*cost*sinp+directionn(ind_y)*cosp
  direction(ind_z)=directionn(ind_z)*cost-directionn(ind_x)*sint
  package(ind_I)%dir = direction

  ! Now put the photon to the corresponding grid cell
  ! Determine the cell index where is the photon 
  ! This works only for regular grids!!!!
  ! write(*,*) 'init_photsphere: calling find_dyn_cell1, pack_index = ', I
  CALL find_dyn_cell1(package(ind_I)%pos,ind_cell_numb)
  IF ((inputflux .EQ. 0) ) THEN
   IF(homogeneous) THEN
    CALL freq_from_planck(freq, T_eff)   ! here the frequency is sampled from a Planck law
   ELSE
    cur_mgi = dyn_cell(ind_cell_numb)%model_index
    cur_Teff = model_grid(cur_mgi)%T
    L_star = L_star + 4.D0*const_pi*(R_star)**2*const_stefbolz*T_eff**4/n_pack
    if(cur_mgi > n_modelgrid) then
     cur_Teff = -T_eff
    end if
    ! write(64,*) cur_Teff
    CALL freq_from_planck(freq, cur_Teff)
   END IF
   IF(ind_I > tot_saved_packets) package(ind_I)%freq_rf = freq
  ELSE IF ((inputflux .EQ. 1 .OR. inputflux == 2) .AND. (ind_I==1)) THEN
  CALL freq_from_file(n_pack,frequencies) ! frequency is sampled using an existing emergent flux
  DO ind_J = tot_saved_packets + 1, n_pack
   package(ind_J)%freq_rf = frequencies(ind_J)
  END DO
 END IF
  ! write(*,*) 'init_photsphere: init cell numb = ', ind_cell_numb
  IF(ind_cell_numb > SIZE(dyn_cell)) THEN
   write(*,*) 'init_photsphere: wrong cell number'
   CALL abort()
  END IF
  package(ind_I)%cell_numb = ind_cell_numb

  ! Flag the packet as an active r-pkt and allow all kind of cell crossings
  package(ind_I)%active     = 1
  package(ind_I)%typ        = type_rpkt
  package(ind_I)%n_interactions = 0
  package(ind_I)%next_cross = NONE
  package(ind_I)%virtual = .FALSE.

  IF(enable_diffusion == 1) THEN
   cur_mgi = dyn_cell(ind_cell_numb)%model_index
   is_diff = model_grid(cur_mgi)%is_difapp
   IF(is_diff) THEN
    package(ind_I)%typ = type_dpkt
   END IF
  END IF

  ! Assign rf energy and frequency to the packet
  package(ind_I)%e_rf = L_star/n_pack  

  ! Now convert the energy and frequency to their cmf values
  ! CALL doppler_factor(I, R_star * direction, direction, D)
  CALL doppler_factor(ind_I, doppler_D)
  package(ind_I)%freq_cmf = package(ind_I)%freq_rf * doppler_D 
  package(ind_I)%e_cmf    = package(ind_I)%e_rf * doppler_D  

  ! Assine 1 to the last_line whith which package is in resonance
  package(ind_I)%last_line = no_line
  package(ind_I)%delta_s = 0.D0
 END IF
END DO

! CALL save_output(ind_savephdistr)
        

END SUBROUTINE init_photsphere

