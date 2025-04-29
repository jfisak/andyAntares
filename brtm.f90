! a basic sbr for the initialization of the backward ray tracing method
!
! INPUT: NONE
! OUTPUT: NONE
! 
SUBROUTINE brtm()

USE MPI
USE types
USE constants
IMPLICIT NONE

DOUBLE PRECISION, DIMENSION(const_dimofspace)                          :: obs_point, ccd_point, ccd_centre
INTEGER                                                 :: cur_vpack
INTEGER, PARAMETER                                      :: Nvpackets = 50000
! number of packet flown into the photosphere
INTEGER                                                 :: n_inside, n_outside
DOUBLE PRECISION, DIMENSION(const_dimofspace)                          :: cur_pos, cur_direction
DOUBLE PRECISION                                        :: wale_start, wale_end
DOUBLE PRECISION                                        :: nu_min, nu_max, ran_freq
DOUBLE PRECISION                                        :: ran2

DOUBLE PRECISION, DIMENSION(n_nubin)                    :: cur_spectrum, freqs

LOGICAL                                                 :: procout=.false.

! definition of a detector
INTEGER                                                 :: det_nu, det_nv, cur_ccd, det_tot_nuv
INTEGER                                                 :: det_cur_nu, det_cur_nv, cur_nu
DOUBLE PRECISION                                        :: det_lu, det_lv
DOUBLE PRECISION, DIMENSION(const_dimofspace)                          :: det_vec_u, det_vec_v
DOUBLE PRECISION                                        :: det_cell_wu, det_cell_wv
DOUBLE PRECISION, DIMENSION(const_dimofspace)                          :: uvmin
DOUBLE PRECISION, ALLOCATABLE                           :: det_matrix(:,:), det_spectra(:,:), det_sending(:,:)

INTEGER                                                 :: my_ccd_start, my_ccd_end
INTEGER                                                 :: N_single, N_zbytek

DOUBLE PRECISION                                        :: ccdc_phi, ccdc_rad, ccdc_theta
DOUBLE PRECISION                                        :: obs_ccd_dist

DOUBLE PRECISION, DIMENSION(n_nubin)                    :: specflux, sendflux
DOUBLE PRECISION                                        :: delta_nu, delta_e, freq
INTEGER                                                 :: ind_I, nubin, pack_index
DOUBLE PRECISION                                        :: rand_u, rand_v
INTEGER                                                 :: N_tot_zbytek
LOGICAL                                                 :: ccd_mode
INTEGER                                                 :: n_crossed, n_vpacks
INTEGER, PARAMETER                                      :: max_n_crossed = 2000000




write(99,*) '___________________________________________________________________'
write(99,*) '___________________________________________________________________'
write(99,*) '______________BACKWARD RAY TRACING METHOD__________________________'
write(99,*) '___________________________________________________________________'
write(99,*) '___________________________________________________________________'

! to start calculations we have to deallocate the package array firstly
DEALLOCATE(package)

ccd_mode = .false.

wale_start = 6000   ! in Angstroms
wale_end = 7000   ! in Angstroms

nu_max = const_c / (wale_start * 1.D-8)
nu_min = const_c / (wale_end * 1.D-8)

! a temporary definition of a detector
! number of points in each CCD chip
det_nu = 100
det_nv = 100
det_tot_nuv = det_nu * det_nv

! a size of a detector
det_lu = 5.0
det_lv = 5.0

delta_nu = (nu_max - nu_min) / n_nubin

sendflux(:) = 0.D0
DO ind_I= 1, n_nubin 
 freqs(ind_I) = nu_min + (ind_I - 1) * delta_nu
END DO

ALLOCATE(det_sending(det_nu, det_nv), det_matrix(det_nu, det_nv))!, det_spectra(det_tot_nuv, n_nubin))
det_sending(:,:) = 0.D0
! observing point
! obs_point = (/ -R_inf  ,  0.D0,  0.D0 /)
obs_ccd_dist = 0.2*sqrt(det_lu**2+det_lv**2)
! ccd_centre = (/ -R_inf/2.0,  R_inf/2.D0,  R_inf/4.D0 /)
ccd_centre = (/ -0.25*R_inf,  -.25*R_inf,  -.25*R_inf /)
! write(48,*) ccd_centre
obs_point = ccd_centre + obs_ccd_dist * ccd_centre/norm2(ccd_centre)
! write(48,*) obs_point

ccdc_rad = sqrt(ccd_centre(ind_x)**2+ccd_centre(ind_y)**2+ccd_centre(ind_z))
ccdc_theta = acos(ccd_centre(ind_z)/ccdc_rad)
! write(*,*) 'brtm: ccdc_rad = ', ccdc_rad, ' ccdc_theta = ', ccdc_theta
! phi is more complicated to calculate
IF(ccd_centre(ind_x) > 0.0) THEN
 ccdc_phi = atan(ccd_centre(ind_y)/ccd_centre(ind_x))
ELSE IF(ccd_centre(ind_x) == 0.0) THEN
 IF(ccd_centre(ind_y) > 0.0) THEN
  ccdc_phi = const_pi/2.0
 ELSE IF(ccd_centre(ind_y) < 0.0) THEN
  ccdc_phi = -const_pi/2.0
 END IF
ELSE IF(ccd_centre(ind_x) < 0.0) THEN
 IF(ccd_centre(ind_y) >= 0.0) THEN
  ccdc_phi = atan(ccd_centre(ind_y)/ccd_centre(ind_x)) + const_pi
 ELSE IF(ccd_centre(ind_y) < 0.0) THEN
  ccdc_phi = atan(ccd_centre(ind_y)/ccd_centre(ind_x)) - const_pi
 END IF
END IF

! write(*,*) 'brtm: ccdc_rad = ', ccdc_rad/R_star, ' ccdc_phi = ', ccdc_phi, ' ccdc_theta = ', ccdc_theta



! !!! only a temporary solution !!!
! a calculation of the vectors u and v
det_vec_u = (/ cos(ccdc_theta) * cos(ccdc_phi), cos(ccdc_theta) * sin(ccdc_phi), -sin(ccdc_theta) /)
det_vec_v = (/ -sin(ccdc_phi), cos(ccdc_phi), 0.D0 /)
! write(31,*) ccd_centre, det_vec_u
! write(31,*) ccd_centre, det_vec_v
! write(31,*) ccd_centre, ccd_centre
! a size of a single cell
det_cell_wu = det_lu / DBLE(det_nu)
det_cell_wv = det_lv / DBLE(det_nv)

! lower coordinates of a ccd chip
uvmin = ccd_centre - 0.5D00 * (det_vec_u * det_lu + det_vec_v * det_lv)
! write(48,*) uvmin

! 
#if mpi == 1
 N_single = det_tot_nuv/n_tasks
 N_zbytek = det_tot_nuv - n_tasks * N_single
 IF(N_zbytek /= 0) N_tot_zbytek = (N_zbytek + 1) * (N_single + 1) + N_zbytek
 IF(my_rank <= N_zbytek - 1) THEN
  my_ccd_start = my_rank * (N_single + 1) + my_rank
  my_ccd_end = (my_rank + 1) * (N_single + 1) + N_single
 ELSE IF(N_zbytek == 0) THEN
  my_ccd_start = my_rank * N_single + 1
  my_ccd_end = (my_rank + 1) * N_single
 ELSE IF(my_rank > N_zbytek - 1) THEN
  my_ccd_start = N_tot_zbytek + 1 + (my_rank - N_zbytek) * N_single + (my_rank - N_zbytek)
  my_ccd_end = N_tot_zbytek + 1 + (my_rank - N_zbytek + 1) * N_single + (my_rank - N_zbytek)
 END IF
 IF(my_rank == n_tasks - 1) THEN
  my_ccd_end = det_tot_nuv
 END IF
#else
 my_ccd_start = 1
 my_ccd_end = det_tot_nuv
#endif
! write(*,*) 'brtm: N_single = ', N_single, ' N_zbytek = ', N_zbytek
! write(*,*) 'brtm: my_ccd_start = ', my_ccd_start, ' my_ccd_end = ', my_ccd_end
! then one by one we will be sending packets through the CCD chip
cur_ccd = 0
n_crossed = 0
n_vpacks = 0
DO 
 ! conditions to stop the loop according to the current computing mode
 IF(ccd_mode) THEN ! going along the ccd detector one point by one
  IF(cur_ccd == 0) THEN
   cur_ccd = my_ccd_start
  ELSE IF(cur_ccd >= my_ccd_start .and. cur_ccd < my_ccd_end) THEN
   cur_ccd = cur_ccd + 1
  ELSE IF(cur_ccd == my_ccd_end) THEN
   exit
  END IF
  rand_u = ran2(idum) * det_cell_wu
  rand_v = ran2(idum) * det_cell_wv
  ccd_point = uvmin + det_cell_wu * det_vec_u * (det_cur_nu - 1 + rand_u) + &
    & det_cell_wv * det_vec_v * (det_cur_nv - 1 + rand_v)
  det_cur_nv = INT((cur_ccd - 1)/det_nu) + 1
  det_cur_nu = INT(cur_ccd - (det_cur_nv -1) * det_nu)
 ELSE ! random positions in the detector, the calculation is stopped until a critical
      ! number of packets crossing the photosphere is reached
  IF(n_crossed >= max_n_crossed) THEN
   EXIT
  END IF
 END IF
 ! write(*,*) 'brtm: cur_ccd = ', cur_ccd
 
 
 ! write(*,*) 'brtm: ccd_point = ', ccd_point(ind_y) - ccd_centre(ind_y), ccd_point(ind_z) - ccd_centre(ind_z)

 cur_direction = (ccd_point - obs_point)/norm2(ccd_point - obs_point)
 ! write(46,*) obs_point, ccd_point - obs_point
 ! write(47,*) ccd_point

 n_inside = 0
 n_outside = 0
 
 ALLOCATE(package(Nvpackets + 1))
 
 write(*,*) 'brtm: n_crossed = ', n_crossed
 
 DO cur_vpack = 1, Nvpackets
  n_vpacks = n_vpacks + 1
  IF(.not. ccd_mode) THEN
   rand_u = ran2(idum) * det_lu
   rand_v = ran2(idum) * det_lv
   ! write(*,*) 'brtm: rand_u = ', rand_u, ' rand_v = ', rand_v
   ccd_point = uvmin + det_vec_u * rand_u + &
     &  det_vec_v * rand_v
   det_cur_nv = FLOOR(rand_u/det_lu * det_nu) + 1
   det_cur_nu = FLOOR(rand_v/det_lv * det_nv) + 1
   ! write(*,*) 'brtm: ccd_point = ', (ccd_point - uvmin)/det_cell_wu
   ! write(*,*) 'brtm: det_cur_nv = ', det_cur_nv, ' det_cur_nu = ', det_cur_nu
  END IF
 
  if(procout) write(*,*) 'brtm: processing the v-packet: cur_vpack = ', cur_vpack
  ! a basic initialisation of a packet
  package(cur_vpack)%pos = obs_point
  package(cur_vpack)%dir = (ccd_point - obs_point)/norm2(ccd_point - obs_point)
 
  ran_freq = nu_min + (nu_max - nu_min) * ran2(idum)
  package(cur_vpack)%freq_rf = ran_freq
  
  ! basic properties of a packet
  package(cur_vpack)%active = 1
  package(cur_vpack)%typ = type_rpkt
  package(cur_vpack)%n_interactions = 0
  package(cur_vpack)%next_cross = NONE
  package(cur_vpack)%e_rf = 1.e10
  package(cur_vpack)%last_line = no_line
  package(cur_vpack)%delta_s = 0.D0
  package(cur_vpack)%virtual = .TRUE.
  
  ! calling a sbr to process a v-packet
  CALL do_vpackage(cur_vpack)
  
  ! here we must find out, whether the packet is inside the photosphere, or outside the grid
  cur_pos = package(cur_vpack)%pos
  
  ! counting number of packets flown into the photosphere
  IF(norm2(cur_pos) < R_star) THEN
   n_inside = n_inside + 1
   package(cur_vpack)%typ = type_photosphere
   n_crossed = n_crossed + 1
  ELSE
   n_outside = n_outside + 1
   package(cur_vpack)%typ = type_escaped
  END IF
  
  IF(.not. ccd_mode) THEN
   det_sending(det_cur_nu, det_cur_nv) = DBLE(n_crossed)/DBLE(n_vpacks)
  END IF
 
 END DO ! a loop over virtual packets
 ! write(*,*) 'brtm: det_cur_nv = ', det_cur_nv, ' det_cur_nu = ', det_cur_nu, ' n_inside = ', n_inside, ' n_outside = ', n_outside
 
 IF(ccd_mode) THEN
  det_sending(det_cur_nu, det_cur_nv) = DBLE(n_inside)/DBLE(n_vpacks)
 END IF
 ! CALL do_brtm_spectrum(Nvpackets, nu_min, nu_max, freqs, cur_spectrum, ccd_centre, det_nu, det_nv)
 DO pack_index = 1, Nvpackets
  ! And take all which actually escaped
  ! write(*,*) 'do_spectrum: pack_index = ', pack_index, ' typ = ', package(pack_index)%typ
  IF(package(pack_index)%typ == type_photosphere) THEN
   freq = package(pack_index)%freq_rf
   ! Only bin those packets which are in the allowed frequency range
   IF ((freq > nu_min) .AND. (freq < nu_max)) THEN
    nubin = floor( (freq - nu_min) / delta_nu ) + 1
    ! put the star to 100 parsecs
    delta_e = package(pack_index)%e_rf * 4.0/const_pi * (norm2(obs_point)**2 * det_nu * det_nv)
    ! write(*,*) 'do_spectrum: e_rf = ', package(pack_index)%e_rf
    sendflux(nubin) = sendflux(nubin) + delta_e
   ENDIF
  END IF
 END DO
 
 ! saving output
 ! saving a temporary spectrum into a variable
 ! DO cur_specpoint = 1, n_nubin
 !  det_spectra(cur_ccd, cur_specpoint) = cur_spectrum(cur_specpoint)
 ! END DO

 ! calculation the absorbed/sent ratio

 ! cleaning procedures
 DEALLOCATE(package)

END DO ! loop over detector cells

IF(n_tasks > 1) THEN
 ! DO ind_I = 1, det_nu
  CALL MPI_ALLREDUCE(det_sending(1:det_nu,1:det_nv), det_matrix(1:det_nu,1:det_nv), det_nu * det_nv, MPI_DOUBLE, MPI_SUM, &
  mpi_comm_world, ierr)
 ! END DO
 CALL MPI_ALLREDUCE(sendflux(:), specflux(:), n_nubin, MPI_DOUBLE, MPI_SUM, &
  mpi_comm_world, ierr)
 det_matrix(1:det_nu,1:det_nv) = det_matrix(1:det_nu,1:det_nv)/n_tasks
END IF
IF(my_rank == 0) THEN
 OPEN(449,FILE='ccd_matrix.dat')
 DO cur_nu = 1, det_nu
  write(449,*) det_matrix(cur_nu,:)
 END DO
 CLOSE(449)
 
 OPEN(450, FILE='ccd_spectrum.dat')
  DO ind_I = 1, n_nubin
   write(450, *) 1.D8 * const_c/freqs(ind_I), specflux(ind_I)
  END DO
 CLOSE(450)
END IF




END SUBROUTINE brtm
