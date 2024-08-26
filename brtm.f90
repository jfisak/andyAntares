! a basic sbr for the initialization of the backward ray tracing method
SUBROUTINE brtm()

USE types
USE constants
IMPLICIT NONE

DOUBLE PRECISION, DIMENSION(3)                          :: obs_point, ccd_point, ccd_centre
INTEGER                                                 :: cur_vpack
INTEGER, PARAMETER                                      :: Nvpackets = 500
! number of packet flown into the photosphere
INTEGER                                                 :: n_inside, n_outside
DOUBLE PRECISION, DIMENSION(3)                          :: cur_pos, cur_direction
DOUBLE PRECISION                                        :: wale_start, wale_end
DOUBLE PRECISION                                        :: nu_min, nu_max, ran_freq
DOUBLE PRECISION                                        :: ran2

DOUBLE PRECISION, DIMENSION(n_nubin)                    :: cur_spectrum, freqs

LOGICAL                                                 :: procout=.false.

! definition of a detector
INTEGER                                                 :: det_nu, det_nv, cur_ccd, det_tot_nuv
INTEGER                                                 :: det_cur_nu, det_cur_nv, cur_nu
DOUBLE PRECISION                                        :: det_lu, det_lv
DOUBLE PRECISION, DIMENSION(3)                          :: det_vec_u, det_vec_v
DOUBLE PRECISION                                        :: det_cell_wu, det_cell_wv
DOUBLE PRECISION, DIMENSION(3)                          :: uvmin
DOUBLE PRECISION, ALLOCATABLE                           :: det_matrix(:,:), det_spectra(:,:)

INTEGER                                                 :: my_ccd_start, my_ccd_end
INTEGER                                                 :: N_single, N_zbytek

DOUBLE PRECISION                                        :: ccdc_phi, ccdc_rad, ccdc_theta
DOUBLE PRECISION                                        :: obs_ccd_dist

write(99,*) '___________________________________________________________________'
write(99,*) '___________________________________________________________________'
write(99,*) '______________BACKWARD RAY TRACING METHOD__________________________'
write(99,*) '___________________________________________________________________'
write(99,*) '___________________________________________________________________'

! to start calculations we have to deallocate the package array firstly
DEALLOCATE(package)

wale_start = 200   ! in Angstroms
wale_end = 400   ! in Angstroms

nu_max = light_speed / (wale_start * 1.D-8)
nu_min = light_speed / (wale_end * 1.D-8)

! a temporary definition of a detector
! number of points in each CCD chip
det_nu = 20
det_nv = 20
det_tot_nuv = det_nu * det_nv

! a size of a detector
det_lu = 5.0
det_lv = 5.0


ALLOCATE(det_matrix(det_nu, det_nv), det_spectra(det_tot_nuv, n_nubin))
! observing point
! obs_point = (/ -R_inf  ,  0.D0,  0.D0 /)
obs_ccd_dist = 0.2*sqrt(det_lu**2+det_lv**2)
! ccd_centre = (/ -R_inf/2.0,  R_inf/2.D0,  R_inf/4.D0 /)
ccd_centre = (/ -R_inf,  0.D0,  0.D0 /)
write(48,*) ccd_centre
obs_point = ccd_centre + obs_ccd_dist * ccd_centre/norm2(ccd_centre)
write(48,*) obs_point

ccdc_rad = sqrt(ccd_centre(1)**2+ccd_centre(2)**2+ccd_centre(3))
ccdc_theta = acos(ccd_centre(3)/ccdc_rad)
write(*,*) 'brtm: ccdc_rad = ', ccdc_rad, ' ccdc_theta = ', ccdc_theta
! phi is more complicated to calculate
IF(ccd_centre(1) > 0.0) THEN
 IF(ccd_centre(2) >= 0.0) THEN
  ccdc_phi = atan(ccd_centre(1)/ccd_centre(1))
 ELSE IF(ccd_centre(2) < 0.0) THEN
  ccdc_phi = atan(ccd_centre(1)/ccd_centre(1)) + 2.0 * pi
 END IF
ELSE IF(ccd_centre(1) == 0.0) THEN
 IF(ccd_centre(2) > 0.0) THEN
  ccdc_phi = pi/2.0
 ELSE IF(ccd_centre(2) < 0.0) THEN
  ccdc_phi = 3.0 * pi/2.0
 END IF
ELSE IF(ccd_centre(1) < 0.0) THEN
 ccdc_phi = atan(ccd_centre(1)/ccd_centre(1)) + pi
END IF

write(*,*) 'brtm: ccdc_rad = ', ccdc_rad/R_star, ' ccdc_phi = ', ccdc_phi, ' ccdc_theta = ', ccdc_theta



! !!! only a temporary solution !!!
! a calculation of the vectors u and v
det_vec_u = (/ cos(ccdc_theta) * cos(ccdc_phi), cos(ccdc_theta) * sin(ccdc_phi), -sin(ccdc_theta) /)
det_vec_v = (/ -sin(ccdc_phi), cos(ccdc_phi), 0.D0 /)
write(31,*) ccd_centre, det_vec_u
write(31,*) ccd_centre, det_vec_v
write(31,*) ccd_centre, ccd_centre
! a size of a single cell
det_cell_wu = det_lv / DBLE(det_nu)
det_cell_wv = det_lu / DBLE(det_nv)

! lower coordinates of a ccd chip
uvmin = ccd_centre - 0.5D00 * (det_vec_u * det_lu + det_vec_v * det_lv)
write(48,*) uvmin

! 
#if mpi == 1
 N_single = det_tot_nuv/n_tasks
 N_zbytek = det_tot_nuv - n_tasks * N_single
 IF(my_rank <= N_zbytek - 1) THEN
  my_ccd_start = my_rank * (N_single + 1) + 1
  my_ccd_end = my_rank * (N_single + 1) + N_single
 ELSE IF(N_zbytek == 0) THEN
  my_ccd_start = my_rank * (N_single + 1) + 1
  my_ccd_end = my_rank * (N_single + 1) + N_single
 ELSE
  my_ccd_start = N_zbytek * (N_single + 1) + (my_rank - N_zbytek - 1) * N_single + 1
  my_ccd_end = N_zbytek * (N_single + 1) + (my_rank - N_zbytek - 1) * N_single + N_single +1
 END IF
#else
 my_ccd_start = 1
 my_ccd_end = det_tot_nuv
#endif
write(*,*) 'brtm: N_single = ', N_single, ' N_zbytek = ', N_zbytek
write(*,*) 'brtm: my_ccd_start = ', my_ccd_start, ' my_ccd_end = ', my_ccd_end
! then one by one we will be sending packets through the CCD chip
DO cur_ccd = my_ccd_start, my_ccd_end
 ! write(*,*) 'brtm: cur_ccd = ', cur_ccd
 det_cur_nv = INT((cur_ccd - 1)/det_nu) + 1
 det_cur_nu = INT(cur_ccd - (det_cur_nv -1) * det_nu)
 
 ccd_point = uvmin + det_cell_wu * det_vec_u * (det_cur_nu - 0.5D0) + &
   & det_cell_wv * det_vec_v * (det_cur_nv - 0.5D0)
 
 write(*,*) 'brtm: ccd_point = ', ccd_point(2) - ccd_centre(2), ccd_point(3) - ccd_centre(3)

 cur_direction = (ccd_point - obs_point)/norm2(ccd_point - obs_point)
 write(46,*) obs_point, ccd_point - obs_point
 write(47,*) ccd_point


! END DO  
!  ! ccd_point = 
! 
! STOP 'brtm: testing'
!  
! DO 
 n_inside = 0
 n_outside = 0
 
 ALLOCATE(package(Nvpackets + 1))
 
 
 DO cur_vpack = 1, Nvpackets
 
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
  ELSE
   n_outside = n_outside + 1
   package(cur_vpack)%typ = type_escaped
  END IF
  
 
 END DO ! a loop over virtual packets
 write(*,*) 'brtm: det_cur_nv = ', det_cur_nv, ' det_cur_nu = ', det_cur_nu, ' n_inside = ', n_inside, ' n_outside = ', n_outside
 
 CALL do_brtm_spectrum(Nvpackets, nu_min, nu_max, freqs, cur_spectrum)
 
 ! saving output
 ! saving a temporary spectrum into a variable
 ! DO cur_specpoint = 1, n_nubin
 !  det_spectra(cur_ccd, cur_specpoint) = cur_spectrum(cur_specpoint)
 ! END DO

 ! calculation the absorbed/sent ratio
 det_matrix(det_cur_nu, det_cur_nv) = DBLE(n_inside)/DBLE(det_tot_nuv)

 ! cleaning procedures
 DEALLOCATE(package)

END DO ! loop over detector cells

OPEN(449,FILE='ccd_matrix.dat')
DO cur_nu = 1, det_nu
 write(449,*) det_matrix(cur_nu,:)
END DO
CLOSE(449)

END SUBROUTINE brtm
