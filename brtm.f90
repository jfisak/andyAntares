! a basic sbr for the initialization of the backward ray tracing method
SUBROUTINE brtm()

USE types
USE constants
IMPLICIT NONE

DOUBLE PRECISION, DIMENSION(3)                          :: obs_point, ccd_point, ccd_centre
INTEGER                                                 :: cur_vpack
INTEGER, PARAMETER                                      :: Nvpackets = 50000
! number of packet flown into the photosphere
INTEGER                                                 :: n_inside, n_outside
DOUBLE PRECISION, DIMENSION(3)                          :: cur_pos
DOUBLE PRECISION                                        :: wale_start, wale_end
DOUBLE PRECISION                                        :: nu_min, nu_max, ran_freq
DOUBLE PRECISION                                        :: ran2

INTEGER                                                 :: I

DOUBLE PRECISION, DIMENSION(n_nubin)                    :: cur_spectrum, freqs

! definition of a detector
INTEGER                                                 :: det_nu, det_nv, cur_ccd, det_tot_nuv
INTEGER                                                 :: det_cur_nu, det_cur_nv
DOUBLE PRECISION                                        :: det_lu, det_lv
DOUBLE PRECISION, DIMENSION(3)                          :: det_vec_u, det_vec_v
DOUBLE PRECISION                                        :: det_cell_wu, det_cell_wv
DOUBLE PRECISION, DIMENSION(3)                          :: uvmin

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
det_nu = 10
det_nv = 10
det_tot_nuv = det_nu * det_nv

! observing point
obs_point = (/ -2*R_inf  ,  0.D0,  0.D0 /)
ccd_centre = (/ -2*R_inf - 1.5D1,  0.D0,  0.D0 /)

! a size of a detector
det_lu = 5
det_lv = 5

! !!! only a temporary solution !!!
! a calculation of the vectors u and v
det_vec_u = (/ 0, 1, 0 /)
det_vec_v = (/ 0, 0, 1 /)

! a size of a single cell
det_cell_wu = det_lv / DBLE(det_nu)
det_cell_wv = det_lu / DBLE(det_nv)

! lower coordinates of a ccd chip
uvmin = ccd_centre - 0.5D00 * (det_vec_u * det_lu + det_vec_v * det_lv)
! then one by one we will be sending packets through the CCD chip
DO cur_ccd = 1, det_tot_nuv
 det_cur_nv = INT((cur_ccd - 1)/det_nu) + 1
 det_cur_nu = cur_ccd - (det_cur_nv -1) * det_nu
 
 ccd_point = uvmin + det_cell_wu * det_vec_u * (det_cur_nu + 0.5D0) + &
   & det_cell_wv * det_vec_v * (det_cur_nv + 0.5D0)
 
 ! write(*,*) 'brtm: ccd_point = ', ccd_point(2)/det_lu, ccd_point(3)/det_lv

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
 
  ! write(*,*) 'brtm: processing the v-packet: cur_vpack = ', cur_vpack
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
 
 CALL do_brtm_spectrum(Nvpackets, nu_min, nu_max, freqs, cur_spectrum)
 
 write(*,*) 'brtm: nx, ny = ', det_cur_nu, det_cur_nv, ' n_inside = ', n_inside, 'n_outside = ', n_outside
 
 ! DO I = 1, n_nubin
 !  write(39,*) light_speed/ (freqs(I) * 1.D-8), cur_spectrum(I)
 ! END DO

 DEALLOCATE(package)

END DO ! loop over detector cells

END SUBROUTINE brtm
