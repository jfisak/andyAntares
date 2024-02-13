! a basic sbr for the initialization of the backward ray tracing method
SUBROUTINE brtm()

USE types
USE constants
IMPLICIT NONE

DOUBLE PRECISION, DIMENSION(3)                          :: obs_point, ccd_point
INTEGER                                                 :: cur_vpack
INTEGER, PARAMETER                                      :: Nvpackets = 1000000
! number of packet flown into the photosphere
INTEGER                                                 :: n_inside, n_outside
DOUBLE PRECISION, DIMENSION(3)                          :: cur_pos
DOUBLE PRECISION                                        :: wale_start, wale_end
DOUBLE PRECISION                                        :: nu_min, nu_max, ran_freq
DOUBLE PRECISION                                        :: ran2

DOUBLE PRECISION, DIMENSION(n_nubin)                    :: cur_spectrum, freqs

INTEGER                                                 :: I

wale_start = 200   ! in Angstroms
wale_end = 20000   ! in Angstroms

nu_max = light_speed / (wale_start * 1.D-8)
nu_min = light_speed / (wale_end * 1.D-8)

DEALLOCATE(package)

write(99,*) '___________________________________________________________________'
write(99,*) '___________________________________________________________________'
write(99,*) '______________BACKWARD RAY TRACING METHOD__________________________'
write(99,*) '___________________________________________________________________'
write(99,*) '___________________________________________________________________'

! observing point
obs_point = (/ -2*R_inf  ,  0.D0,  0.D0 /)
                                    
ccd_point = (/ -1.9*R_inf,  0.D0,  0.D0 /)

n_inside = 0
n_outside = 0

ALLOCATE(package(Nvpackets + 1))
! CCD chip


DO cur_vpack = 1, Nvpackets

 write(*,*) 'brtm: processing the v-packet: cur_vpack = ', cur_vpack
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

write(*,*) 'brtm: n_inside = ', n_inside, 'n_outside = ', n_outside

DO I = 1, n_nubin
 write(39,*) freqs(I), cur_spectrum(I)
END DO

END SUBROUTINE brtm
