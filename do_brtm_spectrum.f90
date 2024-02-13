SUBROUTINE do_brtm_spectrum(n_vpacks, nu_min, nu_max, freqs, specflux)

USE types
USE constants
IMPLICIT NONE

INTEGER                                                 :: n_vpacks, pack_index
DOUBLE PRECISION                                        :: nu_min, nu_max
DOUBLE PRECISION, DIMENSION(n_nubin)                    :: specflux, freqs
DOUBLE PRECISION                                        :: delta_nu, freq, delta_e

INTEGER                                                 :: I, nubin


delta_nu = (nu_max - nu_min) / n_nubin

DO I= 1, n_nubin 
 specflux(I) = 0.D0
 freqs(I) = nu_min + (I - 1) * delta_nu
 ! escs = 0
 !write(99,*) i, spectrum(i)%freq, spectrum(i)%flux
END DO


DO pack_index = 1, n_vpacks
 ! And take all which actually escaped
 write(*,*) 'do_spectrum: pack_index = ', pack_index, ' typ = ', package(pack_index)%typ
 IF (package(pack_index)%typ == type_photosphere) THEN
  freq = package(pack_index)%freq_rf
  ! Only bin those packets which are in the allowed frequency range
  IF ((freq > nu_min) .AND. (freq < nu_max)) THEN
   nubin = floor( (freq - nu_min) / delta_nu ) + 1
   ! put the star to 100 parsecs
   !delta_e = (package(pack_index)%e_rf / delta_nu)! / (4.D0 * pi * (1.D2 * parsec)**2)
   delta_e = 100
   ! write(*,*) 'do_spectrum: e_rf = ', package(pack_index)%e_rf
   specflux(nubin) = specflux(nubin) + delta_e
  ENDIF
 END IF
END DO





END SUBROUTINE do_brtm_spectrum
