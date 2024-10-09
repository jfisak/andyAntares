SUBROUTINE do_spectrum(n_pack)

USE types
USE constants

IMPLICIT NONE    

INTEGER                               :: I, n_pack, pack_index, nubin
DOUBLE PRECISION                      :: delta_nu, delta_e, freq, lambda, flambda, planck, ls_A
DOUBLE PRECISION, DIMENSION(n_nubin)  :: specflux, redspecflux, freqs
INTEGER, DIMENSION(n_nubin)           :: escs
CHARACTER(LEN=80)                     :: outputspecfile
INTEGER                               :: distribution = 0
DOUBLE PRECISION                      :: wale_start, wale_end 
DOUBLE PRECISION                      :: nu_max, nu_min

wale_start = 500   ! in Angstroms
wale_end = 20000   ! in Angstroms

nu_max = const_c / (wale_start * 1.D-8)
nu_min = const_c / (wale_end * 1.D-8)

write(outputspecfile, "(A, A13)") trim(outputfolder), "/spectrum.dat"
! Set up the frequency grid to extract spectrum
write(99,*) 'SETUP FREQ GRID'
delta_nu = (nu_max - nu_min) / n_nubin
DO I= 1, n_nubin 
   specflux(I) = 0.D0
   SELECT CASE(distribution)
   CASE(0)
    freqs(I) = nu_min + (I - 1) * delta_nu
   CASE(1)
    ! freqs(I) = 
   END SELECT
   escs = 0
   !write(99,*) i, spectrum(i)%freq, spectrum(i)%flux
END DO
  
! Loop over all packets
write(99,*) 'BIN PACKETS' 
DO pack_index = 1, n_pack
 ! And take all which actually escaped
 ! write(*,*) 'do_spectrum: pack_index = ', pack_index, ' typ = ', package(pack_index)%typ
 IF (package(pack_index)%typ .EQ. type_escaped) THEN
  freq = package(pack_index)%freq_rf
  ! Only bin those packets which are in the allowed frequency range
  IF ((freq .GT. nu_min) .AND. (freq .LT. nu_max)) THEN
   nubin = floor( (freq - nu_min) / delta_nu ) + 1
   ! put the star to 100 parsecs
   delta_e = (package(pack_index)%e_rf / delta_nu)! / (4.D0 * const_pi * (1.D2 * const_pc)**2)
   ! write(*,*) 'do_spectrum: e_rf = ', package(pack_index)%e_rf
   specflux(nubin) = specflux(nubin) + delta_e
   escs(nubin) = escs(nubin) + 1
  ENDIF
 END IF
END DO


! DO I = 1, n_nubin
!     frequency = spectrum(I)%freq
!     planck =( 2.D0 * const_h * frequency**3 / const_c**2  ) * (  1.D0 / ( EXP( (const_h * frequency) / (const_kB * T_eff) ) - 1.D0 )  )
!     WRITE(19,*)  frequency, spectrum(I)%flux, spectrum(I)%esc, planck
! END DO

 write(99,*) 'WRITE TO FILE'
 
 ls_A = const_c * 1.D8

#if mpi==1
 ! IF(my_rank /= 0) THEN
 CALL MPI_REDUCE(specflux, redspecflux, n_nubin, &
  MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
 ! END IF
 IF(my_rank == 0) THEN
  DO I = 1, n_nubin
   ! write(99,*) 'do_spectrum: specflux(I) = ', specflux(I)
   specflux(I) = redspecflux(I) / DBLE(n_tasks)
  END DO
#endif
 OPEN (UNIT=19, FILE=outputspecfile)     
  DO I = 1, n_nubin  
   lambda = ls_A / freqs(I)
   flambda = specflux(I) * ( ls_A / lambda**2 )
   planck = ( 2.D0 * const_h * ls_A**2 / lambda**5 ) * &
    ( 1.D0 / ( EXP( const_h * ls_A / (lambda * const_kB * T_eff) ) - 1.D0) )
   ! write(99,*) 'do_spectrum: lambda = ', lambda, ' flux = ', flambda
   write(19,*) lambda, flambda, planck, flambda/planck, escs(I)
  END DO
 CLOSE(19)
#if mpi==1
 END IF
#endif


END SUBROUTINE do_spectrum
