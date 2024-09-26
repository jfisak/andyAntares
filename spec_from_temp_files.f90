PROGRAM spec_from_temp_files

IMPLICIT NONE

! definition of photon
TYPE photon 
   INTEGER                         :: cell_numb, pack_numb, active
   DOUBLE PRECISION                :: e_cmf, e_rf, freq_cmf, freq_rf, delta_s
   INTEGER                         :: typ, next_cross, last_line
   INTEGER                         :: n_interactions
   DOUBLE PRECISION, DIMENSION(3)  :: pos, dir 
   INTEGER                         :: l_ele, l_ion, l_lev
END TYPE photon
INTEGER, PARAMETER                 :: type_escaped=-99 

TYPE(photon), ALLOCATABLE          :: package(:)
TYPE(photon)                                 :: cur_package
INTEGER                            :: n_file, tot_saved_packets, cur_pack, rank_num
INTEGER                                 :: reading_packets
CHARACTER(LEN=100)                       :: outputfolder, temp_file_name
CHARACTER(30)                            :: temp_filename = 'temp_packet'
CHARACTER(3)                            :: n_ranks
INTEGER                                 :: nranks
LOGICAL                                 :: file_exists
INTEGER                                 :: n_packs
! variables fo spectra creation
INTEGER, PARAMETER                    :: n_nubin = 1000
INTEGER                               :: I, pack_index, nubin
DOUBLE PRECISION                      :: delta_nu, delta_e, freq, lambda, flambda, planck, ls_A
DOUBLE PRECISION, DIMENSION(n_nubin)  :: specflux, redspecflux, freqs
INTEGER, DIMENSION(n_nubin)           :: escs
CHARACTER(LEN=80)                       :: outputspecfile
CHARACTER(LEN=5)                        :: Teff, min_wale, max_wale
INTEGER                                 :: T_eff, minwale, maxwale
DOUBLE PRECISION                        :: nu_max, nu_min
INTEGER, DIMENSION (9)                  :: TT
! physical constants
  DOUBLE PRECISION, PARAMETER        :: pi=3.1415926535897932D+00,&
                                        const_me_g=9.109534D-28,&
                                        const_mp_g=1.6726485D-24,&
                                        const_sigma_e=6.6516D-25,&
                                        h=6.626176D-27,&
                                        light_speed=2.99792458D+10,&
                                        e_charge=4.803242D-10,&
                                        ftran=0.6407D+00, &       
                                        nio=4.5655967D+14,&
                                        const=1.D-04,&
                                        vel_ter=920.0D+05,&
                                        r_sun=695990.D+05,&
                                        beta=2.11638D+00,  &   
                                        BOLK=1.380662D-16,&
                                        m_sun=1.989D+33,&
                                        sigma =5.6704D-05 !ergcm^(-2)s(-1)K(-4) !D. H. Cohen et al.2012
  DOUBLE PRECISION, PARAMETER        :: parsec=30.857D17, const_ev = 1.60217646D-12, saha_const=2.0706839D-16, b = 1.D0



n_file = 0
tot_saved_packets = 0
cur_pack = 0
rank_num = 0

!_____________________________________________________________________________________________________________________________
! READING IMPORTANT VARIABLES TYPED BY USER
! reads name of output folder
CALL GET_ENVIRONMENT_VARIABLE("OUTPUTFO", outputfolder)
! name of output folder
IF(outputfolder == '') THEN
 STOP 'variable OUTPUTFO is not set'
END IF
! number of ranks
CALL GET_ENVIRONMENT_VARIABLE("NRANKS", n_ranks)
IF(n_ranks == '') THEN
 STOP 'variable NRANKS is not set'
END IF
READ(n_ranks,*) nranks
! effective temperature
CALL GET_ENVIRONMENT_VARIABLE("TEFF", Teff)
IF(Teff == '') THEN
 STOP 'variable TEFF is not set'
END IF
READ(Teff,*) T_eff
! maximal wawelength
CALL GET_ENVIRONMENT_VARIABLE("MXWW", max_wale)
IF(max_wale == '') THEN
 maxwale = 7000
ELSE
 READ(maxwale,*) max_wale
END IF
! minimal wawelength
CALL GET_ENVIRONMENT_VARIABLE("MNWW", min_wale)
IF(min_wale == '') THEN
 minwale = 900
ELSE
 READ(minwale,*) min_wale
END IF
! calculation of maximal and minimal frequency
nu_max = light_speed / (DBLE(minwale) * 1.D-8)
nu_min = light_speed / (DBLE(maxwale) * 1.D-8)
write(*,*) 'nu_max = ', nu_max, ' nu_min = ', nu_min

! calculates number of packets
DO rank_num = 0, nranks - 1
 DO
  n_file = n_file + 1
  write(temp_file_name,"(A, A1, A, I3.3, A1, I6.6, A4)") TRIM(outputfolder), '/', TRIM(temp_filename), rank_num, "_", n_file, ".dat"
  !write(temp_file_name,"(A, I3.3, A1, I6.6, A4)") TRIM(temp_filename), rank_num, "_", n_file, ".dat"
  inquire(FILE=temp_file_name, EXIST=file_exists)
  ! write(*,*) 'find_unfinished_run: rank_num = ', rank_num, ' n_file = ', n_file
  ! write(*,*) 'find_unfinished_run: filename = ', temp_file_name, ' E? = ', file_exists
  IF(file_exists) THEN
   OPEN(1, FORM="unformatted", FILE=temp_file_name)
    DO
     READ(1, iostat=reading_packets) cur_package
     IF(reading_packets /= 0) EXIT
     n_packs = n_packs + 1
     ! write(*,*) n_packs
    END DO
   CLOSE(1)
  ELSE
   EXIT
  END IF
 END DO
 n_file = 0
END DO

write(*,*) 'allocating array package with dimension = ', n_packs
ALLOCATE(package(n_packs))

DO rank_num = 0, nranks - 1
 DO
  n_file = n_file + 1
  write(temp_file_name,"(A, A1, A, I3.3, A1, I6.6, A4)") TRIM(outputfolder), '/', TRIM(temp_filename), rank_num, "_", n_file, ".dat"
  !write(temp_file_name,"(A, I3.3, A1, I6.6, A4)") TRIM(temp_filename), rank_num, "_", n_file, ".dat"
  inquire(FILE=temp_file_name, EXIST=file_exists)
  ! write(*,*) 'find_unfinished_run: rank_num = ', rank_num, ' n_file = ', n_file
  ! write(*,*) 'find_unfinished_run: filename = ', temp_file_name, ' E? = ', file_exists
  IF(file_exists) THEN
   OPEN(21, FORM="unformatted", FILE=temp_file_name)
   ! OPEN(21, FILE=temp_file_name)
    DO
     READ(21, iostat=reading_packets) cur_package
     ! write(*,*) cur_package
     IF(reading_packets /= 0) EXIT
     cur_pack = cur_pack + 1
     IF(cur_pack > SIZE(package)) STOP 'number of loaded packets is larger than number of packets for a computation'
     package(cur_pack) = cur_package
     tot_saved_packets = tot_saved_packets + 1
     ! write(*,*) tot_saved_packets
    END DO
    ! write(*,*) 'find_unfinished_run: tot_saved_packets = ', tot_saved_packets
   CLOSE(21)
  END IF
  IF(.NOT. file_exists) EXIT
 END DO
 n_file = 0
END DO

! finally create the temporary spectrum
delta_nu = (nu_max - nu_min) / n_nubin
DO I= 1, n_nubin 
   freqs(I) = nu_min + (I - 1) * delta_nu
   specflux(I) = 0.D0
   escs = 0
   !write(99,*) i, spectrum(i)%freq, spectrum(i)%flux
END DO
  
! Loop over all packets
DO pack_index = 1, n_packs
 ! And take all which actually escaped
 ! write(*,*) 'do_spectrum: pack_index = ', pack_index, ' typ = ', package(pack_index)%typ
 IF (package(pack_index)%typ .EQ. type_escaped) THEN
  freq = package(pack_index)%freq_rf
  ! Only bin those packets which are in the allowed frequency range
  IF ((freq .GT. nu_min) .AND. (freq .LT. nu_max)) THEN
   nubin = floor( (freq - nu_min) / delta_nu ) + 1
   ! write(*,*) 'nubin = ', nubin
   ! put the star to 100 parsecs
   delta_e = (package(pack_index)%e_rf / delta_nu) / (4.D0 * pi * (1.D2 * parsec)**2)
   specflux(nubin) = specflux(nubin) + delta_e
   escs(nubin) = escs(nubin) + 1
  END IF
 END IF
END DO


! DO I = 1, n_nubin
!     frequency = spectrum(I)%freq
!     planck =( 2.D0 * h * frequency**3 / light_speed**2  ) * (  1.D0 / ( EXP( (h * frequency) / (BOLK * T_eff) ) - 1.D0 )  )
!     WRITE(19,*)  frequency, spectrum(I)%flux, spectrum(I)%esc, planck
! END DO

 ls_A = light_speed * 1.D8

 CALL DATE_AND_TIME(VALUES = TT)
 write(outputspecfile,"(A, A5, i4.4, I2.2, I2.2, I2.2, I2.2, I2.2, A4)") &
  TRIM(outputfolder), '/spec', TT(1), TT(2), TT(3), TT(5), TT(6), TT(7) ,".dat"
 OPEN (UNIT=19, FILE=outputspecfile)     
  DO I = 1, n_nubin  
   lambda = ls_A / freqs(I)
   flambda = specflux(I) * ( ls_A / lambda**2 )
   planck = ( 2.D0 * h * ls_A**2 / lambda**5 ) * &
    ( 1.D0 / ( EXP( h * ls_A / (lambda * BOLK * DBLE(T_eff)) ) - 1.D0) )
   ! write(99,*) 'do_spectrum: lambda = ', lambda, ' flux = ', flambda
   write(19,*) lambda, flambda, planck, flambda/planck, escs(I)
  END DO
 CLOSE(19)


END PROGRAM spec_from_temp_files
