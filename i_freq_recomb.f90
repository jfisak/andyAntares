! this SBR will sample a random frequency of the packet for recombination
! it is sampled from integrals over b-f emission coefficient
! this procedure should be called by i_package if and only if n_points
! (number of photoionization data) does exist, because if it does not
! exist, the total rate is equal to zero thus no recombination deactivation
! is possible to happen
!
! INPUT: indexe(INT): the element index
!        indexi(INT): the ion index
!        indexl(INT): the level index
!        pack_index(INT): the index of the packet
! OUTPUT: ran_frequency(DBLE): a generated frequency
! 
SUBROUTINE i_freq_recomb(indexe, indexi, indexl, pack_index, ran_frequency)
USE types
USE constants
IMPLICIT NONE

INTEGER                                 :: indexl, pack_index
INTEGER                                 :: indexe, indexi
INTEGER                                 :: n_points
! model grid information
INTEGER                                 :: act_mgi
DOUBLE PRECISION                        :: temp
DOUBLE PRECISION                       :: ran2
! loop variables
INTEGER                                 :: I
INTEGER                                 :: get_package_model_index
! variables for integration
DOUBLE PRECISION                        :: act_freq, act_value, act_sum
DOUBLE PRECISION                        :: ran_frequency
DOUBLE PRECISION, ALLOCATABLE           :: exps(:), freqs(:), css(:), ints(:)
! integral value
DOUBLE  PRECISION                       :: int_value
DOUBLE PRECISION                        :: rand_z
INTEGER                                 :: act_point
! linear interpolation
DOUBLE PRECISION                        :: freq1, freq2, css1, css2, ali, bli
DOUBLE PRECISION                        :: freqt
INTEGER                                 :: Istart

n_points = SIZE(elements(indexe)%ions(indexi - 1)%levels(indexl)%photcros(1,:))
IF(n_points == 0) THEN
 ! write(*,*) 'i_package: n_points = 0'
 STOP 'i_freq_recomb: cannot generate a new r-packet frequency'
END IF
IF(n_points /= 0) THEN
 ALLOCATE(freqs(n_points), css(n_points),exps(n_points),ints(n_points))
END IF
freqs(:) = elements(indexe)%ions(indexi -1)%levels(indexl)%photcros(1,:)
css(:) = elements(indexe)%ions(indexi -1)%levels(indexl)%photcros(2,:)

freqt = (MINVAL(elements(indexe)%ions(indexI)%levels(:)%exci_energy) - &
 elements(indexe)%ions(indexi - 1)%levels(indexl)%exci_energy) / const_h
! write(*,*) 'i_freq_recomb: phfreq1 = ', elements(indexe)%ions(indexi - 1)%levels(indexl)%phfreq,&
! ' phfreq2 = ', freqt
! ran_frequency = freqt
! RETURN
Istart = 0
DO I = 1, n_points
 IF(freqs(I) > freqt) THEN
  Istart = I
  EXIT
 END IF
END DO

IF(Istart == 0) THEN
 write(*,*) 'i_freq_recomb: indexe = ', indexe, ' indexi = ', indexi
 write(*,*) ' indexl = ', indexl
 write(*,*) ' initial point was not found'
 STOP
END IF


! calculation of temperature
act_mgi = get_package_model_index(pack_index)
temp = model_grid(act_mgi)%t

! random number
rand_z = ran2(idum)
! write(*,*) 'i_freq_recomb: rand_z = ', rand_z

! saving field of exponentials, it will speed up the calculation procedure
DO I=1,n_points
 act_freq = freqs(I)
 exps(I) = exp(-( const_h * act_freq ) / ( const_kB * temp ))
 ! write(*,*) 'exps(I) = ', exps(I)
END DO
! calculation of the integral value
act_sum = 0.D0
ints(Istart) = 0.D0
DO I=Istart,n_points - 1
 act_value = 2.D0  * const_h / const_c**2 *&
  (css(I + 1) * freqs(I + 1)**3 * exps( I + 1) + css(I) * freqs(I)**3 * exps(I)) &
  * (freqs(I + 1) - freqs(I)) / 2.D0
 ints(I + 1) = act_sum + act_value
 act_sum = act_sum + act_value
 ! write(*,*) 'i_freq_recomb: act_sum = ', act_sum
END DO
! this is the right side of equation for the random frequency calculation
int_value = act_sum * rand_z
! write(*,*) 'i_freq_recomb: int_value = ', int_value
DEALLOCATE(exps)
 ! write(*,*) 'i_freq_recomb: integral value = ', int_value
! we can find frequency now
!  it is calculated in this way:
! we calculate integral from larger frequencies to lower frequencies
! and if the summed up value is larger than value z*int
! we found interval of frequency, where the new frequency is placed
act_sum = 0.D0
DO I=Istart + 1,n_points
 ! write(*,*) 'i_freq_recomb: act_sum = ', act_sum
 IF(ints(I) >= int_value) THEN
  act_point = I
  ! write(*,*) 'i_freq_recomb: act_point = ', I
  EXIT
 END IF
END DO
! now classical linear interpolation
!CALL lin_int(n_points, css, freqs, act_point, int_value, ran_frequency)
! if the frequency is out of range of the frequency interval
! the total rate will be equal to zero
css1 = ints(act_point)
css2 = ints(act_point - 1)
freq1 = freqs(act_point)
freq2 = freqs(act_point - 1)
ali = (freq1 - freq2) / (css1 - css2)
bli = (freq2 * css1 - freq1 * css2) / (css1 - css2)
ran_frequency = ali * int_value + bli
! write(77,*) ran_frequency
! write(*,*) 'i_freq_recomb: rand_z = ', rand_z
! write(*,*) 'i_freq_recomb: act_point = ', act_point, ' n_points = ', n_points
! write(*,*) 'i_freq_recomb: css1 = ', css1, ' css2 = ', css2,&
!  ' freq1 = ', freq1, ' freq2 = ', freq2, ' ali = ', ali, ' bli = ', bli
! write(*,*) 'i_freq_recomb: ran_frequency = ', ran_frequency
IF(ran_frequency < 0.D0) THEN
 write(*,*) 'i_freq_recomb: ran_frequency < 0'
 STOP
END IF


END SUBROUTINE i_freq_recomb
