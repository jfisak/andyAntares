SUBROUTINE k_freq_fb(pack_index, act_proc, ran_freq, actikrates)
USE types
USE rates_k
IMPLICIT NONE

! input
INTEGER                                         :: pack_index, act_proc
DOUBLE PRECISION                                :: ran_freq
! ion
INTEGER                                         :: indexe, indexi, indexl
! photoionization data
INTEGER                                         :: nfreq
DOUBLE PRECISION, ALLOCATABLE                   :: freq(:), cross(:)
DOUBLE PRECISION                                :: ran2
! initial point and frequency
INTEGER                                         :: initPoint
DOUBLE PRECISION                                :: init_freq
! integral calculations
DOUBLE PRECISION, ALLOCATABLE                   :: func(:), intval(:)
DOUBLE PRECISION                                :: summ
DOUBLE PRECISION                                :: integral1
INTEGER                                         :: actIndex
! propagation grid informations
INTEGER                                         :: cur_mgi, get_package_model_index
DOUBLE PRECISION                                :: el_dens, temp, x, act_pop
! random numbers
DOUBLE PRECISION                                :: ran_num
! loop variables
INTEGER                                         :: I
! linear interpolation
DOUBLE PRECISION                                :: ali, bli, int1, int2, func1, func2
TYPE(krates)                                    :: actikrates

! informations about ion
! write(*,*) 'k_freq_fb: act_proc = ', act_proc, ' allocated? Lcfb = ', ALLOCATED(actikrates%Lcool_fbind)
indexe = INT(actikrates%Lcool_fbE(1, act_proc))
indexi = INT(actikrates%Lcool_fbE(2, act_proc))
indexl = INT(actikrates%Lcool_fbE(3, act_proc))
initPoint = INT(actikrates%Lcool_fbE(5, act_proc))
! write(*,*) 'k_freq_fb: initPoint = ', initPoint
! getting the photoionization cross section
! nfreq cannot be equal to zero, because a process with a zero rate could
! not be chosen in the previous step
nfreq = SIZE(elements(indexe)%ions(indexi)%levels(indexl)%photcros(1,:))
! write(*,*) 'k_freq_fb: el = ', indexe, ' ion = ', indexi, ' lev = ', indexl, ' n = ', nfreq,&
!  ' initPoint = ', initPoint, ' rate = ', actikrates%Lcool_fbE(4, act_proc)
ALLOCATE(freq(nfreq), cross(nfreq))
ALLOCATE(func(nfreq), intval(nfreq))
freq(:) = elements(indexe)%ions(indexi)%levels(indexl)%photcros(1,:)
cross(:) = elements(indexe)%ions(indexi)%levels(indexl)%photcros(2,:)
! physical informations in the propagation cell
! actual model grid index
cur_mgi = get_package_model_index(pack_index)
! electron density
el_dens = model_grid(cur_mgi)%e_dens
! temperature
temp = model_grid(cur_mgi)%T
! population
CALL populations(indexe, indexi, indexl, cur_mgi, act_pop)
! random number
ran_num = ran2(idum)
! ! calculation of initial frequency
init_freq = (MINVAL(elements(indexe)%ions(indexi + 1)%levels(:)%exci_energy) - &
       elements(indexe)%ions(indexi)%levels(indexl)%exci_energy) / h
! write(*,*) 'k_freq_fb: el = ', indexe, ' ion = ', indexi, ' lev = ', indexl, &
!  ' init_freq = ', init_freq, ' energy = ', init_freq * h / e_v
! STOP 'k_freq_fb: testing'
! getting the first point
!initPoint = 0
DO I = 1, nfreq
  ! write(*,*) 'cool_ionization: freq = ', init_freq, ' freq(', I, ') = ', freq(I)
 IF(init_freq < freq(I)) THEN
  initPoint = I
  EXIT
 END IF
END DO
! write(*,*) ' k_freq_fb: initPoint = ', initPoint, ' nfreq = ', nfreq
IF(initPoint == 0) STOP 'k_freq_fb: initPoint = 0'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculation of the integral
! going backwards to find frequency \nu which satisfies
! the equality ((4.37) (Kromer 2009))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculation of the first integral
! firstly filling the arrays
DO I = 1, nfreq
 x = ( h * freq(I)) / ( BOLK * temp )
 func(I) = cross(I) * freq(I)**3.0 * exp(-x)
END DO
DO I = 1, nfreq - 1
 intval(I) = (func(I + 1) + func(I)) * (freq(I + 1) - freq(I))
 ! write(*,*) 'k_freq_fb: intval(', I, ') = ', intval(I)
END DO
intval(nfreq) = 0.D0
! the integral calculation using the trapezoid rule
summ = 0.D0
DO I = initPoint, nfreq - 1
 summ = summ + intval(I)
END DO
! the left side of eq. 
integral1 = ran_num * summ
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculation of the second integral
summ = 0.D0
DO I = 1, nfreq - 1
 summ = summ + intval(nfreq - I)
 IF(summ > integral1) THEN
  actIndex = nfreq - I
  ! write(*,*) 'k_freq_fb: actIndex = ', actIndex, ' nfreq = ', nfreq, ' I = ', I
  EXIT
 END IF
END DO
! and linear interpolation
func1 = freq(actIndex)
func2 = freq(actIndex + 1)
int1 = summ
int2 = summ - intval(actIndex)
ali = (func1 - func2) / (int1 - int2)
bli = (func2 * int1 - func1 * int2) / (int1 - int2)
! finally we get random frequency
ran_freq = ali * integral1 + bli
! write(*,*) 'k_freq_fb: intval = ', intval(actIndex + 1)
! write(*,*) 'k_freq_fb: func1 = ', func1, ' func2 = ', func2, 'int1 = ', int1, &
!  ' int2 = ', int2, ' ali = ', ali, ' bli = ', bli
IF(ran_freq < 0.D0) THEN
 write(*,*) 'k_freq_fb: func1 = ', func1, ' func2 = ', func2, 'int1 = ', int1, &
  ' int2 = ', int2, ' ali = ', ali, ' bli = ', bli
 STOP
END IF
! write(*,*) 'k_freq_fb: ran_freq = ', ran_freq


END SUBROUTINE k_freq_fb
