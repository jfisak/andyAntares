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
! initial point and frequency
INTEGER                                         :: initPoint
DOUBLE PRECISION                                :: init_freq
! integral calculations
DOUBLE PRECISION, ALLOCATABLE                   :: func(:), intval(:)
DOUBLE PRECISION                                :: summ
DOUBLE PRECISION                                :: integral1, integral2
INTEGER                                         :: actIndex
! propagation grid informations
INTEGER                                         :: cur_mgi, get_package_model_index
DOUBLE PRECISION                                :: el_dens, temp, x, act_pop
! random numbers
DOUBLE PRECISION                                :: random, ran_num
! loop variables
INTEGER                                         :: I
! linear interpolation
DOUBLE PRECISION                                :: ali, bli, int1, int2, func1, func2
INTEGER                                         :: actPoint
TYPE(krates)                                    :: actikrates

! informations about ion
indexe = actikrates%Lcool_fbind(1, act_proc)
indexi = actikrates%Lcool_fbind(2, act_proc)
indexl = actikrates%Lcool_fbind(3, act_proc)
! getting the photoionization cross section
! nfreq cannot be equal to zerou, because a process with a zero rate could
! not be chosen in the previous step
nfreq = SIZE(elements(indexe)%ions(indexi)%levels(indexl)%photcros(1,:))
ALLOCATE(freq(nfreq), cross(nfreq))
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
ran_num = random()
! calculation of initial frequency
init_freq = (elements(indexe)%ions(indexi + 1)%levels(1)%exci_energy - &
       elements(indexe)%ions(indexi)%levels(indexl)%exci_energy) / h
! getting the first point
DO I = 1, nfreq
 IF(init_freq < freq(I)) THEN
  initPoint = I
  EXIT
 END IF
END DO
ALLOCATE(func(nfreq - initPoint + 1), intval(nfreq - 1))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculation of the first integral
! firstly filling the arrays
DO I = initPoint, nfreq
 x = ( h * freq(I)) / ( BOLK * temp )
 func(I - initPoint + 1) = cross(I) * freq(I)**3.0 * exp(-x)
END DO
! the integral calculation using the trapezoid rule
summ = 0.D0
DO I = initPoint, nfreq - 1
 summ = (func(I + 1) + func(I)) * (freq(I + 1) - freq(I))
 intval(I) = summ
END DO
integral1 = ran_num * summ
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculation of the second integral
! going backwards to find frequency \nu which satisfies
! the equality ((4.37) (Kromer 2009))
summ = 0.D0
DO I = 1, nfreq - 1
 summ = summ + intval(nfreq - 1 - I)
 IF(summ > integral1) THEN
  actIndex = nfreq - 1 - I
  EXIT
 END IF
END DO
! and linear interpolation
int1 = freq(actIndex)
int2 = freq(actIndex + 1)
func1 = intval(actIndex)
func2 = intval(actIndex + 1)
ali = (func1 - func2) / (int1 - int2)
bli = (func2 * int1 - func1 * int2) / (int1 - int2)
! finally we get random frequency
ran_freq = ali * integral1 + bli
write(*,*) 'k_freq_fb: ran_freq = ', ran_freq


END SUBROUTINE k_freq_fb
