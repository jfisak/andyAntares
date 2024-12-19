! generates random frequency for a free-bound transition
!
! INPUT: pack_index(INT): the index of the packet
!        act_proc(INT): index of the current procedure
!        actikrates(krates): the rates for the corresponding transitions
! OUTPUT: ran_freq(DBLE): a random frequency
!
SUBROUTINE k_freq_fb(pack_index, act_proc, ran_freq, actikrates)
USE types
USE constants
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
INTEGER                                         :: ind_I
! linear interpolation
DOUBLE PRECISION                                :: ali, bli, int1, int2, func1, func2
TYPE(krates)                                    :: actikrates

! informations about ion
indexe = INT(actikrates%Lcool_fbE(ind_element, act_proc))
indexi = INT(actikrates%Lcool_fbE(ind_ion, act_proc))
indexl = INT(actikrates%Lcool_fbE(ind_level, act_proc))
initPoint = INT(actikrates%Lcool_fbE(ind_initpoint, act_proc))
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
       elements(indexe)%ions(indexi)%levels(indexl)%exci_energy) / const_h
! write(*,*) 'k_freq_fb: el = ', indexe, ' ion = ', indexi, ' lev = ', indexl, &
!  ' init_freq = ', init_freq, ' energy = ', init_freq * const_h / e_v
! STOP 'k_freq_fb: testing'
! getting the first point
!initPoint = 0
DO ind_I = 1, nfreq
  ! write(*,*) 'cool_ionization: freq = ', init_freq, ' freq(', ind_I, ') = ', freq(ind_I)
 IF(init_freq < freq(ind_I)) THEN
  initPoint = ind_I
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
DO ind_I = 1, nfreq
 x = ( const_h * freq(ind_I)) / ( const_kB * temp )
 func(ind_I) = cross(ind_I) * freq(ind_I)**3.0 * exp(-x)
END DO
DO ind_I = 1, nfreq - 1
 intval(ind_I) = (func(ind_I + 1) + func(ind_I)) * (freq(ind_I + 1) - freq(ind_I))
 ! write(*,*) 'k_freq_fb: intval(', ind_I, ') = ', intval(ind_I)
END DO
intval(nfreq) = 0.D0
! the integral calculation using the trapezoid rule
summ = 0.D0
DO ind_I = initPoint, nfreq - 1
 summ = summ + intval(ind_I)
END DO
! the left side of eq. 
integral1 = ran_num * summ
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculation of the second integral
summ = 0.D0
DO ind_I = 1, nfreq - 1
 summ = summ + intval(nfreq - ind_I)
 IF(summ > integral1) THEN
  actIndex = nfreq - ind_I
  ! write(*,*) 'k_freq_fb: actIndex = ', actIndex, ' nfreq = ', nfreq, ' ind_I = ', ind_I
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


END SUBROUTINE k_freq_fb
