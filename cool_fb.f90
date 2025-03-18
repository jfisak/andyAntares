!_________________________ k-packet rates ____________________________________________________
! calculation of fb cooling rates for k-packets
!
! Zfb -- the total cooling rate
!
! indexi -- index of an initial state
! indexi - 1 -- index of a final state
!
! the loops over ions are thus indexed 2, 3, nions
!
!_____________________________________________________________________________________________
SUBROUTINE cool_fb(pack_index, Zfb, actikrates)
USE types
USE constants
USE rates_k
IMPLICIT NONE

! input
INTEGER                                 :: pack_index
! output
DOUBLE PRECISION                        :: Zfb
! ion informations
INTEGER                                 :: indexe, indexi, indexl
INTEGER                                 :: n_ions, n_levels
! physical parameters
DOUBLE PRECISION                        :: el_dens, temp
INTEGER                                 :: cur_mgi
! photcross info
INTEGER                                 :: nfreq
DOUBLE PRECISION                        :: init_freq
DOUBLE PRECISION, ALLOCATABLE           :: crossfreq(:), cross(:), func(:), func2(:)
! index
INTEGER                                 :: ind_I
! loop variable
INTEGER                                 :: ind_J
INTEGER                                 :: act_rate
! integral calculation
DOUBLE PRECISION                        :: actInt, alphEspont, alphaSpont, summ, factor_x
INTEGER                                 :: actPoint
INTEGER                                 :: get_package_model_index
DOUBLE PRECISION                        :: tot_pop, uppper_en, lower_en
DOUBLE PRECISION                        :: sfactor
TYPE(krates)                            :: actikrates



! number of computed rates
act_rate = 0
! physical informations in the propagation cell
! actual model grid index
cur_mgi = get_package_model_index(pack_index)
! electron density
el_dens = model_grid(cur_mgi)%e_dens
! temperature
temp = model_grid(cur_mgi)%T

Zfb = 0.D0
DO indexe = 1, n_elements
 n_ions = SIZE(elements(indexe)%ions)
 DO indexi = 2, n_ions
  n_levels = SIZE(elements(indexe)%ions(indexi - 1)%levels)
  DO indexl = 1, n_levels
   IF(ALLOCATED(elements(indexe)%ions(indexi - 1)%levels(indexl)%photcros)) THEN
    nfreq = SIZE(elements(indexe)%ions(indexi - 1)%levels(indexl)%photcros(1,:))
   ELSE
    nfreq = 0
   END IF
   ! setting indexe and indexi
   act_rate = act_rate + 1
   actikrates%Lcool_fbE(ind_element, act_rate) = indexe
   actikrates%Lcool_fbE(ind_ion, act_rate) = indexi - 1
   actikrates%Lcool_fbE(ind_level, act_rate) = indexl
   IF(nfreq /= 0) THEN
    ALLOCATE(crossfreq(nfreq), cross(nfreq))
    crossfreq(:) = elements(indexe)%ions(indexi - 1)%levels(indexl)%photcros(1,:)
    cross(:) = elements(indexe)%ions(indexi - 1)%levels(indexl)%photcros(2,:)
    ! number density of a ground state of ion indexi + 1, indexe
    ! frequency
    init_freq = (MINVAL(elements(indexe)%ions(indexi)%levels(:)%exci_energy) - &
             elements(indexe)%ions(indexi - 1)%levels(indexl)%exci_energy) / const_h
    ! write(*,*) 'cool_fb: el = ', indexe, ' ion = ', indexi - 1, ' lev = ', indexl,&
    !  ' init_freq = ', init_freq, ' rate = ', actikrates%Lcool_fbE(act_rate)
    ! looking for initial point
    actPoint = 0
    DO ind_I = 1, nfreq
     IF(init_freq < crossfreq(ind_I)) THEN
      ! write(*,*) 'cool_ionization: freq = ', init_freq, ' crossfreq(', I, ') = ', crossfreq(I)
      actPoint = ind_I
      EXIT
     END IF
    END DO
    actikrates%Lcool_fbE(ind_initpoint, act_rate) = actPoint
    ! write(*,*) 'cool_ionization: el = ', indexe, ' ion = ', indexi - 1, ' lev = ', indexl,&
    !  ' n = ', nfreq, ' initp = ', actPoint
    ! if the initial point's frequency is too large we will not be able to
    ! calculate the integral, which is equal to zero in this case
    IF(actPoint == 0) THEN
     actikrates%Lcool_fbE(ind_rate,act_rate) = 0.D0
     ! write(*,*) 'cool_ionization: el = ', indexe, ' ion = ', indexi - 1, ' lev = ', indexl,&
     !  ' n = ', nfreq, ' initp = ', actPoint, ' rate = ', actikrates%Lcool_fbE(4, act_rate)
     DEALLOCATE(crossfreq, cross)
     CYCLE
    END IF
    ! calculation of the integral alpha E spont after Kromer(), Eq. (4.34)
    ! filling the arrays
    ALLOCATE(func(nfreq - actPoint + 1), func2(nfreq - actPoint + 1))
    DO ind_J = actPoint, nfreq
     ! write(*,*) 'cool_ionization: ind_J = ', ind_J, ' al(crfr) = ', ALLOCATED(crossfreq)
     factor_x = ( const_h * crossfreq(ind_J) ) / ( const_kB * temp)
     func(ind_J - actPoint + 1) = cross(ind_J) / (const_h * init_freq) * const_h * crossfreq(ind_J)**3.0 / const_c**2.0 * &
     exp(-factor_x)
     ! write(*,*) 'cool_ionization: cross = ', cross(ind_J), ' init_freq = ', init_freq
     func2(ind_J - actPoint + 1) = cross(ind_J) / (const_h * crossfreq(ind_J)) * const_h * crossfreq(ind_J)**3.0 / const_c**2.0 * &
     exp(-factor_x)
    END DO ! calculation of the integral
    ! calculation of integral using the trapezoid rule
    summ = 0.D0
    DO ind_J = 1, SIZE(func) - 1
     actInt = (func(ind_J) + func(ind_J + 1) ) * (crossfreq(ind_J + 1) - crossfreq(ind_J))
     summ = summ + actInt
     ! write(*,*) 'cool_fb: summ = ', summ
    END DO
    alphEspont = 4.D0 * const_pi * summ ! / const_c**2 / init_freq * summ
    ! write(*,*) 'cool_fb: summ = ', summ, ' alphEspont = ', alphEspont, ' init_freq = ', init_freq
    ! alpha spont after Kromer() eq. (4.35)
    summ = 0.D0
    DO ind_J = 1, SIZE(func) - 1
     actInt = (func2(ind_J)  + func2(ind_J + 1) ) * (crossfreq(ind_J + 1) - crossfreq(ind_J))
     summ = summ + actInt
    END DO
    alphaSpont = 4.D0 * const_pi * summ ! / const_c**2 * summ
    ! write(*,*) 'cool_ionization: alphaSpont = ', alphaSpont, 'alphEspont = ', alphEspont
    CALL saha_factor(indexe, indexi, indexl, temp,  sfactor)
    ! population of the given ion
    tot_pop = model_grid(cur_mgi)%grid_comp(indexe)%grid_ion(indexi)%tot_pop
    uppper_en = MINVAL(elements(indexe)%ions(indexi)%levels(:)%exci_energy)
    lower_en = elements(indexe)%ions(indexi - 1)%levels(indexl)%exci_energy
    actikrates%Lcool_fbE(ind_rate,act_rate) = tot_pop * el_dens * sfactor * &
     (alphEspont - alphaSpont) * (uppper_en - lower_en)
    Zfb = Zfb + actikrates%Lcool_fbE(ind_rate,act_rate)
    ! write(*,*) 'cool_fb: eldens = ', el_dens, ' tot_pop = ', tot_pop
    ! write(*,*) 'cool_fb: uppper_en = ', uppper_en, ' lower_en = ', lower_en
    ! write(*,*) 'cool_fb: sfactor = ', sfactor, ' ales = ', alphEspont, ' als = ',  alphaSpont 
    ! temporary solution
    ! Zfb = 0.D0
    ! write(*,*) 'cool_fb: Zfb = ', Zfb, 'actikrates%Lcool_fbE(', act_rate, ') = ', actikrates%Lcool_fbE(act_rate)
    DEALLOCATE(crossfreq, cross)
    DEALLOCATE(func, func2)
   ELSE
    actikrates%Lcool_fbE(ind_rate,act_rate) = 0.D0
   END IF ! npoints == 0
   !write(*,*) 'cool_fb: actikrates%Lcool_fbE(', act_rate, ') = ', actikrates%Lcool_fbE(act_rate), &
   ! ' Zfb = ', Zfb
   ! write(*,*) 'cool_ionization: el = ', indexe, ' ion = ', indexi - 1, ' lev = ', indexl,&
   !  ' n = ', nfreq, ' initp = ', actPoint, ' rate = ', actikrates%Lcool_fbE(4, act_rate)
  END DO ! levels
 END DO ! ions
END DO ! elements

!STOP 'cool_fb'

END SUBROUTINE cool_fb
