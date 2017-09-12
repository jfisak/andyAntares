SUBROUTINE cool_fb(pack_index, Zfb)
USE types
USE rates
IMPLICIT NONE

! input
INTEGER                                 :: pack_index
! output
DOUBLE PRECISION                        :: Zfb
! ion informations
INTEGER                                 :: indexe, indexi, indexl
INTEGER                                 :: n_ions, n_levels
INTEGER                                 :: n_phcs
! physical parameters
DOUBLE PRECISION                        :: el_dens, temp
INTEGER                                 :: cur_mgi
! photcross info
INTEGER                                 :: nfreq
DOUBLE PRECISION                        :: init_freq
DOUBLE PRECISION, ALLOCATABLE           :: crossfreq(:), cross(:), func(:)
! index
INTEGER                                 :: I
! loop variable
INTEGER                                 :: J
INTEGER                                 :: act_rate
! integral calculation
DOUBLE PRECISION                        :: actInt, integral, summ, x
INTEGER                                 :: actPoint
INTEGER                                 :: get_package_model_index
DOUBLE PRECISION                        :: tot_pop, uppper_en, lower_en



IF(.NOT. ASSOCIATED(Lcool_fbE)) THEN
 n_phcs = 0
 DO indexe = 1, n_elements
  n_ions = SIZE(elements(indexe)%ions)
  DO indexi = 2, n_ions 
   n_levels = SIZE(elements(indexe)%ions(indexi)%levels)
   DO indexl = 1, n_levels
    n_phcs = n_phcs + 1
   END DO ! levels
  END DO ! ions
 END DO ! elements
 ALLOCATE(Lcool_fbE(n_phcs), Lcool_fbind(2,n_phcs))
END IF

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
  n_levels = SIZE(elements(indexe)%ions(indexi)%levels)
  DO indexl = 1, n_levels
   ! for this case we have an ion without electrons thus photoionization is impossible
   IF(indexi == indexe + 1) THEN
    nfreq = 0
   ELSE
    nfreq = SIZE(elements(indexe)%ions(indexi)%levels(indexl)%photcros(1,:))
   END IF
   ! setting indexe and indexi
   act_rate = act_rate + 1
   Lcool_fbind(1, act_rate) = indexe
   Lcool_fbind(2, act_rate) = indexi
   IF(nfreq /= 0) THEN
    ALLOCATE(crossfreq(nfreq), cross(nfreq))
    crossfreq(:) = elements(indexe)%ions(indexi)%levels(indexl)%photcros(1,:)
    cross(:) = elements(indexe)%ions(indexi)%levels(indexl)%photcros(2,:)
    ! number density of a ground state of ion indexi + 1, indexe
    ! frequency
    init_freq = (elements(indexe)%ions(indexi + 1)%levels(1)%exci_energy - &
             elements(indexe)%ions(indexi)%levels(indexl)%exci_energy) / h
    ! looking for initial point
    DO I = 1, nfreq
     IF(init_freq < crossfreq(I)) THEN
      !write(*,*) 'cool_ionization: freq = ', freq, ' crossfreq(', I, ') = ', crossfreq(I)
      actPoint = I
      EXIT
     END IF
    END DO
    ! if the initial point's frequency is too large we will not be able to
    ! calculate the integral, which is equal to zero in this case
    IF(actPoint == 0) THEN
     Lcool_fbE(act_rate) = 0.D0
     CYCLE
    END IF
    ! calculation of the integral
    ! filling the arrays
    ALLOCATE(func(nfreq - actPoint + 1))
    DO J = actPoint, nfreq
     x = ( h * crossfreq(J) ) / ( BOLK * temp)
     func(J - actPoint + 1) = cross(J) * crossfreq(J)**3.0 * exp(-x)
    END DO ! calculation of the integral
    ! calculation of integral using the trapezoid rule
    summ = 0.D0
    DO J = 1, SIZE(func) - 1
     actInt = (func(J) * cross(J) + func(J + 1) * cross(J + 1)) * &
      (crossfreq(J + 1) - crossfreq(J))
     summ = summ + actInt
    END DO
    integral = 4.D0 * pi / light_speed**2 / init_freq * summ
    ! population of the given ion
    tot_pop = model_grid(cur_mgi)%grid_comp(indexe)%grid_ion(indexi)%tot_pop
    uppper_en = elements(indexe)%ions(indexi)%levels(1)%exci_energy
    lower_en = elements(indexe)%ions(indexi - 1)%levels(indexl)%exci_energy
    Lcool_fbE(act_rate) = tot_pop * integral * (uppper_en - lower_en)
    Zfb = Zfb + Lcool_fbE(act_rate)
    !write(*,*) 'cool_fb: Zfb = ', Zfb, 'Lcool_fbE(', act_rate, ') = ', Lcool_fbE(act_rate)
    DEALLOCATE(func, crossfreq, cross)
   ELSE
    Lcool_fbE(act_rate) = 0.D0
   END IF ! npoints == 0
   !write(*,*) 'cool_fb: Lcool_fbE(', act_rate, ') = ', Lcool_fbE(act_rate), &
   ! ' Zfb = ', Zfb
  END DO ! levels
 END DO ! ions
END DO ! elements

!STOP 'cool_fb'

END SUBROUTINE cool_fb
