SUBROUTINE cool_fb(pack_index, Zfb, actikrates)
USE types
USE rates_k
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
TYPE(krates)                            :: actikrates


! indexi - 1 -- index of an initial state
! indexi -- index of a final state

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
   actikrates%Lcool_fbE(1, act_rate) = indexe
   actikrates%Lcool_fbE(2, act_rate) = indexi - 1
   actikrates%Lcool_fbE(3, act_rate) = indexl
   IF(nfreq /= 0) THEN
    ALLOCATE(crossfreq(nfreq), cross(nfreq))
    crossfreq(:) = elements(indexe)%ions(indexi - 1)%levels(indexl)%photcros(1,:)
    cross(:) = elements(indexe)%ions(indexi - 1)%levels(indexl)%photcros(2,:)
    ! number density of a ground state of ion indexi + 1, indexe
    ! frequency
    init_freq = (MINVAL(elements(indexe)%ions(indexi)%levels(:)%exci_energy) - &
             elements(indexe)%ions(indexi - 1)%levels(indexl)%exci_energy) / h
    ! write(*,*) 'cool_fb: el = ', indexe, ' ion = ', indexi - 1, ' lev = ', indexl,&
    !  ' init_freq = ', init_freq, ' rate = ', actikrates%Lcool_fbE(act_rate)
    ! looking for initial point
    DO I = 1, nfreq
     IF(init_freq < crossfreq(I)) THEN
      ! write(*,*) 'cool_ionization: freq = ', init_freq, ' crossfreq(', I, ') = ', crossfreq(I)
      actPoint = I
      EXIT
     END IF
    END DO
    actikrates%Lcool_fbE(5, act_rate) = actPoint
    ! write(*,*) 'cool_ionization: el = ', indexe, ' ion = ', indexi - 1, ' lev = ', indexl,&
    !  ' n = ', nfreq, ' initp = ', actPoint
    ! if the initial point's frequency is too large we will not be able to
    ! calculate the integral, which is equal to zero in this case
    IF(actPoint == 0) THEN
     actikrates%Lcool_fbE(4,act_rate) = 0.D0
     ! write(*,*) 'cool_ionization: el = ', indexe, ' ion = ', indexi - 1, ' lev = ', indexl,&
     !  ' n = ', nfreq, ' initp = ', actPoint, ' rate = ', actikrates%Lcool_fbE(4, act_rate)
     CYCLE
    END IF
    ! calculation of the integral
    ! filling the arrays
!    ALLOCATE(func(nfreq - actPoint + 1))
!    DO J = actPoint, nfreq
!     x = ( h * crossfreq(J) ) / ( BOLK * temp)
!     func(J - actPoint + 1) = cross(J) * crossfreq(J)**3.0 * exp(-x)
!    END DO ! calculation of the integral
!    ! calculation of integral using the trapezoid rule
!    summ = 0.D0
!    DO J = 1, SIZE(func) - 1
!     actInt = (func(J) * cross(J) + func(J + 1) * cross(J + 1)) * &
!      (crossfreq(J + 1) - crossfreq(J))
!     summ = summ + actInt
!    END DO
!    integral = 4.D0 * pi / light_speed**2 / init_freq * summ
!    ! population of the given ion
    tot_pop = model_grid(cur_mgi)%grid_comp(indexe)%grid_ion(indexi)%tot_pop
    uppper_en = MINVAL(elements(indexe)%ions(indexi)%levels(:)%exci_energy)
    lower_en = elements(indexe)%ions(indexi - 1)%levels(indexl)%exci_energy
    ! actikrates%Lcool_fbE(4,act_rate) = tot_pop * integral * (uppper_en - lower_en)
    actikrates%Lcool_fbE(4,act_rate) = tot_pop * el_dens * (uppper_en - lower_en)
    Zfb = Zfb + actikrates%Lcool_fbE(4,act_rate)
    ! temporary solution
    ! Zfb = 0.D0
    !write(*,*) 'cool_fb: Zfb = ', Zfb, 'actikrates%Lcool_fbE(', act_rate, ') = ', actikrates%Lcool_fbE(act_rate)
    DEALLOCATE(crossfreq, cross)
   ELSE
    actikrates%Lcool_fbE(4,act_rate) = 0.D0
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
