!______________________ i-packet ionization and recombination rates _____________________
! precalculation of integrals for ionization and recombination rates
! 
! for a temperature grid with temperatures 
SUBROUTINE i_ion_recomb(division_type)

USE types
USE constants
USE rates_i
IMPLICIT NONE

INTEGER                                         :: division_type
DOUBLE PRECISION                                :: Tmin, Tmax
INTEGER                                         :: ntints
! cross section data
DOUBLE PRECISION, ALLOCATABLE                   :: freq(:), cross(:)  
DOUBLE PRECISION, ALLOCATABLE                   :: func1(:), func2(:)!, func3(:)
DOUBLE PRECISION                                :: freqt
! loop index
INTEGER                                         :: I, act_int
INTEGER                                         :: Istart
DOUBLE PRECISION                                :: cur_temp
INTEGER                                         :: index_temp
INTEGER                                         :: indexe, indexi, indexl
INTEGER                                         :: n_ions, n_levels, npoints
DOUBLE PRECISION                                :: flux, flux_function
! x = h\nu / (k_B * T)
DOUBLE PRECISION                                :: x
DOUBLE PRECISION                                :: sfactor
DOUBLE PRECISION                                :: summ1, summ2! , summ3

INTEGER                                         :: nlns, nluns
INTEGER                                         :: actual_state
INTEGER                                         :: J, K


! temperature grid
! starting and ending point
Tmin = MINVAL(model_grid(1:n_modelgrid)%t)
Tmax = MAXVAL(model_grid(1:n_modelgrid)%t)

! write(*,*) 'i_ion_recomb: Tmin = ', Tmin, ' Tmax = ', Tmax
! STOP 'i_ion_recomb: testing'
! write(*,*) 'i_ion_recomb: t = ', model_grid(:)%t
IF(Tmax == Tmin) THEN
 ntints = 1
ELSE
ntints = INT((Tmax - Tmin)/200)
END IF
ALLOCATE(i_temps(ntints))
ALLOCATE(iints(n_photcrossect))
DO I = 1, n_photcrossect
 ALLOCATE(iints(I)%gammaijk(ntints))
 ALLOCATE(iints(I)%alphaijk(ntints))
 ! ALLOCATE(iints(I)%etaijk(ntints))
END DO
! write(*,*) 'i_ion_recomb: ntints = ', ntints

SELECT CASE(division_type)
! linear division of intervals
CASE(1)
 IF(ntints > 1) THEN
  DO I = 1, ntints
   i_temps(I) = ( Tmax - Tmin ) / DBLE(ntints - 1) * DBLE(I) + (DBLE(ntints) * Tmin - Tmax) / (DBLE(ntints - 1))
  END DO
 ELSE IF(ntints == 1) THEN
  i_temps(1) = Tmin
 END IF
! logarithmic division of intervals
CASE(2)
CASE DEFAULT
 STOP 'i_ion_recomb: this choice is not known'
END SELECT

act_int = 0
! 1.) ionization rates
DO indexe = 1, n_elements
 n_ions = SIZE(elements(indexe)%ions)
 DO indexi = 1, n_ions - 1
  n_levels = SIZE(elements(indexe)%ions(indexi)%levels)
  DO indexl = 1, n_levels
   ! number of rows in the photcross data
   IF(ALLOCATED(elements(indexe)%ions(indexi)%levels(indexl)%photcros)) THEN
    npoints = SIZE(elements(indexe)%ions(indexi)%levels(indexl)%photcros(1,:))
   ELSE
    CYCLE
   END IF
   IF(npoints == 0) CYCLE
   ! write(*,*) 'i_ion_recomb: npoints = ', npoints
   ! func1 -- \gamma_{i, j, k}
   ! func2 -- \alpha^{spont.}_{i, j, k}
   ALLOCATE(freq(npoints), cross(npoints))
   act_int = act_int + 1
   freq(1:npoints) = elements(indexe)%ions(indexi)%levels(indexl)%photcros(1,1:npoints)
   cross(1:npoints) = elements(indexe)%ions(indexi)%levels(indexl)%photcros(2,1:npoints)
   freqt = (MINVAL(elements(indexe)%ions(indexi + 1)%levels(:)%exci_energy) - &
    elements(indexe)%ions(indexi)%levels(indexl)%exci_energy) / h
   ! looking for starting point
   Istart = 0
   DO I = 1, npoints
    IF(freq(I) >= freqt) THEN
     Istart = I
     ALLOCATE(func1(npoints - Istart + 1))
     ALLOCATE(func2(npoints - Istart + 1))
     ! ALLOCATE(func3(npoints - Istart + 1))
     EXIT
    END IF
   END DO
   IF(Istart == 0) THEN
    Istart = npoints
   END IF
   ! write(*,*) 'i_ion_recomb: Istart = ', Istart, ' npoints = ', npoints
   ! the integral calculation
   DO index_temp = 1, ntints
    cur_temp = i_temps(index_temp)
    DO I = Istart, npoints 
     flux = flux_function(0,freq(I), cur_temp, R_star)
     x = h * freq(I) / (BOLK * cur_temp)
    ! func(I) = cross(I) * flux / ( h * freq(I))
     func1(I - Istart + 1) = flux * cross(I) / (h * freq(I)) * (1.0 - exp(-x))
     func2(I - Istart + 1) = cross(I) / (h * freq(I)) * &
      (2.0 * h * freq(I)**3.0) / light_speed**2.0 * exp(-x)
     ! func3(I - Istart + 1) = cross(I) * &
     !  (2.0 * h * freq(I)**3.0) / light_speed**2.0 * exp(-x)
    ! write(*,*) 'i_ion_recomb: J = ', I - Istart + 1
    ! write(*,*)  'flux = ', flux, ' cross(I) = ', cross(I), ' func1(J) = ', func1(I - Istart + 1), &
    !  ' func2(J) = ', func2(I - Istart + 1), ' h * freq = ', h * freq(I)
    END DO ! integral calculation
    summ1 = 0.D0
    summ2 = 0.D0
    ! summ3 = 0.D0
    ! write(*,*) 'i_ion_recomb: S(func) = ', SIZE(func1), ' dim = ', npoints - Istart + 1
    DO I = 1, npoints - Istart
     summ1 = summ1 + (func1(I) + func1(I + 1)) / 2.D0 * (freq(I + 1) - freq(I))
     summ2 = summ2 + (func2(I) + func2(I + 1)) / 2.D0 * (freq(I + 1) - freq(I))
     ! summ3 = summ3 + (func3(I) + func3(I + 1)) / 2.D0 * (freq(I + 1) - freq(I))
    END DO
    ! write(*,*) 'i_ion_recomb: temp = ', cur_temp
    CALL saha_factor(indexe, indexi + 1, indexl, cur_temp, sfactor)
    iints(act_int)%indexe = indexe
    iints(act_int)%indexi = indexi
    iints(act_int)%indexl = indexl
    iints(act_int)%gammaijk(index_temp) = 4.D0 * pi * summ1 
    iints(act_int)%alphaijk(index_temp) = 4.D0 * pi * summ2 * sfactor
    ! iints(act_int)%etaijk(index_temp) = summ3
    ! write(*,*) 'i_ion_recomb: gamma = ', iints(act_int)%gammaijk(index_temp)
    ! write(*,*) 'i_ion_recomb: act_int = ', act_int, ' index_temp = ', index_temp, ' alpha = ', iints(act_int)%alphaijk(index_temp)
    ! write(*,*) 'i_ion_recomb: T = ', cur_temp, ' gamma = ', int1, ' phc = ', int2
    ! write(*,*) 'i_ion_recomb: T = ', cur_temp, ' summ1 = ', summ1, ' summ2 = ', summ2
    ! write(77,*) cur_temp, summ1 * sfactor, summ2 * sfactor
   END DO ! loop over temperatures
   elements(indexe)%ions(indexi)%levels(indexl)%phfreqi = act_int
   ! write(*,*) 'i_ion_recomb: act_int = ', act_int
   DEALLOCATE(freq, cross, func1, func2)
  ! DEALLOCATE(freq, cross, func1, func2, func3)
  END DO ! loop over ionic levels
 END DO ! loop over ions
END DO ! loop over elements

! calculation of possible upper and lower transitions for every energy level
DO indexe = 1, n_elements
 n_ions = SIZE(elements(indexe)%ions)
 DO indexi = 1, n_ions
  n_levels = SIZE(elements(indexe)%ions(indexi)%levels)
  DO actual_state = 1, n_levels
   ! number of transitions up and down
   nlns = 0
   nluns = 0
    DO I = 1, ntransitions
     IF(linelist(I)%indexe == indexe .AND. linelist(I)%indexi == indexi) THEN
      ! transitions to a lower level
      IF(linelist(I)%upper == actual_state) THEN
       nlns = nlns + 1
      END IF
      ! transitions to a upper level
      IF(linelist(I)%lower == actual_state) THEN
       nluns = nluns + 1
      END IF
     END IF
    ! write(*,*) 'do_ipackage: ', linelist(I)%indexe, linelist(I)%indexi, linelist(I)%lower, linelist(I)%upper
    END DO
    ! STOP 'do_ipackage: testing'
    ! write(*,*) 'do_ipackage: number of found transitions: ', nlns, nluns
    ALLOCATE(elements(indexe)%ions(indexi)%levels(actual_state)%linetransitions(nlns))
    ALLOCATE(elements(indexe)%ions(indexi)%levels(actual_state)%lineuptransitions(nluns))
    ! we will save these possible transitions into an array
    J = 0
    K = 0
    DO I = 1, ntransitions
     IF(linelist(I)%indexe == indexe .AND. linelist(I)%indexi == indexi) THEN
      ! transitions to a lower level
      IF(linelist(I)%upper == actual_state) THEN
       J = J + 1
       elements(indexe)%ions(indexi)%levels(actual_state)%linetransitions(J) = I
      END IF
      ! transitions to a upper level
      IF(linelist(I)%lower == actual_state) THEN
       K = K + 1
       elements(indexe)%ions(indexi)%levels(actual_state)%lineuptransitions(K) = I
      END IF
     END IF
    END DO
   ! write(*,*) 'i_ion_recomb: el = ', indexe, ' ion = ', indexi, ' l = ', actual_state, ' nlns = ', nlns, &
   ! ' linetransitions = ', elements(indexe)%ions(indexi)%levels(actual_state)%linetransitions
  END DO ! loop over levels
 END DO ! loop over ions
END DO ! loop over elements

! STOP 'i_ion_recomb: testing'


END SUBROUTINE i_ion_recomb
