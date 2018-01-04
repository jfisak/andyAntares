! calculation of absorption coefficient in continuum
SUBROUTINE r_kappa_cont(pack_index, kappa, actirrates)
USE types
USE rates_r

IMPLICIT NONE

! input variables
INTEGER                                         :: pack_index, n_cont
! 
INTEGER                                         :: current_mgi
DOUBLE PRECISION                                :: freq
! physical parameters
DOUBLE PRECISION                                :: electron_density
! element information
INTEGER                                         :: indexe, indexi, indexl
! loop variables
INTEGER                                         :: I
! actual frequency
DOUBLE PRECISION                                :: act_freq
! kappa coefficients
DOUBLE PRECISION                                :: thomson
! output variables
DOUBLE PRECISION                                :: kappa
! linear interpolation
DOUBLE PRECISION                                :: ali, bli, freq1, freq2, func1, func2
INTEGER                                         :: actPoint
DOUBLE PRECISION                                :: cross_sect
INTEGER                                         :: n_sigma
! ion information
INTEGER                                         :: n_ions, n_levels
DOUBLE PRECISION                                :: sigma_rs_he
INTEGER                                         :: get_package_model_index
INTEGER                                         :: act_continuum
DOUBLE PRECISION                                :: act_pop
INTEGER, PARAMETER                              :: n_thompson = 1
TYPE(rrates)                                    :: actirrates
DOUBLE PRECISION                                :: temp

!write(*,*) 'r_kappa_cont: dim(lcont) = ', SIZE(actirrates%Lcont)
!calculation of basic variables
current_mgi = get_package_model_index(pack_index)
electron_density = model_grid(current_mgi)%e_dens
temp = model_grid(current_mgi)%t
freq = package(pack_index)%freq_cmf
IF(current_mgi .EQ. n_modelgrid + 2) electron_density = 0.D0
!IF(.NOT. ALLOCATED(Lcont)) THEN
!END IF
!write(*,*) 'r_kappa_cont: n_photcrossect = ', n_photcrossect
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!  DIFFERENT OPACITY SOURCES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Thomson scattering
thomson = sigma_e * electron_density
act_continuum = 1
actirrates%Lcont(act_continuum) = thomson
kappa = thomson
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! bound-free processes
! loop over all possible cross sections
! computing number of cross sections
n_elements = SIZE(elements(:))
DO indexe = 1, n_elements
 n_ions = SIZE(elements(indexe)%ions)
 DO indexi = 1, n_ions - 1
  n_levels = SIZE(elements(indexe)%ions(indexi)%levels)
  DO indexl = 1, n_levels
   IF(.NOT. ALLOCATED(elements(indexe)%ions(indexi)%levels(indexl)%photcros)) CONTINUE
   n_sigma = SIZE(elements(indexe)%ions(indexi)%levels(indexl)%photcros(1,:))
!   n_sigma = 0
   actPoint = 0
   ! finding the propper index in saved photcross data
   DO I = 1, n_sigma
    act_freq = elements(indexe)%ions(indexi)%levels(indexl)%photcros(1,I)
    IF(freq < act_freq) THEN
     actPoint = I
     EXIT
    END IF
   END DO
   ! special cases
   ! frequency is lower than the first point
   IF(actPoint == 0 .OR. actPoint == 1) THEN ! did we find the valid data?
    ! no
    cross_sect = 0.D0
   ELSE
    ! yes
    ! now we can compute cross section from the data via linear interpolation
    freq1 = elements(indexe)%ions(indexi)%levels(indexl)%photcros(1, actPoint - 1)
    freq2 = elements(indexe)%ions(indexi)%levels(indexl)%photcros(1, actPoint)
    func1 = elements(indexe)%ions(indexi)%levels(indexl)%photcros(2, actPoint - 1)
    func2 = elements(indexe)%ions(indexi)%levels(indexl)%photcros(2, actPoint)
    ali = (func1 - func2) / (freq1 - freq2)
    bli = (func2 * freq1 - func1 * freq2) / (freq1 - freq2)
    cross_sect = ali * freq + bli
    !print*, 'r_kappa_cont: cross_sect = ', cross_sect, ' act_continuum = ', act_continuum!, ' dim(Lcont) = ', SIZE(actirrates%Lcont)
    CALL populations(indexe, indexi, indexl, current_mgi, act_pop)
    act_continuum = act_continuum + 1
    ! valid only for LTE approximation
    actirrates%Lcont(act_continuum) = cross_sect * act_pop * (1-exp(-(h * freq)/(BOLK * temp)))
    !actirrates%Lcont(act_continuum) = 0.D0
    write(*,*) 'r_kappa_cont: indexe = ', indexe, ' indexi = ', indexi, &
     ' indexl = ', indexl, ' actirrates%Lcont(', act_continuum, ') = ', actirrates%Lcont(act_continuum),&
     ' population = ', act_pop
    kappa = kappa + actirrates%Lcont(act_continuum)
   END IF ! finding valid data
  END DO ! levels
 END DO ! ions
END DO ! elements

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Rayleigh scattering
! by helium




END SUBROUTINE r_kappa_cont
