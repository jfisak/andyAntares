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
! free-free variables
DOUBLE PRECISION, PARAMETER                     :: ffconst = 3.69255D8
DOUBLE PRECISION                                :: alphaff, gff
DOUBLE PRECISION                                :: kappaff

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
actirrates%Lcont(1, act_continuum) = 0
actirrates%Lcont(2, act_continuum) = 0 
actirrates%Lcont(3, act_continuum) = 0
actirrates%Lcont(4, act_continuum) = thomson
kappa = thomson
! write(*,*) 'r_kappa_cont: thomson = ', thomson
! write(*,*) 'r_kappa_cont: electron_density = ', electron_density, ' sigma_e = ', sigma_e
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
   IF(.NOT. ALLOCATED(elements(indexe)%ions(indexi)%levels(indexl)%photcros)) CYCLE
   n_sigma = SIZE(elements(indexe)%ions(indexi)%levels(indexl)%photcros(1,:))
   IF(n_sigma == 0) CYCLE
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
   ! no
   IF(actPoint == 0 .OR. actPoint == 1) THEN ! did we find the valid data?
    act_continuum = act_continuum + 1
    actirrates%Lcont(1, act_continuum) = indexe
    actirrates%Lcont(2, act_continuum) = indexi 
    actirrates%Lcont(3, act_continuum) = indexl
    actirrates%Lcont(4, act_continuum) = 0.D0
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
    CALL populations(indexe, indexi, indexl, current_mgi, act_pop)
    ! write(*,*)  'r_kappa_cont: cross_sect = ', cross_sect, ' act_continuum = ', act_continuum, &
    !  'act_pop = ', act_pop, ' 1-exp() = ', (1-exp(-(h * freq)/(BOLK * temp)))
    act_continuum = act_continuum + 1
    ! valid only for LTE approximation
    actirrates%Lcont(1, act_continuum) = indexe
    actirrates%Lcont(2, act_continuum) = indexi 
    actirrates%Lcont(3, act_continuum) = indexl
    actirrates%Lcont(4, act_continuum) = cross_sect * act_pop * (1-exp(-(h * freq)/(BOLK * temp)))
    !write(*,*) 'r_kappa_cont: act_continuum = ', act_continuum, &
    ! ' actirrates%Lcont(1, act_continuum) = ', actirrates%Lcont(1, act_continuum), &
    ! ' actirrates%Lcont(2, act_continuum) = ', actirrates%Lcont(2, act_continuum), &
    ! ' actirrates%Lcont(3, act_continuum) = ', actirrates%Lcont(3, act_continuum), &
    ! ' actirrates%Lcont(4, act_continuum) = ', actirrates%Lcont(4, act_continuum)
     kappa = kappa + actirrates%Lcont(4, act_continuum)
   END IF ! finding valid data
  END DO ! levels
 END DO ! ions
END DO ! elements

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Free-free opacity
act_continuum = act_continuum + 1
actirrates%Lcont(1, act_continuum) = 0
actirrates%Lcont(2, act_continuum) = 0
actirrates%Lcont(3, act_continuum) = 0
kappaff = 0.D0
DO indexe = 1, n_elements
 n_ions = SIZE(elements(indexe)%ions)
 DO indexi = 2, n_ions
  act_pop = model_grid(current_mgi)%grid_comp(indexe)%grid_ion(indexi)%tot_pop
  ! calculation of alpha_ff
  CALL gauntff(freq, temp, gff)
  alphaff = ffconst * DBLE((indexi - 1))**2 * gff /sqrt(temp) / freq ** 3.0 
  ! write(*,*) 'r_kappa_cont: ffconst = ', ffconst, ' i - 1 ^2 = ', DBLE((indexi - 1)**2), &
  !  ' gff = ', gff, ' sT = ', sqrt(temp), 'freq = ', freq
  kappaff = kappaff + electron_density * act_pop * alphaff * &
   (1.E0 - exp(-(h * freq) / (BOLK * temp)))
  ! write(*,*) 'r_kappa_cont: electron_density = ', electron_density, ' act_pop = ', &
  !  act_pop, ' alphaff = ', alphaff, ' 1-exp() = ', 1.E0 - exp(-(h * freq) / (BOLK * temp))
  ! write(*,*) 'r_kappa_cont: indexe = ', indexe, ' indexi = ', indexi, ' kappaff = ', kappaff
 END DO
END DO
kappa = kappa + kappaff
! write(*,*) 'r_kappa_cont: kappaff = ', kappaff, 'thomson = ', thomson
! write(*,*) 'r_kappa_cont: kappa = ', kappa
actirrates%Lcont(4, act_continuum) = kappaff
! STOP 'r_kappa_cont: testing'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Rayleigh scattering
! by helium




END SUBROUTINE r_kappa_cont
