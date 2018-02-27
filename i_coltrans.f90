! this subroutine calculates collisional rates for the given energy
! level
SUBROUTINE i_coltrans(approx, pack_index, level, nlns, linetransitions, nluns, &
lineuptransitions, population, Zdown, Zup, Zcoll, actirates)
USE types
USE rates_i
IMPLICIT NONE
! input variables
! used approximation for the collisional term calculation
INTEGER                                 :: approx
! line -- number of line in the linelist field
INTEGER                                 :: pack_index, level, nlns, nluns
INTEGER, DIMENSION(nlns)                :: linetransitions
INTEGER, DIMENSION(nluns)               :: lineuptransitions
DOUBLE PRECISION                        :: population, low_pop
! constans
DOUBLE PRECISION, PARAMETER             :: c0 = 5.465D-11
DOUBLE PRECISION, PARAMETER             :: IH = 13.6 * e_v
DOUBLE PRECISION, PARAMETER             :: coll_const = 14.5
! indexes
INTEGER                                 :: act_line
! atomic data
INTEGER                                 :: element_index, ion_index
! value of collision rate
DOUBLE PRECISION                        :: actVal
INTEGER                                 :: current_mgi
! physical parameters of the cell
DOUBLE PRECISION                        :: el_temperature, electron_density, &
                                           temperature
INTEGER                                 :: get_package_model_index
! loop variable
INTEGER                                 :: I
! a value of gamma function
DOUBLE PRECISION                        :: gf
! parameters of transitions
DOUBLE PRECISION                        :: freq, osc_str
DOUBLE PRECISION                        :: exci_energy_l, exci_energy_u
DOUBLE PRECISION                        :: x
! output variables
! total rate
DOUBLE PRECISION                        :: Zdown, Zup, Zcoll
! the rates for the given transitions
DOUBLE PRECISION                        :: lcoll
DOUBLE PRECISION, PARAMETER             :: times = 1.D0
TYPE(irates)                           :: actirates

!IF(times /= 1.D0) THEN
! CALL warning('collisional rates are multiplied by a non-one factor')
!END IF
SELECT CASE(approx)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! van Regemorter approximation
CASE(1)
 current_mgi = get_package_model_index(pack_index)
 ! electron density
 electron_density = model_grid(current_mgi)%e_dens
 ! --
 temperature = model_grid(current_mgi)%T
 el_temperature = temperature
 ! number of lines we are interested in
 nlns = SIZE(linetransitions)

 ! we have to know which element and ion we are calculating data for
 ! we will use the knowledge of lines and assume that at it is
 ! possible at least one transition upwards or downwards and from
 ! the first element we get the element and the ion informations
 IF(SIZE(linetransitions) /= 0) THEN
  element_index = linelist(linetransitions(1))%indexe
  ion_index = linelist(linetransitions(1))%indexi
 ELSE IF(SIZE(lineuptransitions) /= 0) THEN
  element_index = linelist(lineuptransitions(1))%indexe
  ion_index = linelist(lineuptransitions(1))%indexi
 END IF 
 ! initialization of rate values
 Zup = 0.D0
 Zdown = 0.D0
 Zcoll = 0.D0
 DO I = 1, nlns
  ! important physical quantities
  act_line = linetransitions(I)
  osc_str = linelist(act_line)%f_ul
  exci_energy_l = elements(element_index)%ions(ion_index)%levels(linelist(act_line)%lower)%exci_energy
  exci_energy_u = elements(element_index)%ions(ion_index)%levels(linelist(act_line)%upper)%exci_energy
  ! frequency of transition
  freq = linelist(act_line)%freq
  x = (h * freq) / (BOLK * temperature)
  ! gamma function
  CALL gamma_function(x, act_line, gf)
  ! value of the collision coefficient c_{i, j, k -> i', j, k}
  actVal = electron_density * c0 * (temperature)**(1.0/2.0) * &
   coll_const * (IH / (h * freq)) * osc_str * ((h * freq) / (BOLK * el_temperature)) * &
   exp(-(h * freq) / (BOLK * el_temperature)) * gf
!  actVal = actVal * (linelist(act_line)%upper - linelist(act_line)%lower)
 ! internal downward jump
  CALL populations(element_index, ion_index, linelist(act_line)%lower, current_mgi, low_pop)
  actirates%Lma_int_docoll(I) = low_pop * actVal * exci_energy_l
  Zdown = Zdown + actirates%Lma_int_docoll(I)
 ! collisional deexcitation
  lcoll = low_pop * actVal * (exci_energy_u - exci_energy_l)
  ! write(*,*) 'i_coltrans: lcoll = ', lcoll, ' low_pop = ', low_pop
  Zcoll = Zcoll + lcoll
!        (linelist(act_line)%upper - linelist(act_line)%lower)
 ! write(*,*) 'i_coltrans: low_pop = ', low_pop, ' el_temperature = ', el_temperature, &
 !  ' electron_density = ', electron_density, ' temperature = ', temperature
 END DO
 ! print*, 'collisional_rates: lcoll = ', lcoll
! Ztot = pop * Ztot
! print*, 'collisional_rates: Ztot = ', Ztot
 ! upward jumps
 DO I = 1, nluns
  act_line = lineuptransitions(I)
  ! important physical quantities
  osc_str = linelist(act_line)%f_ul
  exci_energy_l = elements(element_index)%ions(ion_index)%levels(linelist(act_line)%lower)%exci_energy
  ! frequency of transition
  freq = linelist(act_line)%freq
  x = (h * freq) / (BOLK * temperature)
  ! gamma function
  CALL gamma_function(x, act_line, gf)
  ! value of the collision coefficient
  actVal = population * electron_density * c0 * (temperature)**(1.0/2.0) * &
        coll_const * (IH / (h * freq)) * osc_str * &
        ((h * freq) / (BOLK * el_temperature)) * &
        exp(-(h * freq) / (BOLK * el_temperature)) * gf
  actirates%Lma_int_upcoll(I) = actVal * exci_energy_l
  ! write(*,*) 'i_coltrans: Lma_int_upcoll = ', actirates%Lma_int_upcoll(I), ' population = ', population
 END DO
CASE DEFAULT
 STOP 'collisional_rates: this approximation is not known'
END SELECT


END SUBROUTINE i_coltrans
