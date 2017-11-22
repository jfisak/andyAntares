SUBROUTINE i_radtrans(nlns, linetransitions, nluns, lineuptransitions, population, &
Zintdown, Zintup, Zrad)!, actirates)
USE types
USE rates
IMPLICIT NONE

! input variables
INTEGER                                 :: nlns, nluns
INTEGER, DIMENSION(nlns)                :: linetransitions
INTEGER, DIMENSION(nluns)               :: lineuptransitions
DOUBLE PRECISION                        :: stat_weight, exci_energy_u, exci_energy_l
DOUBLE PRECISION                        :: exci_energy
DOUBLE PRECISION                        :: population
INTEGER                                 :: element_index, ion_index
INTEGER                                 :: act_line
DOUBLE PRECISION                        :: actVal
INTEGER                                 :: I
! output variables
DOUBLE PRECISION                        :: Zintdown, Zintup, Zrad
INTEGER                         :: OMP_GET_THREAD_NUM, my_rank

my_rank = OMP_GET_THREAD_NUM()

Zintdown = 0.D0
Zrad= 0.D0
Zintup = 0.D0
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
ELSE
 Zintdown = 0.D0
 Zintup = 0.D0
 Zrad = 0.D0
END IF 
 !write(*,*) 'ALLOCATED: i_radtrans: my_rank = ', my_rank, ' recrad = ', size(actirates%Lma_recrad), &
 !           ' intrecrad = ', size(actirates%Lma_int_recrad), ' reccol = ', size(actirates%Lma_reccol), &
 !           ' intreccol = ', size(actirates%Lma_int_reccol)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! internal downward jump and radiative deexcitation
DO I = 1, nlns
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! 1.) internal downward jump
 act_line = linetransitions(I)
 ! the basic variables
 stat_weight = elements(element_index)%ions(ion_index)%levels(linelist(act_line)%upper)%stat_waight
 exci_energy_u = elements(element_index)%ions(ion_index)%levels(linelist(act_line)%upper)%exci_energy
 exci_energy_l = elements(element_index)%ions(ion_index)%levels(linelist(act_line)%lower)%exci_energy
 ! internal downward jump
 ! calculation of a rate coefficient
 actVal = population * linelist(act_line)%A_ul * exci_energy_l
 !print*, 'i_radtrans: stat_waight, exci_energy, population, linelist(act_line)%A_ul, actVal', &
 !      stat_weight, exci_energy, population, linelist(act_line)%A_ul, actVal
 actirates%Lma_int_dorad(I) = actVal * stat_weight
 Zintdown = Zintdown + actirates%Lma_int_dorad(I)
! write(*,*) 'i_radtrans: Zintdown = ', Zintdown
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! radiative deexcitation
 !print*, 'do_ipackage: population = ', population
 actVal = population * linelist(act_line)%A_ul * (exci_energy_u - exci_energy_l)
 Zrad = Zrad + stat_weight * actVal
END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! internal upward jump
DO I = 1, nluns
 act_line = lineuptransitions(I)
 exci_energy_l = elements(element_index)%ions(ion_index)%levels(linelist(act_line)%lower)%exci_energy
 stat_weight = elements(element_index)%ions(ion_index)%levels(linelist(act_line)%lower)%stat_waight
 exci_energy = elements(element_index)%ions(ion_index)%levels(linelist(act_line)%lower)%exci_energy
 ! internal jump up
 actVal = population * linelist(act_line)%A_ul * exci_energy_l * stat_weight
 actirates%Lma_int_uprad(I) = actVal
 Zintup = Zintup + actirates%Lma_int_uprad(I)
! write(*,*) 'i_radtrans: act_line = ', act_line, ' stat_weight = ', stat_weight, &
!  ' exci_energy_l = ', exci_energy_l, ' population = ', population, ' linelist(act_line)%A_ul = ', &
!  linelist(act_line)%A_ul
! write(*,*) 'i_radtrans: Zintup = ', Zintup
END DO

END SUBROUTINE i_radtrans
