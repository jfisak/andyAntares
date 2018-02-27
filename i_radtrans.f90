!_________________________ i-packet rates ____________________________________________________
! calculation of radiative transitions rates of an i-packets
!
! * Zintdown -- total rate of internal downward jumps
! * Zrad -- total rate of radiative deactivations of a macro atom
! * Zintup -- total rate of internal upward jumps
SUBROUTINE i_radtrans(current_mgi, nlns, linetransitions, nluns, lineuptransitions, up_pop, &
Zintdown, Zintup, Zrad, actirates)
USE types
USE rates_i
IMPLICIT NONE

! input variables
INTEGER                                 :: nlns, nluns
INTEGER                                 :: current_mgi
INTEGER, DIMENSION(nlns)                :: linetransitions
INTEGER, DIMENSION(nluns)               :: lineuptransitions
DOUBLE PRECISION                        :: stat_weight_l, stat_weight_u
DOUBLE PRECISION                        :: exci_energy_u, exci_energy_l
DOUBLE PRECISION                        :: up_pop, low_pop
! beta lu calculation
DOUBLE PRECISION                        :: taulu, betalu
DOUBLE PRECISION                        :: Bul, Blu, Jlu
DOUBLE PRECISION                        :: flux_function
INTEGER                                 :: element_index, ion_index
INTEGER                                 :: act_line
DOUBLE PRECISION                        :: actVal
INTEGER                                 :: I
! output variables
DOUBLE PRECISION                        :: Zintdown, Zintup, Zrad
INTEGER                         :: OMP_GET_THREAD_NUM, my_rank
TYPE(irates)      :: actirates

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
 !       ' intdorad = ', size(actirates%Lma_int_dorad)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! internal downward jump and radiative deexcitation
DO I = 1, nlns
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! 1.) internal downward jump
 act_line = linetransitions(I)
 ! the basic variables
 stat_weight_u = elements(element_index)%ions(ion_index)%levels(linelist(act_line)%upper)%stat_waight
 stat_weight_l = elements(element_index)%ions(ion_index)%levels(linelist(act_line)%lower)%stat_waight
 exci_energy_u = elements(element_index)%ions(ion_index)%levels(linelist(act_line)%upper)%exci_energy
 exci_energy_l = elements(element_index)%ions(ion_index)%levels(linelist(act_line)%lower)%exci_energy
 ! internal downward jump
 ! calculation of a rate coefficient
 CALL populations(element_index, ion_index, I, current_mgi, low_pop)
 ! Einstein Blu coefficient
 Blu = 4.0 * pi / (h * linelist(act_line)%freq) * linelist(act_line)%A_ul
 ! optical depth
 taulu = low_pop * Blu * h * light_speed / (4.0 * pi) *&
  (1.D0 - (stat_weight_l * up_pop) / (stat_weight_u * up_pop))
 ! probability of escape of the packet after scattering in line
 betalu = 1 / taulu * (1 - exp(-taulu))
 actVal = up_pop * betalu * Blu
 ! write(*,*) 'i_radtrans: exci_energy, up_pop, Blu, actVal', exci_energy_u, up_pop, Blu, actVal
 actirates%Lma_int_dorad(I) = actVal * exci_energy_l
 IF(actirates%Lma_int_dorad(I) < 0.D0) STOP 'i_radtrans: Lma_int_dorad < 0'
 Zintdown = Zintdown + actirates%Lma_int_dorad(I)
 ! write(*,*) 'i_radtrans: Zintdown = ', Zintdown
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! radiative deexcitation
 !print*, 'do_ipackage: up_pop = ', up_pop
 actVal = up_pop * betalu * linelist(act_line)%A_ul
 Zrad = Zrad + actVal * (exci_energy_u - exci_energy_l)
 IF(exci_energy_u - exci_energy_l < 0) STOP 'i_radtrans: exci_energy_u - exci_energy_l < 0'
 ! write(*,*) 'i_radtrans: Zrad = ', Zrad
END DO
 ! STOP 'i_radtrans, testing'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! internal upward jump
low_pop = up_pop
DO I = 1, nluns
 act_line = lineuptransitions(I)
 exci_energy_l = elements(element_index)%ions(ion_index)%levels(linelist(act_line)%lower)%exci_energy
 exci_energy_u = elements(element_index)%ions(ion_index)%levels(linelist(act_line)%upper)%exci_energy
 stat_weight_l = elements(element_index)%ions(ion_index)%levels(linelist(act_line)%lower)%stat_waight
 stat_weight_u = elements(element_index)%ions(ion_index)%levels(linelist(act_line)%upper)%stat_waight
! CALL populations(element_index, ion_index, linelist(act_line)%lower, current_mgi, low_pop)
 CALL populations(element_index, ion_index, linelist(act_line)%upper, current_mgi, up_pop)
 Jlu = flux_function(0, linelist(act_line)%freq, model_grid(current_mgi)%T)
 ! calculation of Blu and Bul
 Blu = 4.0 * pi / (h * linelist(act_line)%freq) * linelist(act_line)%A_ul
 Bul = stat_weight_l / stat_weight_u * Blu
 ! internal jump up
 actVal = (low_pop * Blu - up_pop * Bul) * betalu * exci_energy_l * Jlu
 IF(actVal < 0.D0) STOP 'i_radtrans: (l_pop * Blu - u_pop * Bul) < 0'
 ! write(*,*) 'i_radtrans: nB - nB = ', (low_pop * Blu - up_pop * Bul)
 ! write(*,*) 'i_radtrans: low_pop = ', low_pop, ' up_pop = ', up_pop, ' up_pop / low_pop = ', up_pop / low_pop
 actirates%Lma_int_uprad(I) = actVal
 ! write(*,*) 'i_radtrans: Lma_int_uprad = ', actVal
 Zintup = Zintup + actirates%Lma_int_uprad(I)
 ! write(*,*) 'i_radtrans: act_line = ', act_line, ' stat_weight = ', stat_weight, &
 !  ' exci_energy_l = ', exci_energy_l, ' up_pop = ', up_pop, ' linelist(act_line)%A_ul = ', &
 !  linelist(act_line)%A_ul
 ! write(*,*) 'i_radtrans: Zintup = ', Zintup
END DO

END SUBROUTINE i_radtrans
