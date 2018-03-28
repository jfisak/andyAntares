!_________________________ i-packet rates ____________________________________________________
! calculation of radiative transitions rates of an i-packets
!
! * Zintdown -- total rate of internal downward jumps
! * Zrad -- total rate of radiative deactivations of a macro atom
! * Zintup -- total rate of internal upward jumps
SUBROUTINE i_radtrans(current_mgi, indexe, indexi, indexl, act_pop, Zintdown, Zintup, Zrad, actirates)
USE types
USE rates_i
IMPLICIT NONE

! input variables
INTEGER                                 :: nlns, nluns
INTEGER                                 :: indexe, indexi, indexl
INTEGER                                 :: current_mgi
INTEGER, ALLOCATABLE                    :: linetransitions(:)
INTEGER, ALLOCATABLE                    :: lineuptransitions(:)
DOUBLE PRECISION                        :: stat_weight_l, stat_weight_u
DOUBLE PRECISION                        :: exci_energy_u, exci_energy_l
DOUBLE PRECISION                        :: act_pop, up_pop, low_pop
! beta lu calculation
DOUBLE PRECISION                        :: taulu, betalu
DOUBLE PRECISION                        :: Bul, Blu, Jlu
DOUBLE PRECISION                        :: flux_function
INTEGER                                 :: act_line
DOUBLE PRECISION                        :: actVal
INTEGER                                 :: I
! output variables
DOUBLE PRECISION                        :: Zintdown, Zintup, Zrad
INTEGER                         :: OMP_GET_THREAD_NUM, my_rank
TYPE(irates)      :: actirates
INTEGER                                 :: dummypackage, n_pack_d
DOUBLE PRECISION                        :: constant

my_rank = OMP_GET_THREAD_NUM()
n_pack_d = SIZE(package)
dummypackage = n_pack_d - n_dummy_packs + my_rank + 1

constant = (pi * e_charge**2)/( me_g * light_speed)

Zintdown = 0.D0
Zrad= 0.D0
Zintup = 0.D0
nlns = SIZE(elements(indexe)%ions(indexi)%levels(indexl)%linetransitions)
nluns = SIZE(elements(indexe)%ions(indexi)%levels(indexl)%lineuptransitions)
ALLOCATE(linetransitions(nlns), lineuptransitions(nluns))
linetransitions = elements(indexe)%ions(indexi)%levels(indexl)%linetransitions
lineuptransitions = elements(indexe)%ions(indexi)%levels(indexl)%lineuptransitions
! we have to know which element and ion we are calculating data for
! we will use the knowledge of lines and assume that at it is
! possible at least one transition upwards or downwards and from
! the first element we get the element and the ion informations

up_pop = act_pop
! write(36,*) '***************************************************************************************'
! write(36,*) 'i_radtrans: nlns = ', nlns
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! internal downward jump and radiative deexcitation
DO I = 1, nlns
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! 1.) internal downward jump
 act_line = linetransitions(I)
 ! the basic variables
 stat_weight_u = elements(indexe)%ions(indexi)%levels(linelist(act_line)%upper)%stat_waight
 stat_weight_l = elements(indexe)%ions(indexi)%levels(linelist(act_line)%lower)%stat_waight
 exci_energy_u = elements(indexe)%ions(indexi)%levels(linelist(act_line)%upper)%exci_energy
 exci_energy_l = elements(indexe)%ions(indexi)%levels(linelist(act_line)%lower)%exci_energy
 ! internal downward jump
 ! calculation of a rate coefficient
 CALL populations(indexe, indexi, linelist(act_line)%lower, current_mgi, low_pop)
 ! Einstein Blu coefficient
 Blu = light_speed**2.0 / (2.0 * h * linelist(act_line)%freq**3.0) * stat_weight_u / stat_weight_l &
  * linelist(act_line)%A_ul
 ! optical depth
 taulu = low_pop * Blu * h * light_speed * (R_inf / V_inf) / (4.0 * pi) * &
  (1.D0 - (stat_weight_l * up_pop) / (stat_weight_u * low_pop))
 betalu = 1.D0 / taulu * (1.D0 - exp(- taulu))
 actVal = up_pop * betalu * linelist(act_line)%A_ul
 actirates%Lma_int_dorad(I) = actVal * exci_energy_l
 IF(actirates%Lma_int_dorad(I) < 0.D0) STOP 'i_radtrans: Lma_int_dorad < 0'
 Zintdown = Zintdown + actirates%Lma_int_dorad(I)
 actirates%Lma_rad(I) = actVal * (exci_energy_u - exci_energy_l)
 Zrad = Zrad + actirates%Lma_rad(I)
 ! write(*,*) 'i_radtrans: Lrad = ', actirates%Lma_rad(I), ' Lintdown = ', actirates%Lma_int_dorad(I)
 ! write(*,*) 'i_radtrans: e_l = ', exci_energy_l, ' e_u - u_l = ', exci_energy_u - exci_energy_l
 IF(exci_energy_u - exci_energy_l < 0) STOP 'i_radtrans: exci_energy_u - exci_energy_l < 0'
END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! internal upward jump
! write(36,*) 'i_radtrans: nluns = ', nluns
low_pop = act_pop
DO I = 1, nluns
 act_line = lineuptransitions(I)
 exci_energy_l = elements(indexe)%ions(indexi)%levels(linelist(act_line)%lower)%exci_energy
 exci_energy_u = elements(indexe)%ions(indexi)%levels(linelist(act_line)%upper)%exci_energy
 stat_weight_l = elements(indexe)%ions(indexi)%levels(linelist(act_line)%lower)%stat_waight
 stat_weight_u = elements(indexe)%ions(indexi)%levels(linelist(act_line)%upper)%stat_waight
! CALL populations(indexe, indexi, linelist(act_line)%lower, current_mgi, low_pop)
 CALL populations(indexe, indexi, linelist(act_line)%upper, current_mgi, up_pop)
 Jlu = flux_function(0, linelist(act_line)%freq, model_grid(current_mgi)%T)
 ! calculation of Blu and Bul
 Blu = light_speed ** 2.0 / (2.0 * h * linelist(act_line)%freq**3.0) * stat_weight_u / stat_weight_l &
  * linelist(act_line)%A_ul
 Bul = stat_weight_l / stat_weight_u * Blu
! internal jump up
! write(*,*) 'i_radtrans: pack_index = ', pack_index, ' freq = ', linelist(act_line)%freq
! CALL move_package(dummypackage, l_dist)
! CALL velo(dummypackage, vel_vec)
! IF(vec_length(vel_vec) /= 0.D0) THEN
 taulu = low_pop * Blu * h * light_speed * (R_inf / V_inf) / (4.0 * pi) * &
  (1.D0 - (stat_weight_l * up_pop) / (stat_weight_u * low_pop))
 betalu = 1.D0 / taulu * (1.D0 - exp(- taulu))
 actVal = (low_pop * Blu - up_pop * Bul) * betalu  * Jlu
 IF(actVal < 0.D0) STOP 'i_radtrans: (l_pop * Blu - u_pop * Bul) < 0'
 actirates%Lma_int_uprad(I) = actVal * exci_energy_l
 ! write(36,*) 'i_radtrans: nB - nB = ', (low_pop * Blu - up_pop * Bul)
 ! write(36,*) 'i_radtrans: betalu = ', betalu, ' Jlu = ', Jlu
 ! write(36, *) 'i_radtrans: ', linelist(act_line)%lower, '->',&
 !  linelist(act_line)%upper, ' Lup = ', actirates%Lma_int_uprad(I),&
 !  ' wale = ', 1.D8 * light_speed / linelist(act_line)%freq
 Zintup = Zintup + actirates%Lma_int_uprad(I)
 ! write(*,*) 'i_radtrans: wale = ', 1.D8 * light_speed / linelist(act_line)%freq
 ! write(*,*) 'i_radtrans: low_pop = ', low_pop, ' up_pop = ', up_pop, ' up_pop / low_pop = ', up_pop / low_pop
 ! write(*,*) 'i_radtrans: act_line = ', act_line, ' exci_energy_l = ', exci_energy_l, ' up_pop = ',&
 !  up_pop, ' linelist(act_line)%A_ul = ', linelist(act_line)%A_ul
! ELSE
!  actirates%Lma_int_uprad(I) = 0.D0
! END IF
 ! write(*,*) 'i_radtrans: Zintup = ', Zintup
END DO

END SUBROUTINE i_radtrans
