!_________________________ i-packet rates ____________________________________________________
! calculation of radiative transitions rates of an i-packets
!
! * Zintdown -- total rate of internal downward jumps
! * Zrad -- total rate of radiative deactivations of a macro atom
! * Zintup -- total rate of internal upward jumps
SUBROUTINE i_radtrans(current_mgi, indexe, indexi, indexl, Zintdown, Zintup, Zrad, actirates, pack_index)
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
DOUBLE PRECISION                        :: Bul, Blu, Jlu, Aul
DOUBLE PRECISION                        :: flux_function
INTEGER                                 :: act_line
INTEGER                                 :: lower_level, upper_level
DOUBLE PRECISION                        :: actVal
INTEGER                                 :: I
! output variables
DOUBLE PRECISION                        :: Zintdown, Zintup, Zrad
TYPE(irates)      :: actirates
INTEGER                                 :: dummypackage, pack_index
DOUBLE PRECISION                        :: constanta
DOUBLE PRECISION                        :: ROverV, roverw
DOUBLE PRECISION                        :: fr_line
DOUBLE PRECISION                        :: corrFactor
! stimulated emission as negative absorption
LOGICAL                                 :: stmasnab=.false.

constanta = (pi * e_charge**2)/( me_g * light_speed)

! initialization of total rates
Zintdown = 0.D0
Zrad= 0.D0
Zintup = 0.D0

CALL populations(indexe, indexi, indexl, current_mgi, act_pop)

! number of transitions up and down
nlns = SIZE(elements(indexe)%ions(indexi)%levels(indexl)%linetransitions)
nluns = SIZE(elements(indexe)%ions(indexi)%levels(indexl)%lineuptransitions)
ALLOCATE(linetransitions(nlns), lineuptransitions(nluns))
linetransitions = elements(indexe)%ions(indexi)%levels(indexl)%linetransitions
lineuptransitions = elements(indexe)%ions(indexi)%levels(indexl)%lineuptransitions
! we have to know which element and ion we are calculating data for
! we will use the knowledge of lines

! transitions down
! the actual population is now the upper population
up_pop = act_pop
! write(36,*) '********************************************************'
! write(36,*) 'i_radtrans: nlns = ', nlns
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
dummypackage = SIZE(package)
package(dummypackage) = package(pack_index)
CALL emit_rpackage(dummypackage)
! write(*,*) 'i_radtrans: dir = ', package(dummypackage)%dir
! internal downward jump and radiative deexcitation
DO I = 1, nlns
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! 1.) internal downward jump
 act_line = linetransitions(I)
 lower_level = linelist(act_line)%lower
 upper_level = linelist(act_line)%upper
 ! the basic variables
 stat_weight_u = elements(indexe)%ions(indexi)%levels(upper_level)%stat_waight
 stat_weight_l = elements(indexe)%ions(indexi)%levels(lower_level)%stat_waight
 exci_energy_u = elements(indexe)%ions(indexi)%levels(upper_level)%exci_energy
 exci_energy_l = elements(indexe)%ions(indexi)%levels(lower_level)%exci_energy
 fr_line = linelist(act_line)%freq
 ! population of lower level
 CALL populations(indexe, indexi, lower_level, current_mgi, low_pop)
 corrFactor = 1.D0 - (DBLE(stat_weight_l) * up_pop) / (DBLE(stat_weight_u) * low_pop)
 IF(corrFactor < 0.D0) write(*,*) 'WARNING: correction factor 1 - (gl nu) / (gu nl) < 0'
 !______________________________________________________________________
 ! internal downward jump & radiative deexcitation
 ! calculation of R/V
 ! homologous approximation
 ROverV = roverw()
 ! optical depth
 !  * linelist(act_line)%f_lu * corrFactor
 taulu = light_speed / fr_line * constanta * &
   linelist(act_line)%f_lu * up_pop * ROverV * corrFactor
 betalu = 1.D0 / taulu * (1.D0 - exp(- taulu))
 Blu = 4 * pi**2 * e_charge**2 / (me_g * light_speed * h * fr_line) * linelist(act_line)%f_lu
 Bul = stat_weight_l / stat_weight_u * Blu
 ! Aul = 2 * h * fr_line**3 / light_speed**2 * Bul
 Aul = 8.D0 * fr_line**2 * pi**2 * e_charge**2/ (me_g * light_speed**3) *&
  stat_weight_l / stat_weight_u * linelist(act_line)%f_lu
 Jlu = flux_function(0, linelist(act_line)%freq, model_grid(current_mgi)%T, model_grid(current_mgi)%rwind)
 IF(stmasnab) THEN
  actVal = (low_pop * Blu - up_pop * Bul) * betalu * Jlu
 ELSE
  ! Blu = light_speed ** 2.0 / (2.0 * h * fr_line**3.0) * DBLE(stat_weight_u) / DBLE(stat_weight_l) &
  !  * linelist(act_line)%A_ul
  ! calculation of Blu and Bul
  ! actVal = up_pop * betalu * (Aul + Bul * Jlu)! * corrFactor 
  actVal = Aul * betalu * up_pop
 END IF
 ! write(*,*) 'i_radtrans: actVal = ', actVal, ' e_l = ', exci_energy_l, ' e_u - e_l = ', exci_energy_u - exci_energy_l
 actirates%Lma_int_dorad(I) = actVal * exci_energy_l
 ! actirates%Lma_int_dorad(I) = 2.0 * fr_line **2 / light_speed**2 * stat_weight_l / stat_weight_u * &
 !  linelist(act_line)%f_lu * betalu * exci_energy_l
 IF(actirates%Lma_int_dorad(I) < 0.D0) STOP 'i_radtrans: Lma_int_dorad < 0'
 ! write(*,*) 'i_radtrans: exci_energy_u = ', exci_energy_u
 actirates%Lma_rad(I) = actVal * (exci_energy_u - exci_energy_l)
 ! actirates%Lma_rad(I) = 2.0 * fr_line **2 / light_speed**2 * stat_weight_l / stat_weight_u * &
 !  linelist(act_line)%f_lu * betalu * exci_energy_l
 Zintdown = Zintdown + actirates%Lma_int_dorad(I)
 ! write(*,*) 'i_radtrans: Zintdown = ', Zintdown
 Zrad = Zrad + actirates%Lma_rad(I)
 ! write(*,*) 'i_radtrans: Zrad = ', Zrad
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
 ! calculation of up_pop
 CALL populations(indexe, indexi, linelist(act_line)%upper, current_mgi, up_pop)
 corrFactor = 1.D0 - (DBLE(stat_weight_l) * up_pop) / (DBLE(stat_weight_u) * low_pop)
 IF(corrFactor < 0.D0) write(*,*) 'WARNING: correction factor 1 - (gl nu) / (gu nl) < 0'
 Jlu = flux_function(1, linelist(act_line)%freq, model_grid(current_mgi)%T, model_grid(current_mgi)%rwind)
 ! calculation of Blu and Bul
 ! Blu = light_speed ** 2.0 / (2.0 * h * linelist(act_line)%freq**3.0) * DBLE(stat_weight_u) / DBLE(stat_weight_l) &
 !  * linelist(act_line)%A_ul
 fr_line = linelist(act_line)%freq
 Blu = 4 * pi**2 * e_charge**2 / (me_g * light_speed * h * fr_line) * linelist(act_line)%f_lu
 Bul = DBLE(stat_weight_l) / DBLE(stat_weight_u) * Blu
 Aul = 2 * h * fr_line**3 / light_speed**2 * Bul
! internal jump up
! write(*,*) 'i_radtrans: pack_index = ', pack_index, ' freq = ', linelist(act_line)%freq
! IF(vec_length(vel_vec) /= 0.D0) THEN
 ROverV = roverw()
 ! optical depth
 ! taulu = low_pop * pi * e_v ** 2.0  * ROverV / (me_g * light_speed * linelist(act_line)%freq) &
 !  * linelist(act_line)%f_lu * corrFactor
 taulu = light_speed / linelist(act_line)%freq * constanta * &
  linelist(act_line)%f_lu * low_pop * corrFactor * ROverV
 betalu = 1.D0 / taulu * (1.D0 - exp(- taulu))
 IF(stmasnab) THEN
  actVal = up_pop * betalu * linelist(act_line)%A_ul
 ELSE
  ! actVal = up_pop * betalu * (Aul + Bul * Jlu)
  actVal = (Blu * low_pop - Bul * up_pop) * betalu * Jlu
 END IF
 IF(actVal < 0.D0) STOP 'i_radtrans: (l_pop * Blu - u_pop * Bul) < 0'
 actirates%Lma_int_uprad(I) = actVal * exci_energy_l
 ! actirates%Lma_int_uprad(I) = linelist(act_line)%f_lu / (h * fr_line) * corrFactor * &
 !  betalu * Jlu * exci_energy_l
!  write(*,*) 'i_radtrans: up_pop = ', up_pop, 'low_pop = ', low_pop
!  write(*,*) 'i_radtrans: nB - nB = ', (low_pop * Blu - up_pop * Bul)
!  write(*,*) 'i_radtrans: taulu = ', taulu, ' corrFactor = ', corrFactor
!  write(*,*) 'i_radtrans: betalu = ', betalu, ' Jlu = ', Jlu
 ! write(36, *) 'i_radtrans: ', linelist(act_line)%lower, '->',&
 !  linelist(act_line)%upper, ' Lup = ', actirates%Lma_int_uprad(I),&
 !  ' wale = ', 1.D8 * light_speed / linelist(act_line)%freq
 Zintup = Zintup + actirates%Lma_int_uprad(I)
!  write(*,*) 'i_radtrans: Lma_int_uprad = ', actirates%Lma_int_uprad(I)
!  write(*,*) 'i_radtrans: Zintup = ', Zintup
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
