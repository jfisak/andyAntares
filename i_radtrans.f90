!_________________________ i-packet rates ____________________________________________________
! calculation of radiative transitions rates of an i-packets
!
! * Zintdown -- total rate of internal downward jumps
! * Zrad -- total rate of radiative deactivations of a macro atom
! * Zintup -- total rate of internal upward jumps
SUBROUTINE i_radtrans(current_mgi, indexe, indexi, indexl, Zintdown, Zintup, Zrad, actirates, pack_index)
USE types
USE constants
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

DOUBLE PRECISION                        :: cur_r

constanta = (const_pi * e_charge**2)/( me_g * light_speed)

cur_r = norm2(package(pack_index)%pos)
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
dummypackage = SIZE(package)
package(dummypackage) = package(pack_index)
CALL emit_rpackage(dummypackage)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!______________________________________________________________________
! internal downward jump and radiative deexcitation
DO I = 1, nlns
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
 ROverV = roverw(pack_index, R_star, fr_line)


 taulu = light_speed / fr_line * constanta * &
   linelist(act_line)%f_lu * up_pop * ROverV * corrFactor


 betalu = 1.D0 / taulu * (1.D0 - exp(- taulu))


 Blu = 4 * const_pi**2 * e_charge**2 / (me_g * light_speed * const_h * fr_line) * linelist(act_line)%f_lu
 Bul = stat_weight_l / stat_weight_u * Blu
 Aul = 8.D0 * fr_line**2 * const_pi**2 * e_charge**2/ (me_g * light_speed**3) *&
  stat_weight_l / stat_weight_u * linelist(act_line)%f_lu


 Jlu = flux_function(1, linelist(act_line)%freq, model_grid(current_mgi)%T, cur_r)
 
 !______________________________________________________________________
 ! stimulated emission as negative absorption?
 IF(stmasnab) THEN
  actVal = (low_pop * Blu - up_pop * Bul) * betalu * Jlu
 ELSE
  actVal = Aul * betalu * up_pop
 END IF
 !______________________________________________________________________



 actirates%Lma_int_dorad(I) = actVal * exci_energy_l



 IF(actirates%Lma_int_dorad(I) < 0.D0) STOP 'i_radtrans: Lma_int_dorad < 0'
 
 actirates%Lma_rad(I) = actVal * (exci_energy_u - exci_energy_l)
 Zintdown = Zintdown + actirates%Lma_int_dorad(I)
 Zrad = Zrad + actirates%Lma_rad(I)

 IF(exci_energy_u - exci_energy_l < 0) STOP 'i_radtrans: exci_energy_u - exci_energy_l < 0'
END DO



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! internal upward jump
! write(36,*) 'i_radtrans: nluns = ', nluns
low_pop = act_pop
DO I = 1, nluns

 act_line = lineuptransitions(I)
 lower_level = linelist(act_line)%lower
 upper_level = linelist(act_line)%upper


 exci_energy_l = elements(indexe)%ions(indexi)%levels(lower_level)%exci_energy
 exci_energy_u = elements(indexe)%ions(indexi)%levels(upper_level)%exci_energy
 stat_weight_l = elements(indexe)%ions(indexi)%levels(lower_level)%stat_waight
 stat_weight_u = elements(indexe)%ions(indexi)%levels(upper_level)%stat_waight


 CALL populations(indexe, indexi, upper_level, current_mgi, up_pop)


 corrFactor = 1.D0 - (DBLE(stat_weight_l) * up_pop) / (DBLE(stat_weight_u) * low_pop)

 fr_line = linelist(act_line)%freq

 IF(corrFactor < 0.D0) write(*,*) 'WARNING: correction factor 1 - (gl nu) / (gu nl) < 0'

 Jlu = flux_function(1, fr_line, model_grid(current_mgi)%T, cur_r)


 Blu = 4 * const_pi**2 * e_charge**2 / (me_g * light_speed * const_h * fr_line) * linelist(act_line)%f_lu
 Bul = DBLE(stat_weight_l) / DBLE(stat_weight_u) * Blu
 Aul = 2 * const_h * fr_line**3 / light_speed**2 * Bul


 ROverV = roverw(pack_index, R_star, fr_line)

 taulu = light_speed / fr_line * constanta * &
  linelist(act_line)%f_lu * low_pop * corrFactor * ROverV


 betalu = 1.D0 / taulu * (1.D0 - exp(- taulu))

 IF(stmasnab) THEN
  actVal = up_pop * betalu * linelist(act_line)%A_ul
 ELSE
  actVal = (Blu * low_pop - Bul * up_pop) * betalu * Jlu
 END IF

 ! write(*,*) 'i_radtrans: Blu = ', Blu, ' Bul = ', Bul, ' up_pop = ', up_pop, ' Jlu = ', Jlu


 IF(actVal < 0.D0) STOP 'i_radtrans: (l_pop * Blu - u_pop * Bul) < 0'


 actirates%Lma_int_uprad(I) = actVal * exci_energy_l
 Zintup = Zintup + actirates%Lma_int_uprad(I)

END DO

END SUBROUTINE i_radtrans
