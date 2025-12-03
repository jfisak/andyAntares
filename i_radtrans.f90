! calculation of radiative transitions rates of an i-packets
!
! INPUT: current_mgi(INT): modGrid index
!        indexe(INT): element index
!        indexi(INT): ion index
!        indexl(INT): level index
!        pack_index(INT): index of the packet
! OUTPUT: Zintdown(INT): total rate of internal downward jumps
!         Zintup(INT): total rate of internal upward jumps
!         Zrad(INT): total rate of radiative deactivations of a macro atom
!         actirates(irates): rates for the corresponding transitions
!
SUBROUTINE i_radtrans(current_mgi, indexe, indexi, indexl, Zintdown, Zintup, Zrad, actirates,&
                      pack_index, new_direction)
USE types
USE constants
USE rates_i
USE dummypacket
IMPLICIT NONE

! input variables
INTEGER                                 :: nlns, nluns
INTEGER                                 :: indexe, indexi, indexl
INTEGER                                 :: current_mgi
DOUBLE PRECISION, DIMENSION(const_dimofspace) :: new_direction
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
INTEGER                                 :: ind_I
! output variables
DOUBLE PRECISION                        :: Zintdown, Zintup, Zrad
TYPE(irates)      :: actirates
INTEGER                                 :: pack_index
DOUBLE PRECISION                        :: constanta
DOUBLE PRECISION                        :: ROverV, roverw
DOUBLE PRECISION                        :: fr_line
DOUBLE PRECISION                        :: corrFactor
! stimulated emission as negative absorption
LOGICAL                                 :: stmasnab=.false.

DOUBLE PRECISION                        :: cur_r
INTEGER                                 :: cur_dummypack, dummypack_index

constanta = (const_pi * const_e**2)/( const_me_g * const_c)

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
cur_dummypack = find_free_index()
dummypack_index = cur_dummypack + SIZE(package)
CALL copy_package(pack_index, cur_dummypack)
CALL emit_rpackage(dummypack_index)
new_direction = dummypackage(cur_dummypack)%dir
CALL deactivate_dummy_packet(cur_dummypack)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!______________________________________________________________________
! internal downward jump and radiative deexcitation
DO ind_I = 1, nlns
 act_line = linetransitions(ind_I)
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
 ROverV = roverw(pack_index, 0.D0, fr_line)
 write(*,*) 'i_radtrans: roverw = ', ROverV


 taulu = const_c / fr_line * constanta * &
   linelist(act_line)%f_lu * up_pop * ROverV * corrFactor

write(*,*) 'i_radtrans: taulu = ', taulu

 IF(taulu == 0.D0) THEN
  betalu = 1.D0
 ELSE IF(taulu > 0.D0) THEN
  betalu = 1.D0 / taulu * (1.D0 - exp(- taulu))
 ELSE
  write(*,*) 'i_radtrans: betalu = ', betalu
  STOP 'i_radtrans, betalu < 0'
 END IF


 Blu = 4 * const_pi**2 * const_e**2 / (const_me_g * const_c * const_h * fr_line) * linelist(act_line)%f_lu
 Bul = stat_weight_l / stat_weight_u * Blu
 Aul = 8.D0 * fr_line**2 * const_pi**2 * const_e**2/ (const_me_g * const_c**3) *&
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



 actirates%Lma_int_dorad(ind_I) = actVal * exci_energy_l



 IF(actirates%Lma_int_dorad(ind_I) < 0.D0) STOP 'i_radtrans: Lma_int_dorad < 0'
 
 actirates%Lma_rad(ind_I) = actVal * (exci_energy_u - exci_energy_l)
 Zintdown = Zintdown + actirates%Lma_int_dorad(ind_I)
 Zrad = Zrad + actirates%Lma_rad(ind_I)

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
DO ind_I = 1, nluns

 act_line = lineuptransitions(ind_I)
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


 Blu = 4 * const_pi**2 * const_e**2 / (const_me_g * const_c * const_h * fr_line) * linelist(act_line)%f_lu
 Bul = DBLE(stat_weight_l) / DBLE(stat_weight_u) * Blu
 Aul = 2 * const_h * fr_line**3 / const_c**2 * Bul


 ROverV = roverw(pack_index, R_star, fr_line)

 taulu = const_c / fr_line * constanta * &
  linelist(act_line)%f_lu * low_pop * corrFactor * ROverV


 IF(taulu == 0.D0) THEN
  betalu = 1.D0
 ELSE IF(taulu > 0.D0) THEN
  betalu = 1.D0 / taulu * (1.D0 - exp(- taulu))
 ELSE
  write(*,*) 'i_radtrans: betalu = ', betalu
  STOP 'i_radtrans, betalu < 0'
 END IF

 IF(stmasnab) THEN
  actVal = up_pop * betalu * linelist(act_line)%A_ul
 ELSE
  actVal = (Blu * low_pop - Bul * up_pop) * betalu * Jlu
 END IF

 ! write(*,*) 'i_radtrans: Blu = ', Blu, ' Bul = ', Bul, ' up_pop = ', up_pop, ' Jlu = ', Jlu


 IF(actVal < 0.D0) STOP 'i_radtrans: (l_pop * Blu - u_pop * Bul) < 0'


 actirates%Lma_int_uprad(ind_I) = actVal * exci_energy_l
 Zintup = Zintup + actirates%Lma_int_uprad(ind_I)

END DO

END SUBROUTINE i_radtrans
