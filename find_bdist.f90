SUBROUTINE find_bdist(pack_index, basic_cell_numb, bdist)

USE types

IMPLICIT NONE

! input variables
INTEGER                         :: pack_index, basic_cell_numb
! output variables
DOUBLE PRECISION                :: bdist
! variables for photon
DOUBLE PRECISION, DIMENSION(3)  :: dir, phot_pos
! variables for basic cells
INTEGER                         :: actCell
INTEGER                         :: next_bas_cell, next_bcell
DOUBLE PRECISION, DIMENSION(3)  :: bcorner, bwidth
DOUBLE PRECISION                :: bt1, bt2, bt3, bt4, bt5, bt6

  bcorner = dyn_cell(basic_cell_numb)%corner
  bwidth = dyn_cell(basic_cell_numb)%width
  phot_pos = package(pack_index)%pos
  dir = package(pack_index)%dir
 ! we will calculate parameters t1,...,t6
 IF(dir(1) /= 0) THEN
  bt1 = (bcorner(1) - phot_pos(1))/(dir(1))
  bt4 = (bcorner(1) + bwidth(1) - &
        phot_pos(1))/(dir(1))
 ELSE
  bt1 = 0
  bt4 = 0
 END IF
 IF(dir(2) /= 0) THEN
  bt2 = (bcorner(2) - phot_pos(2))/(dir(2))
  bt5 = (bcorner(2) + bwidth(2) - &
        phot_pos(2))/(dir(2))
 ELSE
  bt2 = 0
  bt5 = 0
 END IF
 IF(dir(3) /= 0) THEN
  bt3 = (bcorner(3) - phot_pos(3))/(dir(3))
  bt6 = (bcorner(3) + bwidth(3) - &
        phot_pos(3))/(dir(3))
 ELSE
  bt3 = 0
  bt6 = 0
 END IF
 print*, 'find_bdist: bt1 = ', bt1, ' bt2 = ', bt2, ' bt3 = ', &
        bt3, ' bt4 = ', bt4, ' bt5 = ', bt5, ' bt6 = ', bt6

  bdist = 1.D99
  ! we are looking for a bound in front of the photon,
  ! so we have to choose solution with t > 0
  IF( (bt1 > 0.E0) .AND. (bt1 < bdist) ) THEN
   bdist = bt1
!   IF(dyn_cell(basic_cell_numb)%indexc(1) == 1) THEN
!    next_bcell = -99
!   ELSE
!    next_bcell = basic_cell_numb + ny_cell * nz_cell
!   END IF
  END IF
  IF( (bt2 > 0.E0) .AND. (bt2 < bdist) ) THEN
   bdist = bt2
!   IF(dyn_cell(basic_cell_numb)%indexc(1) == 1) THEN
!    next_bcell = -99
!   ELSE
!    next_bcell = basic_cell_numb + nz_cell
!   END IF
  END IF
  IF( (bt3 > 0.E0) .AND. (bt3 < bdist)  ) THEN
   bdist = bt3
!   IF(dyn_cell(basic_cell_numb)%indexc(1) == 1) THEN
!    next_bcell = -99
!   ELSE
!    next_bcell = basic_cell_numb + 1
!   END IF
  END IF
  IF( (bt4 > 0.E0) .AND. (bt4 < bdist)  ) THEN
   bdist = bt4
!   IF(dyn_cell(basic_cell_numb)%indexc(1) == nx_cell) THEN
!    next_bcell = -99
!   ELSE
!    next_bcell = basic_cell_numb - ny_cell * nz_cell
!   END IF
  END IF
  IF( (bt5 > 0.E0) .AND. (bt5 < bdist)  ) THEN
   bdist = bt5
!   IF(dyn_cell(basic_cell_numb)%indexc(1) == ny_cell) THEN
!    next_bcell = -99
!   ELSE
!    next_bcell = basic_cell_numb - nz_cell
!   END IF
  END IF
  IF( (bt6 > 0.E0) .AND. (bt6 < bdist)  ) THEN
   bdist = bt6
!   IF(dyn_cell(basic_cell_numb)%indexc(1) == nz_cell) THEN
!    next_bcell = -99
!   ELSE
!    next_bcell = basic_cell_numb - 1
!   END IF
  END IF
  print*, 'find_bdist: next_bcell = ', next_bcell, ' bdist = ', bdist

END SUBROUTINE find_bdist
