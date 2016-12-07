! this subroutine finds a next basic cell
SUBROUTINE find_basic_cell(cell_numb, pack_index, bdist, next_b_cell)

USE types

IMPLICIT NONE

! input variables
INTEGER                         :: cell_numb, pack_index
DOUBLE PRECISION                :: bdist
! output variables
INTEGER                         :: next_b_cell
DOUBLE PRECISION, DIMENSION(3)  :: crossp, width, corner
INTEGER, DIMENSION(6)           :: temp_cell
INTEGER                         :: I, found
!LOGICAL                         :: exists

found = 0
! the cross position calculation
crossp = package(pack_index)%pos + package(pack_index)%dir * bdist

! we can check now if the neigbour cells contain the crosspoint
! x_plus
temp_cell(1) = cell_numb + ny_cell * nz_cell
! x_minus
temp_cell(2) = cell_numb - ny_cell * nz_cell
! y_plus
temp_cell(3) = cell_numb + nz_cell
! y_minus
temp_cell(4) = cell_numb - nz_cell
! z_plus
temp_cell(5) = cell_numb + 1
! z_minus
temp_cell(6) = cell_numb - 1

!inquire(file='neighbours.dat', exist=exists)
!if(exists .EQV. .FALSE.) then
!OPEN(33,file='neighbours.dat')
! DO I = 1,6
!  print*, 'NEIGHBOUR: ', dyn_cell(temp_cell(I))%corner, dyn_cell(temp_cell(I))%width
!  write(33,*) dyn_cell(temp_cell(I))%corner, dyn_cell(temp_cell(I))%width
! END DO
!CLOSE(33)
!end if

DO I = 1, 6
 IF(temp_cell(I) > 0) THEN
  corner = dyn_cell(temp_cell(I))%corner
  width = dyn_cell(temp_cell(I))%width
  IF((crossp(1) == corner(1) .OR. crossp(1) == corner(1) + width(1) .OR. &
     crossp(2) == corner(2) .OR. crossp(2) == corner(2) + width(2) .OR. &
     crossp(3) == corner(3) .OR. crossp(3) == corner(3) + width(3)) .AND. &
     crossp(1) - corner(1) > 0.D0 .AND. abs(crossp(1) - corner(1)) <= width(1) .AND. &
     crossp(2) - corner(2) > 0.D0 .AND. abs(crossp(2) - corner(2)) <= width(2) .AND. &
     crossp(3) - corner(3) > 0.D0 .AND. abs(crossp(3) - corner(3)) <= width(3) ) THEN
     next_b_cell = temp_cell(I)
 !    print*, 'find_basic_cell: found = ', found, ' cell number: ', temp_cell(I)
     found = found + 1
  END IF
 END IF
END DO

!print*, 'find_basic_cell: found = ', found
IF(found == 0) next_b_cell = -99

!print*, 'find_basic_cell: next_b_cell = ', next_b_cell
END SUBROUTINE find_basic_cell
