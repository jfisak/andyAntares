SUBROUTINE change_cell(pack_index, next_cell)

 ! Change the cell number

 USE types
 USE counters

 IMPLICIT NONE    

 INTEGER                           :: pack_index, next_cell
 ! DOUBLE PRECISION                  :: event_dist
 DOUBLE PRECISION, DIMENSION(3)         :: pos, width, corner


! If package escaped the calculation volume (next_cell=-99) then it
! become no-active and update into type_escaped, else it update the
! cell number
!  if(next_cell > 0) write(22,*) next_cell, dyn_cell(next_cell)%corner, dyn_cell(next_cell)%width
! write(*,*) 'change_cell: pack_index = ', pack_index
 IF (next_cell .LT. 0) THEN 
  package(pack_index)%typ = type_escaped
  package(pack_index)%active = 0
  count_des_esca = count_des_esca + 1
! print*, 'change cell: package = ', pack_index, ' package escaped...'
 ELSE
! print*, 'change cell: next_cell = ', next_cell
  package(pack_index)%cell_numb = next_cell
  IF(next_cell > SIZE(dyn_cell)) THEN
   ! write(*,*) 'change_cell: the cell ', next_cell, ' does not exist...'
   CALL abort()
  END IF
  ! check if the next cell si chosen correctly
  ! corner = dyn_cell(next_cell)%corner
  ! width = dyn_cell(next_cell)%width
  ! pos = package(pack_index)%pos
  ! IF(pos(
 END IF
 
END SUBROUTINE change_cell
