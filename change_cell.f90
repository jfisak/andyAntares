SUBROUTINE change_cell(pack_index, next_cell)

  ! Change the cell number

  USE types

  IMPLICIT NONE    

  INTEGER                           :: pack_index, next_cell
  ! DOUBLE PRECISION                  :: event_dist


  ! If package escaped the calculation volume (next_cell=-99) then it
  ! become no-active and update into type_escaped, else it update the
  ! cell number
  IF (next_cell .LT. 0) THEN 
     package(pack_index)%typ = type_escaped
     package(pack_index)%active = 0
  ELSE
     package(pack_index)%cell_numb = next_cell
  END IF
 
END SUBROUTINE change_cell
