SUBROUTINE change_cell(pack_index, next_cell)

  ! Change the cell number

  USE types
  USE counters

  IMPLICIT NONE    

  INTEGER                           :: pack_index, next_cell
  ! DOUBLE PRECISION                  :: event_dist


  ! If package escaped the calculation volume (next_cell=-99) then it
  ! become no-active and update into type_escaped, else it update the
  ! cell number
!  if(next_cell > 0) write(22,*) next_cell, dyn_cell(next_cell)%corner, dyn_cell(next_cell)%width
  IF (next_cell .LT. 0) THEN 
     package(pack_index)%typ = type_escaped
     package(pack_index)%active = 0
     !$OMP ATOMIC
     count_des_esca = count_des_esca + 1
!     print*, 'change cell: package = ', pack_index, ' package escaped...'
  ELSE
!     print*, 'change cell: next_cell = ', next_cell
     IF(next_cell > SIZE(dyn_cell)) THEN
      write(99,*) 'change_cell: pack_index = ', pack_index, &
       ' next_cell = ', next_cell, ' is larger then total number of cells'
      CALL abort()
     END IF
     package(pack_index)%cell_numb = next_cell
  END IF
 
END SUBROUTINE change_cell
