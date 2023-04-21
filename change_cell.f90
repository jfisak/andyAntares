SUBROUTINE change_cell(pack_index, next_cell)

 ! Change the cell number

 USE types
 USE counters

 IMPLICIT NONE    

 INTEGER                           :: pack_index, next_cell
 ! DOUBLE PRECISION                  :: event_dist
 DOUBLE PRECISION, DIMENSION(3)         :: pos, width, corner

INTEGER                                 :: get_package_model_index
INTEGER                                 :: pomocna_bunka, I, old_cell
INTEGER                                 :: dummypackage

DOUBLE PRECISION                        :: round_number

dummypackage = SIZE(package)

! If package escaped the calculation volume (next_cell=-99) then it
! become no-active and update into type_escaped, else it update the
! cell number
!  if(next_cell > 0) write(22,*) next_cell, dyn_cell(next_cell)%corner, dyn_cell(next_cell)%width
! write(*,*) 'change_cell: pack_index = ', pack_index
 pos = package(pack_index)%pos
 old_cell = package(pack_index)%cell_numb
!  write(*,*) 'change_cell: old_cell= ', old_cell 
!  write(*,*) 'change_cell: neighbors = ', dyn_cell(old_cell)%neighbor
!  write(*,*) 'change_cell: pos = ', pos
 ! CALL find_dyn_cell1(pos, next_cell)
 IF (next_cell .LT. 0) THEN 
  package(pack_index)%typ = type_escaped
  package(pack_index)%active = 0
  count_des_esca = count_des_esca + 1
 ELSE
  corner = dyn_cell(next_cell)%corner
  width = dyn_cell(next_cell)%width
  ! write(*,*)  'change cell: next_cell = ', next_cell
  package(pack_index)%cell_numb = next_cell
  IF(next_cell > SIZE(dyn_cell)) THEN
   CALL abort()
  END IF
  ! check if the next cell si chosen correctly
  ! DO I = 1,3
  !  ! pos(I) = round_number(pos(I))
  !  IF((pos(I) < corner(I) .OR. pos(I) > corner(I) + width(I)) .and. pack_index /= dummypackage ) THEN
  !   ! CALL find_dyn_cell1(pos, pomocna_bunka)
  !   write(4,*) pos, dyn_cell(old_cell)%corner, dyn_cell(old_cell)%width, get_package_model_index(pack_index)
  !   write(4,*) pos, corner, width, get_package_model_index(pack_index)
  !   write(*,*) 'I = ', I
  !   write(*,*) (pos(I)-corner(I)), (pos(I)-corner(I)-width(I))
  !   write(*,*) (pos(I)-corner(I))/width(I), (pos(I)-corner(I)-width(I))/width(I)
  !   write(*,*) '***************************************************'
  !   write(*,*) 'change_cell: pack_index = ', pack_index
  !   write(*,*) 'change_cell: next_cell = ', next_cell
  !   STOP 'change_cell: packet is not located inside the propagation cell'
  !  END IF
  ! END DO
 END IF
 
END SUBROUTINE change_cell
