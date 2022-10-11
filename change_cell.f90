SUBROUTINE change_cell(pack_index, next_cell)

 ! Change the cell number

 USE types
 USE counters

 IMPLICIT NONE    

 INTEGER                           :: pack_index, next_cell
 ! DOUBLE PRECISION                  :: event_dist
 DOUBLE PRECISION, DIMENSION(3)         :: pos

INTEGER                                 :: get_package_model_index
INTEGER                                 :: pomocna_bunka, I, old_cell, next_mgi
INTEGER                                 :: dummypackage

DOUBLE PRECISION                        :: round_number

LOGICAL                                 :: is_difap

dummypackage = SIZE(package)

! If package escaped the calculation volume (next_cell=-99) then it
! become no-active and update into type_escaped, else it update the
! cell number
!  if(next_cell > 0) write(22,*) next_cell, dyn_cell(next_cell)%corner, dyn_cell(next_cell)%width
! write(*,*) 'change_cell: pack_index = ', pack_index
pos = package(pack_index)%pos
old_cell = package(pack_index)%cell_numb
IF (next_cell .LT. 0) THEN 
 package(pack_index)%typ = type_escaped
 package(pack_index)%active = 0
 count_des_esca = count_des_esca + 1
ELSE
 package(pack_index)%cell_numb = next_cell
 IF(next_cell > SIZE(dyn_cell)) THEN
  CALL abort()
 END IF
 if(enable_diffusion == 1) then
  next_mgi = dyn_cell(next_cell)%model_index
  is_difap = model_grid(next_mgi)%is_difapp
  IF(is_difap) THEN
   package(pack_index)%typ = type_dpkt
  ELSE
   package(pack_index)%typ = type_rpkt
  END IF
 end if
END IF
 
END SUBROUTINE change_cell
