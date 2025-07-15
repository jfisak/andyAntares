SUBROUTINE change_cell(pack_index, next_cell)

 ! Change the cell number

 USE types
 USE counters
USE constants

 IMPLICIT NONE    

 INTEGER                           :: pack_index, next_cell
 ! DOUBLE PRECISION                  :: event_dist
 DOUBLE PRECISION, DIMENSION(const_dimofspace)         :: pos

INTEGER                                 :: old_cell, next_mgi
INTEGER                                 :: dummypackage
DOUBLE PRECISION, PARAMETER             :: mininum = 1.D4

LOGICAL                                 :: is_difap
DOUBLE PRECISION, DIMENSION(const_dimofspace)           :: corner, upcorner
INTEGER                                                 :: ind_I, pom_pgi

dummypackage = SIZE(package)

! If package escaped the calculation volume (next_cell=-99) then it
! become no-active and update into type_escaped, else it update the
! cell number
! write(*,*) 'change_cell: pack_index = ', pack_index
pos = package(pack_index)%pos
old_cell = package(pack_index)%cell_numb
if(debug == 2) write(*,*) 'change_cell: next_cell = ', next_cell
if(debug == 2) write(*,*) 'change_cell: neighbors = ', dyn_cell(old_cell)%neighbor
IF (next_cell .LT. 0) THEN 
 package(pack_index)%typ = type_escaped
 package(pack_index)%active = 0
 count_des_esca = count_des_esca + 1
 IF(debug == 2) write(*,*) 'change_cell: packet escaped'
 IF(debug == 2) write(*,*) 'change_cell: ||r|| = ', norm2(package(pack_index)%pos)/R_inf
ELSE
 package(pack_index)%cell_numb = next_cell
 IF(debug == 2) THEN
  corner = dyn_cell(next_cell)%corner
  upcorner = dyn_cell(next_cell)%upcorner
  write(*,*) 'change_cell: change of cell into: ', next_cell
  write(*,*) 'change_cell: corner = ', dyn_cell(next_cell)%corner
  write(*,*) 'change_cell: pos = ', pos
  write(*,*) 'change_cell: upcorner = ', dyn_cell(next_cell)%upcorner
   DO ind_I = 1, const_dimofspace
    IF(((pos(ind_I) < corner(ind_I) - mininum) .OR. (pos(ind_I) > upcorner(ind_I) + mininum)) &
    .and. pack_index <= SIZE(package)) THEN
     write(*,*) 'change_cell: ind_I = ', ind_I
     CALL find_dyn_cell1(pos, pom_pgi)
     write(*,*) 'change_cell: pom_pgi = ', pom_pgi
     write(*,*) 'change_cell: corner = ', dyn_cell(pom_pgi)%corner
     write(*,*) 'change_cell: pos = ', pos
     write(*,*) 'change_cell:  = ', dyn_cell(pom_pgi)%upcorner
     STOP 'change_cell: packet is not located inside the propagation cell'
    END IF
   END DO
 END IF

 ! test on the propGrid cell index number
 IF(next_cell > SIZE(dyn_cell)) THEN
  write(*,*) 'change_cell: r/R_inf = ', norm2(package(pack_index)%pos)/R_inf
  write(*,*) 'change_cell: next_cell = ', next_cell
  write(*,*) 'change_cell: next_cell > number of propGrid cells, exiting now'
  CALL abort()
 END IF

 ! takes into account if the diffusion approximation takes place
 if(enable_diffusion == 1) then
  next_mgi = dyn_cell(next_cell)%model_index
  is_difap = model_grid(next_mgi)%is_difapp
  IF(is_difap) THEN
   package(pack_index)%typ = type_dpkt
   ! if(package(pack_index)%typ == type_rpkt) package(pack_index)%typ = type_dpkt
   ! if(package(pack_index)%typ == type_vrpkt) package(pack_index)%typ = type_vdpkt
  ELSE
   package(pack_index)%typ = type_rpkt
   ! if(package(pack_index)%typ == type_rpkt) package(pack_index)%typ = type_rpkt
   ! if(package(pack_index)%typ == type_vrpkt) package(pack_index)%typ = type_vrpkt
  END IF
 end if
END IF
 
END SUBROUTINE change_cell
