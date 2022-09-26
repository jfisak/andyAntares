! this sbr will be activated when the plasma in the current propagation cell is optically thick and
! Monte Carlo is not efficient instead of difficult propagation a new packet is radiated at the
! propGrid cell boundary
SUBROUTINE do_dpackage(pack_index)

USE types
IMPLICIT NONE

INTEGER                                 :: pack_index
INTEGER                                 :: cur_pgi
DOUBLE PRECISION, ALLOCATABLE           :: rates(:)
INTEGER                                 :: dapprox
INTEGER                                 :: dummypack, I

DOUBLE PRECISION                        :: loc_sum, Z
DOUBLE PRECISION, DIMENSION(3)          :: width, corner
INTEGER, DIMENSION(3)                   :: cur_dir
INTEGER                                 :: next_mgi
DOUBLE PRECISION                        :: ran2, rand

DOUBLE PRECISION                        :: dist
INTEGER                                 :: next_leak, next_cell

LOGICAL                                 :: active

dapprox = 0

dummypack = SIZE(package)
package(dummypack) = package(pack_index)
! cur_pgi = package(pack_index)%cell_numb
! width = dyn_cell(cur_pgi)%width

! the most stupid approximation: a cell surface is chosen and packet will move to the neighboring cell
! if the other cell is also diffusion this machinery continues
! if the other cell is normal a packet is radiated with a Planck law
SELECT CASE (dapprox)
CASE(0)
 loc_sum = 0.E0
 ALLOCATE(rates(6))

 active = .true.

 DO
  cur_pgi = package(pack_index)%cell_numb
  corner = dyn_cell(cur_pgi)%corner
  width = dyn_cell(cur_pgi)%width

  ! x+
  rates(posx) = width(2) * width(3)
  ! x-
  rates(negx) = rates(1)
  ! y+
  rates(posy) = width(1) * width(3)
  ! y-
  rates(negy) = rates(3)
  ! x+
  rates(posz) = width(2) * width(1)
  ! x-
  rates(negz) = rates(5)

  loc_sum = rates(posx) + rates(negx) + rates(posy) + rates(negy) + rates(posz) + rates(negz)

  Z = loc_sum

  rand = ran2(idum) * Z
  
  loc_sum = 0.E0
  IF(rand > loc_sum .and. rand < loc_sum + rates(posx)) THEN
   next_leak = posx
  END IF
  loc_sum = loc_sum + rates(posx)
  IF(rand > loc_sum .and. rand < loc_sum + rates(negx)) THEN
   next_leak = negx
  END IF
  loc_sum = loc_sum + rates(negx)
  IF(rand > loc_sum .and. rand < loc_sum + rates(posy)) THEN
   next_leak = posy
  END IF
  loc_sum = loc_sum + rates(posy)
  IF(rand > loc_sum .and. rand < loc_sum + rates(negy)) THEN
   next_leak = negy
  END IF
  loc_sum = loc_sum + rates(negy)
  IF(rand > loc_sum .and. rand < loc_sum + rates(posz)) THEN
   next_leak = posz
  END IF
  loc_sum = loc_sum + rates(posz)
  IF(rand > loc_sum .and. rand < loc_sum + rates(negz)) THEN
   next_leak = negz
  END IF

  ! random position
  DO I = 1,3
   package(dummypack)%pos(I) = corner(I) + ran2(idum) * width(I)
  END DO
  

  IF(next_leak == posx) THEN
   cur_dir = (/1, 0, 0 /)
   package(dummypack)%dir = cur_dir
   package(dummypack)%pos(1) = corner(1) + width(1)/2.0
  ELSE IF(next_leak == negx) THEN
   cur_dir = (/-1, 0, 0 /)
   package(dummypack)%dir = cur_dir
   package(dummypack)%pos(1) = corner(1) + width(1)/2.0
  ELSE IF(next_leak == posy) THEN
   cur_dir = (/0, 1, 0 /)
   package(dummypack)%dir = cur_dir
   package(dummypack)%pos(2) = corner(2) + width(2)/2.0
  ELSE IF(next_leak == negy) THEN
   cur_dir = (/0, -1, 0 /)
   package(dummypack)%dir = cur_dir
   package(dummypack)%pos(2) = corner(2) + width(2)/2.0
  ELSE IF(next_leak == posz) THEN
   cur_dir = (/0, 0, 1 /)
   package(dummypack)%dir = cur_dir
   package(dummypack)%pos(3) = corner(3) + width(3)/2.0
  ELSE IF(next_leak == negz) THEN
   cur_dir = (/0, 0, -1 /)
   package(dummypack)%dir = cur_dir
   package(dummypack)%pos(3) = corner(3) + width(3)/2.0
  END IF

  ! what is the next cell in this configuration?
  CALL boundary3(dummypack, dist, next_cell)

  next_mgi = dyn_cell(next_cell)%model_index
  
  IF(next_mgi < n_modelgrid + 1) THEN
   CALL change_cell(pack_index, next_cell)
   active = .false.
   ! we will set up the properties of the r-packet if the next propagation cell is not diffussive
  END IF

 END DO

 










CASE DEFAULT
 write(*,*) 'do_dpackage: the choice dapprox = ', dapprox, ' is not known...'
 STOP
END SELECT

END SUBROUTINE do_dpackage
