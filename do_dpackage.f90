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
DOUBLE PRECISION, DIMENSION(3)          :: width, corner, ran_dir, pack_dir, cur_center
INTEGER, DIMENSION(3)                   :: cur_dir
INTEGER                                 :: next_mgi
DOUBLE PRECISION                        :: ran2, rand

DOUBLE PRECISION                        :: dist
INTEGER                                 :: next_leak, next_cell

LOGICAL                                 :: active

DOUBLE PRECISION                        :: D, freq
DOUBLE PRECISION, DIMENSION(3)          :: cross_pos
DOUBLE PRECISION, DIMENSION(3)          :: new_dir

INTEGER                                 :: pomocna_bunka

LOGICAL                                 :: next_diff

dapprox = 0

dummypack = SIZE(package)
package(dummypack) = package(pack_index)
! width = dyn_cell(cur_pgi)%width

! the most stupid approximation: a cell surface is chosen and packet will move to the neighboring cell
! if the other cell is also diffusion this machinery continues
! if the other cell is normal a packet is radiated with a Planck law
SELECT CASE (dapprox)
CASE(0)
 loc_sum = 0.E0
 ALLOCATE(rates(6))

 active = .true.

 cur_pgi = package(pack_index)%cell_numb
 corner = dyn_cell(cur_pgi)%corner
 width = dyn_cell(cur_pgi)%width
 
IF(debug == 3) THEN
 CALL find_dyn_cell1(package(pack_index)%pos, pomocna_bunka)
 write(*,*) 'do_rpackage I: pack_index = ', pack_index, ' cur_pgi = ', cur_pgi, ' neigbors = ', dyn_cell(cur_pgi)%neighbor
 write(*,*) 'do_rpackage I: pack_index = ', pack_index, ' bunka = ', pomocna_bunka
 
 write(*,*) 'do_rpackage I: cell starting = ', dyn_cell(cur_pgi)%corner/R_sun
 write(*,*) 'do_rpackage I: packet pos = ', package(pack_index)%pos/R_sun
 write(*,*) 'do_rpackage I: cell ending = ', (dyn_cell(cur_pgi)%corner + dyn_cell(cur_pgi)%width)/R_sun

 write(*,*) 'do_rpackage I: direction = ', package(pack_index)%dir
END IF

 ! x+
 rates(posx) = width(2) * width(3)
 ! x-
 rates(negx) = rates(posx)
 ! y+
 rates(posy) = width(1) * width(3)
 ! y-
 rates(negy) = rates(posy)
 ! x+
 rates(posz) = width(2) * width(1)
 ! x-
 rates(negz) = rates(posz)

 DO WHILE(active)

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

  CALL d_choosenextcell(cur_pgi, next_leak, next_cell, cross_pos)

  next_mgi = dyn_cell(next_cell)%model_index
  ! packet can be changed into an r-packet
  IF(next_mgi < n_modelgrid + 1) THEN
   next_diff = model_grid(next_mgi)%is_difapp
   IF(next_diff) THEN
    package(pack_index)%cell_numb = next_cell
    active = .false.
   ELSE
    package(pack_index)%pos = cross_pos
    write(23,*) cross_pos
    package(pack_index)%cell_numb = next_cell
    package(pack_index)%typ = type_rpkt
    
    CALL random_unitvector2(ran_dir)
    IF(next_leak == posx) THEN
     new_dir = (/ ran_dir(3), ran_dir(1), -ran_dir(2)   /)
    ELSE IF(next_leak == negx) THEN
     new_dir = (/ -ran_dir(3), ran_dir(1), ran_dir(2)   /)
    ELSE IF(next_leak == posy) THEN
     new_dir = (/ ran_dir(1), ran_dir(3), -ran_dir(2) /)
    ELSE IF(next_leak == negy) THEN
     new_dir = (/ ran_dir(1) , -ran_dir(3) , ran_dir(2) /)
    ELSE IF(next_leak == posz) THEN
     new_dir = (/ -ran_dir(2) , ran_dir(1),ran_dir(3)/)
    ELSE IF(next_leak == negz) THEN
     new_dir = (/ ran_dir(2), ran_dir(1), -ran_dir(3)/)
    END IF
    package(pack_index)%dir = new_dir

    CALL freq_from_planck(freq)
    package(pack_index)%freq_rf = freq
    CALL doppler_factor(pack_index, D)
    package(pack_index)%freq_cmf = package(pack_index)%freq_rf * D 
    package(pack_index)%e_cmf    = package(pack_index)%e_rf * D  
    package(pack_index)%last_line = no_line
    package(pack_index)%next_cross = next_leak

    active = .false.
   END IF
  ! we have to repeat the choice
  ELSE
   ! another choice must be done...
   write(*,*) 'another choice must be done'
  END IF
  
 END DO

write(*,*) 'do_dpackage: pack_index = ', pack_index, ' moving to next cell = ', next_cell
cur_center = corner + width/2.0
write(22,*) cur_center, dyn_cell(next_cell)%corner, dyn_cell(next_cell)%width

CASE DEFAULT
 write(*,*) 'do_dpackage: the choice dapprox = ', dapprox, ' is not known...'
 STOP
END SELECT

END SUBROUTINE do_dpackage
