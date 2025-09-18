! this sbr will be activated when the plasma in the current propagation cell is optically thick and
! Monte Carlo is not efficient instead of difficult propagation a new packet is radiated at the
! propGrid cell boundary
!
! INPUT: pack_index(INT): index of packet
! OUTPUT: NONE
!
! RETURN point: 1x
!
SUBROUTINE do_dpackage(pack_index)

USE types
USE constants
USE counters
IMPLICIT NONE

INTEGER                                 :: pack_index
INTEGER                                 :: cur_pgi
DOUBLE PRECISION, ALLOCATABLE           :: rates(:)
INTEGER                                 :: dapprox
INTEGER                                 :: dummypack

DOUBLE PRECISION                        :: loc_sum, sum_Z
DOUBLE PRECISION, DIMENSION(const_dimofspace)          :: width, corner, ran_dir
DOUBLE PRECISION, DIMENSION(const_dimofspace)          :: new_width
INTEGER                                 :: next_mgi
DOUBLE PRECISION                        :: ran2, rand

INTEGER                                 :: next_leak, next_cell

LOGICAL                                 :: active

DOUBLE PRECISION                        :: doppler_D, freq
DOUBLE PRECISION, DIMENSION(const_dimofspace)          :: cross_pos
DOUBLE PRECISION, DIMENSION(const_dimofspace)          :: new_dir

DOUBLE PRECISION                        :: area

INTEGER                                 :: pomocna_bunka

LOGICAL                                 :: next_diff

LOGICAL                                 :: procout = .TRUE.

INTEGER                                 :: cur_mgi, get_package_model_index
DOUBLE PRECISION                        :: cur_temp

dapprox = 0

IF(debug == 4) procout = .TRUE.

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

 cur_mgi = get_package_model_index(pack_index)
 cur_temp = model_grid(cur_mgi)%T

 ! write(22,*) corner, width
 
IF(debug == 2) THEN
 CALL find_dyn_cell1(package(pack_index)%pos, pomocna_bunka)
 write(*,*) 'do_rpackage I: pack_index = ', pack_index, ' cur_pgi = ', cur_pgi, ' neigbors = ', dyn_cell(cur_pgi)%neighbor
 write(*,*) 'do_rpackage I: pack_index = ', pack_index, ' bunka = ', pomocna_bunka
 
 write(*,*) 'do_rpackage I: cell starting = ', dyn_cell(cur_pgi)%corner/const_Rsun
 write(*,*) 'do_rpackage I: packet pos = ', package(pack_index)%pos/const_Rsun
 write(*,*) 'do_rpackage I: cell ending = ', (dyn_cell(cur_pgi)%corner + dyn_cell(cur_pgi)%width)/const_Rsun

 write(*,*) 'do_rpackage I: direction = ', package(pack_index)%dir
END IF

 ! x+
 rates(posx) = width(ind_y) * width(ind_z)
 ! x-
 rates(negx) = rates(posx)
 ! y+
 rates(posy) = width(ind_x) * width(ind_z)
 ! y-
 rates(negy) = rates(posy)
 ! x+
 rates(posz) = width(ind_y) * width(ind_x)
 ! x-
 rates(negz) = rates(posz)

 DO WHILE(active)

  loc_sum = rates(posx) + rates(negx) + rates(posy) + rates(negy) + rates(posz) + rates(negz)

  sum_Z = loc_sum

  rand = ran2(idum) * sum_Z
  
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

  ! test if the next cell exists
  IF(next_cell > 0) THEN
   next_mgi = dyn_cell(next_cell)%model_index
  ELSE
   package(pack_index)%cell_numb = next_cell
   package(pack_index)%active = 0
   package(pack_index)%typ = type_escaped
   package(pack_index)%pos = cross_pos
   count_des_esca = count_des_esca + 1
   CALL freq_from_planck(freq, cur_temp)
   package(pack_index)%freq_rf = freq
   CALL random_unitvector2(ran_dir)
   IF(next_leak == posx) THEN
    new_dir = (/ ran_dir(ind_z), ran_dir(ind_x), -ran_dir(ind_y)   /)
   ELSE IF(next_leak == negx) THEN
    new_dir = (/ -ran_dir(ind_z), ran_dir(ind_x), ran_dir(ind_y)   /)
   ELSE IF(next_leak == posy) THEN
    new_dir = (/ ran_dir(ind_x), ran_dir(ind_z), -ran_dir(ind_y) /)
   ELSE IF(next_leak == negy) THEN
    new_dir = (/ ran_dir(ind_x) , -ran_dir(ind_z) , ran_dir(ind_y) /)
   ELSE IF(next_leak == posz) THEN
    new_dir = (/ -ran_dir(ind_y) , ran_dir(ind_x),ran_dir(ind_z)/)
   ELSE IF(next_leak == negz) THEN
    new_dir = (/ ran_dir(ind_y), ran_dir(ind_x), -ran_dir(ind_z)/)
   END IF
   if(procout) write(*,*) 'do_dpackage: a change to r-packet and end...'
   count_d_rad_end = count_d_rad_end + 1
   ! RETURN point
   RETURN
  END IF

  ! packet can be changed into an r-packet or just move to another diffusive cell
  
  IF(next_mgi .ne. photosphere_index .and. next_mgi .ne. vacuum_index .and. &
    next_mgi .ne. outerspace_index) THEN
   next_diff = model_grid(next_mgi)%is_difapp
   ! moving to another diffusive cell
   IF(next_diff) THEN
    package(pack_index)%cell_numb = next_cell
    active = .false.
    if(procout) write(*,*) 'do_dpackage: d-change of cell...'
    count_d_change_cell = count_d_change_cell + 1
   ELSE
    package(pack_index)%pos = cross_pos
    package(pack_index)%cell_numb = next_cell
    package(pack_index)%typ = type_rpkt
    
    CALL random_unitvector2(ran_dir)
    IF(next_leak == posx) THEN
     new_dir = (/ ran_dir(ind_z), ran_dir(ind_x), -ran_dir(ind_y)   /)
    ELSE IF(next_leak == negx) THEN
     new_dir = (/ -ran_dir(ind_z), ran_dir(ind_x), ran_dir(ind_y)   /)
    ELSE IF(next_leak == posy) THEN
     new_dir = (/ ran_dir(ind_x), ran_dir(ind_z), -ran_dir(ind_y) /)
    ELSE IF(next_leak == negy) THEN
     new_dir = (/ ran_dir(ind_x) , -ran_dir(ind_z) , ran_dir(ind_y) /)
    ELSE IF(next_leak == posz) THEN
     new_dir = (/ -ran_dir(ind_y) , ran_dir(ind_x),ran_dir(ind_z)/)
    ELSE IF(next_leak == negz) THEN
     new_dir = (/ ran_dir(ind_y), ran_dir(ind_x), -ran_dir(ind_z)/)
    END IF
    cur_temp = model_grid(cur_mgi)%T
    package(pack_index)%dir = new_dir

    CALL freq_from_planck(freq, cur_temp)
    package(pack_index)%freq_rf = freq
    CALL doppler_factor(pack_index, doppler_D)
    package(pack_index)%freq_cmf = package(pack_index)%freq_rf * doppler_D 
    package(pack_index)%e_cmf    = package(pack_index)%e_rf * doppler_D  
    package(pack_index)%last_line = no_line
    package(pack_index)%next_cross = next_leak

    active = .false.

    if(procout) write(*,*) 'do_dpackage: a radiative change...'
    count_d_radiative = count_d_radiative + 1
   END IF
  ELSE
   ! another choice must be done...
  ! we have to repeat the choice
    new_width = dyn_cell(next_cell)%width
    IF(next_leak == posx) THEN
     area = new_width(ind_y) * new_width(ind_z)
    ELSE IF(next_leak == negx) THEN
     area = new_width(ind_y) * new_width(ind_z)
    ELSE IF(next_leak == posy) THEN
     area = new_width(ind_x) * new_width(ind_z)
    ELSE IF(next_leak == negy) THEN
     area = new_width(ind_x) * new_width(ind_z)
    ELSE IF(next_leak == posz) THEN
     area = new_width(ind_y) * new_width(ind_x)
    ELSE IF(next_leak == negz) THEN
     area = new_width(ind_y) * new_width(ind_x)
    END IF
    
    rates(next_leak) = rates(next_leak) - area
    if(rates(next_leak) < 0.0) STOP 'do_dpackage: rate < 0'

    if(procout) write(*,*) 'do_dpackage: a new choice...'
    count_d_new_choice = count_d_new_choice + 1

  END IF
  
 END DO

CASE DEFAULT
 write(*,*) 'do_dpackage: the choice dapprox = ', dapprox, ' is not known...'
 STOP
END SELECT


END SUBROUTINE do_dpackage
