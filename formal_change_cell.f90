SUBROUTINE formal_change_cell(pack_index, next_cell)

USE types
IMPLICIT NONE


INTEGER                                         :: pack_index, next_cross, next_cell
DOUBLE PRECISION, DIMENSION(3)                  :: direction
INTEGER, PARAMETER                              :: dir_x = 1, dir_y = 2, dir_z = 3
INTEGER                                         :: cur_dir, cur_pgi, cur_n1, cur_n2, ind_neighbor
LOGICAL                                         :: x_neg = .FALSE., y_neg = .FALSE., z_neg = .FALSE.
LOGICAL                                         :: cross_1, cross_2
DOUBLE PRECISION, DIMENSION(3)                  :: cur_n

cur_n = package(pack_index)%dir
write(*,*) 'formal_change_cell: pack_index = ', pack_index
write(*,*) 'formal_change_cell: cur_n = ', cur_n
cur_pgi = package(pack_index)%cell_numb
next_cross = package(pack_index)%next_cross
cross_1 = .false.
cross_2 = .false.


! get the direction of the zero corner
! the following rules are applied
! 1.) x-direction
!     n1 = dir_y
!     n2 = dir_z
! 2.) y-direction
!     n1 = dir_x
!     n2 = dir_z
! 3.) z-direction
!     n1 = dir_y
!     n2 = dir_z
!____________________________________
! 1.) x-direction
IF((next_cross >= 19) .and. (next_cross <= 24)) THEN
 cur_dir = dir_x
 cur_n1 = dir_y
 cur_n2 = dir_z
END IF

! 2.) y-direction
IF((next_cross >= 13) .and. (next_cross <= 18)) THEN
 cur_n1 = dir_x
 cur_dir = dir_y
 cur_n2 = dir_z
END IF

! 3.) z-direction
IF((next_cross >= 7) .and. (next_cross <= 12)) THEN
 cur_n1 = dir_x
 cur_n2 = dir_y
 cur_dir = dir_z
END IF


! get the negative components of the direction vector n
IF(cur_n(1) < 0) x_neg = .TRUE.
IF(cur_n(2) < 0) y_neg = .TRUE.
IF(cur_n(3) < 0) z_neg = .TRUE.

write(*,*) 'formal_change_cell: cur_dir = ', cur_dir, ' dir_x = ', dir_x

! calculation of the crosses we are interested in
IF(cur_dir == dir_x) THEN
 IF(y_neg) cross_1 = .TRUE.
 ind_neighbor = dyn_cell(cur_pgi)%neighbor(negy)
 IF(z_neg) cross_2 = .TRUE.
 ind_neighbor = dyn_cell(cur_pgi)%neighbor(negz)
 IF(y_neg .and. z_neg) THEN
  ind_neighbor = dyn_cell(cur_pgi)%neighbor(negy)
  ind_neighbor = dyn_cell(ind_neighbor)%neighbor(negz)
 END IF
 write(*,*) 'formal_change_cell: ind_neighbor = ', ind_neighbor

 

! zatím budu testovat pouze pro x-ovou osu
! y a z budou zatím spát
! ELSE IF(cur_dir == dir_y) THEN
!  IF(x_neg) cross_1 = .TRUE.
!  ind_neighbor = dyn_cell(cur_pgi)%neighbor(xneg)
!  IF(z_neg) cross_2 = .TRUE.
!  ind_neighbor = dyn_cell(cur_pgi)%neighbor(zneg)
! ELSE IF(cur_dir == dir_z) THEN
!  IF(x_neg) cross_2 = .TRUE.
!  ind_neighbor = dyn_cell(cur_pgi)%neighbor(xneg)
!  IF(y_neg) cross_1 = .TRUE.
!  ind_neighbor = dyn_cell(cur_pgi)%neighbor(yneg)
END IF

! we are not interested in the first quadant, no mistake is presented
write(*,*) 'formal_change_cell: cross_1 = ', cross_1, ' cross_2 = ', cross_2
write(*,*) 'formal_change_cell: (.not. cross_1) = ', (.NOT. cross_1)
write(*,*) 'formal_change_cell: ((.not. cross_1) .and. (.not. cross_2)) = ', ((.not. cross_1) .and. (.not. cross_2))
IF((.not. cross_1) .and. (.not. cross_2)) THEN
 IF((next_cross == edxy_xn) .or. (next_cross == edxz_xn) .or. (next_cross == edyz_xn)) THEN
  next_cross = negx
 ELSE IF((next_cross == edxy_xp) .or. (next_cross == edxz_xp) .or. (next_cross == edyz_xp)) THEN
  next_cross = posx
 ELSE IF((next_cross == edxy_yn) .or. (next_cross == edyz_yn) .or. (next_cross == edxz_yn)) THEN
  next_cross = negy
 ELSE IF((next_cross == edxy_yp) .or. (next_cross == edyz_yp) .or. (next_cross == edxz_yp)) THEN
  next_cross = posy
 ELSE IF((next_cross == edxz_zn) .or. (next_cross == edyz_zn) .or. (next_cross == edxy_zn)) THEN
  next_cross = negz
 ELSE IF((next_cross == edxz_zp) .or. (next_cross == edyz_zp) .or. (next_cross == edxy_zn)) THEN
  next_cross = posz
 END IF
 next_cell = dyn_cell(cur_pgi)%neighbor(next_cross)
 RETURN
END IF

! therefore we look on the neighbours of the propGrid cell

next_cell = ind_neighbor
write(*,*) 'formal_change_cell: a new index of the propGrid cell is: ', ind_neighbor










END SUBROUTINE formal_change_cell
