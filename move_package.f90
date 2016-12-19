SUBROUTINE move_package(pack_index, dist)

  ! Move photon package from the curent position for some distance (update package(pack_index)%pos)

  USE types

  IMPLICIT NONE    

  INTEGER                           :: pack_index, nc
  DOUBLE PRECISION                  :: dist, length, D, vec_length
!
  IF (debug .EQ. 1) THEN 
    ! print*, 'before move (pos, pack_index, pos)', package(pack_index)%pos, pack_index, dist
  END IF
!     print*, 'before move (pos, pack_index, pos)', package(pack_index)%pos, pack_index, dist

  !print*, 'moving package #', pack_index
  IF (debug .EQ. 2) THEN 
     print*, package(pack_index)%pos, pack_index, dist
   ! nc=package(pack_index)%cell_numb
   ! print*, nc, cell(nc)%indexc, cell(nc)%corner, cell(nc)%corner+cell_width
   ! print*, FLOOR(package(pack_index)%pos(1)/cell_width + nx_cell/2) + 1
   ! print*, FLOOR(package(pack_index)%pos(2)/cell_width + ny_cell/2) + 1
   ! print*, FLOOR(package(pack_index)%pos(3)/cell_width + nz_cell/2) + 1
  END IF

!  print*, 'move_package: moving package for dist = ', dist
  ! Calculate the position of package
  package(pack_index)%pos(1) = package(pack_index)%pos(1) + dist * package(pack_index)%dir(1) 
  package(pack_index)%pos(2) = package(pack_index)%pos(2) + dist * package(pack_index)%dir(2) 
  package(pack_index)%pos(3) = package(pack_index)%pos(3) + dist * package(pack_index)%dir(3)
  if(abs(package(pack_index)%pos(1)) < 1e-1) package(pack_index)%pos(1) = 0e0
  if(abs(package(pack_index)%pos(2)) < 1e-1) package(pack_index)%pos(2) = 0e0
  if(abs(package(pack_index)%pos(3)) < 1e-1) package(pack_index)%pos(3) = 0e0

  ! Deactivate packets which travel beyond the photosphere
  ! length=SQRT(package(pack_index)%pos(1)**2 + package(pack_index)%pos(2)**2 + package(pack_index)%pos(3)**2)  
  IF ((vec_length(package(pack_index)%pos) .LT. R_star) .AND. (pack_index .NE. dummypackage)) THEN
!      print*, 'package ', pack_index, ' was destroyed'
      package(pack_index)%active = 0
      destroyed_pack = destroyed_pack + 1
  END IF

  ! Rest frame quantities do not change while propagating without any events, 
  ! but cmf quantities need to be updated
  CALL doppler_factor(pack_index, D)
  package(pack_index)%freq_cmf = package(pack_index)%freq_rf * D
  package(pack_index)%e_cmf = package(pack_index)%e_rf * D
!  print*, 'frequency in frame: ', package(pack_index)%freq_rf, &
!        'frequency in CMF: ', package(pack_index)%freq_cmf
  IF (package(pack_index)%freq_cmf < 0) STOP 'FREQUENCY IS LOWER THAN ZERO!!!'
  IF (debug .EQ. 1) THEN 
     print*, 'after move (pos)', package(pack_index)%pos
  END IF
!     print*, 'after move (pos)', package(pack_index)%pos

!  IF (debug .EQ. 2) THEN 
!     print*, package(pack_index)%pos
!     nc=package(pack_index)%cell_numb
!     print*, nc, cell(nc)%indexc, cell(nc)%corner, cell(nc)%corner+basic_cell_width
!     
!     print*, FLOOR(package(pack_index)%pos(1)/cell_width + dble(nx_cell)/2) + 1
!     print*, FLOOR(package(pack_index)%pos(2)/cell_width + dble(ny_cell)/2) + 1
!     print*, FLOOR(package(pack_index)%pos(3)/cell_width + dble(nz_cell)/2) + 1
!  END IF

  package(pack_index)%delta_s = package(pack_index)%delta_s + dist

END SUBROUTINE move_package
