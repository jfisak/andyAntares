SUBROUTINE move_package(pack_index, dist)

  ! Move photon package from the curent position for some distance (update package(pack_index)%pos)

  USE types
  USE counters

  IMPLICIT NONE    

  INTEGER                           :: pack_index, nc
  DOUBLE PRECISION                  :: dist, length, D, vec_length
  INTEGER                           :: dummypackage
  INTEGER                           :: n_pack_d

!write(*,*)'move_package: Thread rank: ', my_rank
n_pack_d = SIZE(package)
dummypackage = n_pack_d - n_dummy_packs + my_rank + 1
!write(*,*) 'move_package: my_rank = ', my_rank, ' dummypackage = ', dummypackage

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
   ! print*, 'package ', pack_index, ' was destroyed because has come back to the photosphere'
   package(pack_index)%active = 0
   !$OMP ATOMIC
   count_des_phot = count_des_phot + 1
  END IF

  ! Rest frame quantities do not change while propagating without any events, 
  ! but cmf quantities need to be updated
  CALL doppler_factor(pack_index, D)
  package(pack_index)%freq_cmf = package(pack_index)%freq_rf * D
  package(pack_index)%e_cmf = package(pack_index)%e_rf * D
!  print*, 'frequency in frame: ', package(pack_index)%freq_rf, &
!        'frequency in CMF: ', package(pack_index)%freq_cmf
  IF (package(pack_index)%freq_cmf < 0) THEN
   write(*,*) 'move_package: package = ', pack_index, ' prop. cell = ', package(pack_index)%cell_numb
   write(*,*) 'FREQUENCY IS LOWER THAN ZERO!!!'
   !STOP 
  END IF
  package(pack_index)%delta_s = package(pack_index)%delta_s + dist

END SUBROUTINE move_package
