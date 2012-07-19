SUBROUTINE move_package(pack_index, dist)

! Move photon package from the curent position for some distance (abdate package(pack_index)%pos)

  USE types

  IMPLICIT NONE    

    INTEGER                           :: pack_index, nc
    DOUBLE PRECISION                  :: dist, length, D, vec_length
!
    IF (debug .EQ. 1) THEN 
        print*, package(pack_index)%pos, pack_index, dist
    END IF

   IF (debug .EQ. 2) THEN 
      print*, package(pack_index)%pos, pack_index, dist
      nc=package(pack_index)%cell_numb
      print*, nc, cell(nc)%indexc, cell(nc)%corner, cell(nc)%corner+cell_width
      print*, FLOOR(package(pack_index)%pos(1)/cell_width + nx_cell/2) + 1
      print*, FLOOR(package(pack_index)%pos(2)/cell_width + ny_cell/2) + 1
      print*, FLOOR(package(pack_index)%pos(3)/cell_width + nz_cell/2) + 1
    END IF

    package(pack_index)%pos(1) = package(pack_index)%pos(1) + dist * package(pack_index)%dir(1) 
    package(pack_index)%pos(2) = package(pack_index)%pos(2) + dist * package(pack_index)%dir(2) 
    package(pack_index)%pos(3) = package(pack_index)%pos(3) + dist * package(pack_index)%dir(3)

!   Deactivate packets which travel beyond the photosphere
!    length=SQRT(package(pack_index)%pos(1)**2 + package(pack_index)%pos(2)**2 + package(pack_index)%pos(3)**2)  
    IF (vec_length(package(pack_index)%pos) .LT. R_star) package(pack_index)%active = 0


!   Rest frame quantities do not change while propagating without any events, 
!   but cmf quantities need to be updated
    CALL doppler_factor(pack_index, D)
    package(pack_index)%freq_cmf = package(pack_index)%freq_rf * D
    package(pack_index)%e_cmf = package(pack_index)%e_rf * package(pack_index)%freq_cmf/package(pack_index)%freq_rf
   
   IF (debug .EQ. 1) THEN 
       print*, package(pack_index)%pos
   END IF

   IF (debug .EQ. 2) THEN 
      print*, package(pack_index)%pos
      nc=package(pack_index)%cell_numb
      print*, nc, cell(nc)%indexc, cell(nc)%corner, cell(nc)%corner+cell_width

      print*, FLOOR(package(pack_index)%pos(1)/cell_width + nx_cell/2) + 1
      print*, FLOOR(package(pack_index)%pos(2)/cell_width + ny_cell/2) + 1
      print*, FLOOR(package(pack_index)%pos(3)/cell_width + nz_cell/2) + 1
    END IF


!  Check if the package is the inside the simulatio grid 
!   IF ((package(pack_index)%pos(1) .LT. -xmax) .OR. ( package(pack_index)%pos(1) .GT. xmax)) THEN
!       STOP 'ERROR - package is outside of grid of x coordinate' - do this cheking for y and z 
!       coordinates also
!  Transformations law from cmf (update freq and energy) should be added
    
END
