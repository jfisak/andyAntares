SUBROUTINE boundary3(pack_index, dist, next_cell) 
! write nx_cell, ny_cell, nz_cell as global variable

! Calculate the shortest distance to the cell surface which photon will cross and return this 
! distance (dist), the cell the photon will go to (next_cell, -99 in case of photon leaving 
! the simulation grid(volume)), and make the photon remember which surface will cross (update 
! package(pack_index)%last_cross).

  USE types

  IMPLICIT NONE    

    INTEGER                         :: pack_index, nc, next_cell,forbidden, hit_surface
    DOUBLE PRECISION                :: dist_plusx, dist_minx, dist_plusy, dist_miny, dist_plusz
    DOUBLE PRECISION                :: dist_minz, dist, surface_pos 

    hit_surface = 0
 
!   Number of the current cell
    nc = package(pack_index)%cell_numb

    CALL next_dyn_cell3(nc,pack_index,next_cell,dist)
    !print*, 'nc = ', nc, 'pack_index = ', pack_index, 'next_cell = ', next_cell, 'dist = ', dist
!   Calculate the distances to the all cell surfaces from the photon current position along the ray
!   (formula for this can be found in http://www.roe.ac.uk/ifa/postgrad/pedagogy/2009_forgan.pdf 
!    - Fig. 2, An Introduction to Monte Carlo Radiative Transfer, Duncan Forgan)
    IF (debug .EQ. 1) THEN 
        print*, 'boundary: dir ', package(pack_index)%dir
        print*, 'boundary: distances posx, negx, posy, negy, posz, negz ', &
             dist_plusx, dist_minx, dist_plusy, dist_miny, dist_plusz, dist_minz
    END IF
   
END SUBROUTINE boundary3
