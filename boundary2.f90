SUBROUTINE boundary2(pack_index, dist, next_cell) 
! write nx_cell, ny_cell, nz_cell as global variable

! Calculate the shortest distance to the cell surface which photon will cross and return this 
! distance (dist), the cell the photon will go to (next_cell, -99 in case of photon leaving 
! the simulation grid(volume)), and make the photon remember which surface will cross (update 
! package(pack_index)%last_cross).

  USE types

  IMPLICIT NONE    

    INTEGER                         :: pack_index, nc, next_cell,forbidden, hit_surface
    INTEGER                         :: get_package_model_index
    DOUBLE PRECISION                :: dist_plusx, dist_minx, dist_plusy, dist_miny, dist_plusz
    DOUBLE PRECISION                :: dist_minz, dist 
    DOUBLE PRECISION                :: tTest

    hit_surface = 0
 
!   Number of the current cell
    nc = package(pack_index)%cell_numb

!   Calculate the distances to the all cell surfaces from the photon current position along the ray
!   (formula for this can be found in http://www.roe.ac.uk/ifa/postgrad/pedagogy/2009_forgan.pdf 
!    - Fig. 2, An Introduction to Monte Carlo Radiative Transfer, Duncan Forgan)
   IF(package(pack_index)%dir(1) /= 0) THEN
    dist_plusx = (dyn_cell(nc)%corner(1) + dyn_cell(nc)%width(1) - package(pack_index)%pos(1))/package(pack_index)%dir(1)
    dist_minx = (dyn_cell(nc)%corner(1) - package(pack_index)%pos(1))/package(pack_index)%dir(1)
   ELSE
    dist_plusx = 0.E0
    dist_minx = 0.E0
   END IF
   IF(package(pack_index)%dir(2) /= 0) THEN
    dist_plusy = (dyn_cell(nc)%corner(2) + dyn_cell(nc)%width(2) - package(pack_index)%pos(2))/package(pack_index)%dir(2)
    dist_miny = (dyn_cell(nc)%corner(2) - package(pack_index)%pos(2))/package(pack_index)%dir(2)
   ELSE
    dist_plusx = 0.E0
    dist_minx = 0.E0
   END IF
   IF(package(pack_index)%dir(3) /= 0) THEN
    dist_plusz = (dyn_cell(nc)%corner(3) + dyn_cell(nc)%width(3) - package(pack_index)%pos(3))/package(pack_index)%dir(3)     
    dist_minz = (dyn_cell(nc)%corner(3) - package(pack_index)%pos(3))/package(pack_index)%dir(3)
   ELSE
    dist_plusx = 0.E0
    dist_minx = 0.E0
   END IF
        print*, 'boundary: dir ', package(pack_index)%dir
        print*, 'boundary: distances posx, negx, posy, negy, posz, negz ', &
             dist_plusx, dist_minx, dist_plusy, dist_miny, dist_plusz, dist_minz
    IF (debug .EQ. 1) THEN 
        print*, 'boundary: dir ', package(pack_index)%dir
        print*, 'boundary: distances posx, negx, posy, negy, posz, negz ', &
             dist_plusx, dist_minx, dist_plusy, dist_miny, dist_plusz, dist_minz
    END IF

!   Initalisation - large distance bigger than grid size (all calculated idstances should be shorter)
    dist = 1.D99
!   Remember which surface the photon crossed in its last movement (positive or negative direction and which surface)
    forbidden = package(pack_index)%last_cross

!   Find the shortest of above distances
!
    IF ((dist_plusx .GT. 0.D0).AND.(dist_plusx .LT. dist).AND.(forbidden .NE. negx)) THEN
      dist = dist_plusx
      hit_surface = posx
    END IF
    IF (debug .EQ. 1) THEN 
        print*, 'boundary: dist, nextcell, current cell, ny_cell, nz_cell ', &
             dist, next_cell, nc , ny_cell, nz_cell
    END IF

    IF ((dist_minx .GT. 0.D0).AND.(dist_minx .LT. dist).AND.(forbidden .NE. posx)) THEN
       dist = dist_minx
      hit_surface = negx
    END IF
    IF (debug .EQ. 1) THEN 
        print*, 'boundary: dist, nextcell ', dist, next_cell
    END IF

    IF ((dist_plusy .GT. 0.D0).AND.(dist_plusy .LT. dist).AND.(forbidden .NE. negy)) THEN
       dist = dist_plusy      
       hit_surface = posy
    END IF
    IF (debug .EQ. 1) THEN 
        print*, 'boundary: dist, nextcell ', dist, next_cell
    END IF

    IF ((dist_miny .GT. 0.D0).AND.(dist_miny .LT. dist).AND.(forbidden .NE. posy)) THEN
      dist = dist_miny
      hit_surface = negy
    END IF
    IF (debug .EQ. 1) THEN 
        print*, 'boundary: dist, nextcell ', dist, next_cell
    END IF

    IF ((dist_plusz .GT. 0.D0).AND.(dist_plusz .LT. dist).AND.(forbidden .NE. negz)) THEN
      dist = dist_plusz
      hit_surface = posz
    END IF
    IF (debug .EQ. 1) THEN 
        print*, 'boundary: dist, nextcell ', dist, next_cell
    END IF

    IF ((dist_minz .GT. 0.D0).AND.(dist_minz .LT. dist).AND.(forbidden .NE. posz)) THEN
      dist = dist_minz
      hit_surface = negz
    END IF
    IF (debug .EQ. 1) THEN 
        print*, 'boundary: dist, nextcell ', dist, next_cell
    END IF
   ! finding testing t
   tTest = dist + minwidth/1.D8
   IF((dist_minx > 0) .AND. (tTest > dist_minx)) tTest = dist + dist_minx/2.D0
   IF((dist_plusx > 0) .AND. (tTest > dist_plusx)) tTest = dist + dist_plusx/2.D0
   IF((dist_miny > 0) .AND. (tTest > dist_miny)) tTest = dist + dist_miny/2.D0
   IF((dist_plusy > 0) .AND. (tTest > dist_plusy)) tTest = dist + dist_plusy/2.D0
   IF((dist_minz > 0) .AND. (tTest > dist_minz)) tTest = dist + dist_minz/2.D0
   IF((dist_plusz > 0) .AND. (tTest > dist_plusz)) tTest = dist + dist_plusz/2.D0
    ! now will calculate next cell
    CALL next_dyn_cell(pack_index,dist,tTest,next_cell)
    IF (debug .EQ. 3) THEN
      print*, nc, cell(nc)%indexc, &
           FLOOR(package(pack_index)%pos(1)/cell_width + dble(nx_cell)/2) + 1, &
           FLOOR(package(pack_index)%pos(2)/cell_width + dble(ny_cell)/2) + 1, &
           FLOOR(package(pack_index)%pos(3)/cell_width + dble(nz_cell)/2) + 1 
      print*, package(pack_index)%pos, package(pack_index)%dir
      print*, cell(nc)%corner, cell(nc)%corner+cell_width
      print*, dist_plusx, dist_minx, dist_plusy, dist_miny, dist_plusz, dist_minz
      print*, dist, next_cell, hit_surface, forbidden       
   ENDIF



    IF (hit_surface .EQ. 0) THEN
      print*, pack_index
      print*, dist, next_cell, hit_surface       
      print*, forbidden
      print*, nc, dyn_cell(nc)%corner, dyn_cell(nc)%corner+ dyn_cell(nc)%width

      print*, FLOOR(package(pack_index)%pos(1)/cell_width + dble(nx_cell)/2) + 1
      print*, FLOOR(package(pack_index)%pos(2)/cell_width + dble(ny_cell)/2) + 1
      print*, FLOOR(package(pack_index)%pos(3)/cell_width + dble(nz_cell)/2) + 1

       print*, package(pack_index)%dir, package(pack_index)%pos
       print*, dist_plusx, dist_minx, dist_plusy, dist_miny, dist_plusz, dist_minz

       STOP 'ERROR in determining distance to the next cell crossing'
    END IF
   
END SUBROUTINE boundary2
