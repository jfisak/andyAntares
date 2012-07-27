SUBROUTINE boundary(pack_index, dist, next_cell) 
! write nx_cell, ny_cell, nz_cell as global variable

! Calculate the shortest distance to the cell surface which photon will cross and return this 
! distance (dist), the cell the photon will go to (next_cell, -99 in case of photon leaving 
! the simulation grid(volume)), and make the photon remember which surface will cross (update 
! package(pack_index)%last_cross).

  USE types

  IMPLICIT NONE    

    INTEGER                         :: pack_index, nc, next_cell,forbidden, hit_surface
    DOUBLE PRECISION                :: dist_plusx, dist_minx, dist_plusy, dist_miny, dist_plusz
    DOUBLE PRECISION                :: dist_minz, dist 

    hit_surface = 0
 
!   Number of the current cell
    nc = package(pack_index)%cell_numb

!   Calculate the distances to the all cell surfaces from the photon curent position along the ray
!   (formula for this can be found in http://www.roe.ac.uk/ifa/postgrad/pedagogy/2009_forgan.pdf 
!    - Fig. 2, An Introduction to Monte Carlo Radiative Transfer, Duncan Forgan)
    dist_plusx = (cell(nc)%corner(1) + cell_width - package(pack_index)%pos(1))/package(pack_index)%dir(1)
    dist_minx = (cell(nc)%corner(1) - package(pack_index)%pos(1))/package(pack_index)%dir(1)
    dist_plusy = (cell(nc)%corner(2) + cell_width - package(pack_index)%pos(2))/package(pack_index)%dir(2)
    dist_miny = (cell(nc)%corner(2) - package(pack_index)%pos(2))/package(pack_index)%dir(2)
    dist_plusz = (cell(nc)%corner(3) + cell_width - package(pack_index)%pos(3))/package(pack_index)%dir(3)     
    dist_minz = (cell(nc)%corner(3) - package(pack_index)%pos(3))/package(pack_index)%dir(3)
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
      IF (cell(nc)%indexc(1) .EQ. nx_cell) THEN
!         Photon escapes
          next_cell = -99 
      ELSE
         next_cell = nc + ny_cell * nz_cell
 	 package(pack_index)%last_cross = posx
      END IF
    END IF
    IF (debug .EQ. 1) THEN 
        print*, 'boundary: dist, nextcell, current cell, ny_cell, nz_cell ', &
             dist, next_cell, nc , ny_cell, nz_cell
    END IF

    IF ((dist_minx .GT. 0.D0).AND.(dist_minx .LT. dist).AND.(forbidden .NE. posx)) THEN
       dist = dist_minx
      hit_surface = negx
       IF (cell(nc)%indexc(1) .EQ. 1) THEN
           next_cell = -99
       ELSE
           next_cell = nc - ny_cell * nz_cell
 	   package(pack_index)%last_cross = negx
       END IF
    END IF
    IF (debug .EQ. 1) THEN 
        print*, 'boundary: dist, nextcell ', dist, next_cell
    END IF

    IF ((dist_plusy .GT. 0.D0).AND.(dist_plusy .LT. dist).AND.(forbidden .NE. negy)) THEN
       dist = dist_plusy      
       hit_surface = posy
       IF (cell(nc)%indexc(2) .EQ. ny_cell) THEN
          next_cell = -99
       ELSE
          next_cell = nc + nz_cell
 	  package(pack_index)%last_cross = posy
       END IF
    END IF
    IF (debug .EQ. 1) THEN 
        print*, 'boundary: dist, nextcell ', dist, next_cell
    END IF

    IF ((dist_miny .GT. 0.D0).AND.(dist_miny .LT. dist).AND.(forbidden .NE. posy)) THEN
      dist = dist_miny
      hit_surface = negy
      IF (cell(nc)%indexc(2) .EQ. 1) THEN
          next_cell = -99
      ELSE
         next_cell = nc - nz_cell
 	 package(pack_index)%last_cross = negy
      END IF
    END IF
    IF (debug .EQ. 1) THEN 
        print*, 'boundary: dist, nextcell ', dist, next_cell
    END IF

    IF ((dist_plusz .GT. 0.D0).AND.(dist_plusz .LT. dist).AND.(forbidden .NE. negz)) THEN
      dist = dist_plusz
      hit_surface = posz
      IF (cell(nc)%indexc(3) .EQ. nz_cell) THEN
          next_cell = -99
      ELSE
         next_cell = nc + 1
 	 package(pack_index)%last_cross = posz
      END IF
    END IF
    IF (debug .EQ. 1) THEN 
        print*, 'boundary: dist, nextcell ', dist, next_cell
    END IF

    IF ((dist_minz .GT. 0.D0).AND.(dist_minz .LT. dist).AND.(forbidden .NE. posz)) THEN
      dist = dist_minz
      hit_surface = negz
      IF (cell(nc)%indexc(3) .EQ. 1) THEN
          next_cell = -99
      ELSE
         next_cell = nc - 1
 	 package(pack_index)%last_cross = negz
      END IF
    END IF
    IF (debug .EQ. 1) THEN 
        print*, 'boundary: dist, nextcell ', dist, next_cell
    END IF

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
      print*, nc, cell(nc)%indexc, cell(nc)%corner, cell(nc)%corner+cell_width

      print*, FLOOR(package(pack_index)%pos(1)/cell_width + dble(nx_cell)/2) + 1
      print*, FLOOR(package(pack_index)%pos(2)/cell_width + dble(ny_cell)/2) + 1
      print*, FLOOR(package(pack_index)%pos(3)/cell_width + dble(nz_cell)/2) + 1

       print*, package(pack_index)%dir, package(pack_index)%pos
       print*, dist_plusx, dist_minx, dist_plusy, dist_miny, dist_plusz, dist_minz

       STOP 'ERROR in determining distance to the next cell crossing'
    END IF
   
END 
