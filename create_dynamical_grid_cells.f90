! this subroutine will create a dynamical grid cell for the given basic cell
  SUBROUTINE create_dynamical_grid_cells(n_dyncell,max_n_dcell)

   USE types

   IMPLICIT NONE
   
   INTEGER                              :: n_dyncell
   ! number of created dynamic cells and
   ! actual number of grid cell
   INTEGER                              :: max_n_dcell, act_n_dyncell
   INTEGER                              :: I,J
   INTEGER                              :: Npart, np
   ! number of new created cells in cell
   INTEGER                              :: no_dcells = 8
   ! variables for boundaries
   INTEGER                              :: up_bound, bound, newbound
   ! maximal number of particles in one cell
   INTEGER, PARAMETER                   :: maxPart = 2
   DOUBLE PRECISION, DIMENSION(3)       :: corner, cell_width_2
   TYPE(virt_particle), ALLOCATABLE     :: local_particle(:), pom(:)
   TYPE(dyn_grid_cell), ALLOCATABLE     ::  pom2(:)
   ! local variables
   DOUBLE PRECISION, DIMENSION(3)       :: loc_corner, loc_cell_width
   INTEGER                              :: loc_np
   INTEGER                              :: loc_downcell, loc_upcell
    corner(1) = dyn_cell(n_dyncell)%corner(1)
    corner(2) = dyn_cell(n_dyncell)%corner(2)
    corner(3) = dyn_cell(n_dyncell)%corner(3)
    cell_width_2(1) = dyn_cell(n_dyncell)%width(1)
    cell_width_2(2) = dyn_cell(n_dyncell)%width(2)
    cell_width_2(3) = dyn_cell(n_dyncell)%width(3)
    ! at first we have to know, how many particles are
    ! in the given cell
     np = 0
     up_bound = 100
     ! temporary solution
     Npart = SIZE(virtual_particle)
     ALLOCATE(local_particle(up_bound))
    DO J = 1, Npart
     ! is the virtual particle in this cell?
     IF((virtual_particle(J)%pos(1) >= corner(1)) .AND. &
        (virtual_particle(J)%pos(1) < (corner(1) + cell_width_2(1))) .AND. &
        (virtual_particle(J)%pos(2) >= corner(2)) .AND. &
        (virtual_particle(J)%pos(2) < (corner(2) + cell_width_2(2))) .AND. &
        (virtual_particle(J)%pos(3) >= corner(3)) .AND. &
        (virtual_particle(J)%pos(3) < (corner(3) + cell_width_2(3)))) THEN
     ! print*, 'we have a new catched virtual particle :-)'
      np = np + 1
      ! if the field is full, we will have to increase its size
      if(np == up_bound) then
      ! print*, 'creating a larger array local_particle'
       ! define a new upper bound
       newbound = 2 * up_bound
       ALLOCATE(pom(up_bound))
       do I = 1, up_bound
        pom(I) = local_particle(I)
       end do
       DEALLOCATE(local_particle)
       ALLOCATE(local_particle(newbound))
       do I = 1, up_bound
        local_particle(I) = pom(I)
       end do
       DEALLOCATE(pom)
       up_bound = newbound
      end if
      local_particle(np) = virtual_particle(J)
     END IF
    END DO
      print*, 'number of particles: ', np
    ! we have found particles included in this basic cell
    ! now we have to generate brand new dynamical cell
    ! if the number of particles np is equal to one or two
    ! we don't have to allocate any new dynamical cells
    ! number of the dynamic grid cell
    ! firstly for the basic cell division
    !n_dg = 0
    up_bound = SIZE(dyn_cell(:))
    act_n_dyncell = n_dyncell
    DO
    ! local number of particles is in the begining of cycle = 0
     loc_np = 0
      loc_corner(1) = dyn_cell(act_n_dyncell)%corner(1)
      loc_corner(2) = dyn_cell(act_n_dyncell)%corner(2)
      loc_corner(3) = dyn_cell(act_n_dyncell)%corner(3)
      loc_cell_width(1) = dyn_cell(act_n_dyncell)%width(1)
      loc_cell_width(2) = dyn_cell(act_n_dyncell)%width(3)
      loc_cell_width(3) = dyn_cell(act_n_dyncell)%width(2)
      loc_upcell = dyn_cell(act_n_dyncell)%up_cell
      loc_downcell = dyn_cell(act_n_dyncell)%down_cell
     ! how many virtual particles is there in this subcell
     DO J = 1, np
      IF((local_particle(J)%pos(1) .GE. loc_corner(1)) .AND. &
         (local_particle(J)%pos(1) < (loc_corner(1) + loc_cell_width(1))) .AND. &
         (local_particle(J)%pos(2) .GE. loc_corner(2)) .AND. &
         (local_particle(J)%pos(2) < (loc_corner(2) + loc_cell_width(2))) .AND. &
         (local_particle(J)%pos(3) .GE. loc_corner(3)) .AND. &
         (local_particle(J)%pos(3) < (loc_corner(3) + loc_cell_width(3)))) THEN
       loc_np = loc_np + 1
      END IF
     END DO
     ! now we have to decide what to do on the basement of number of
     ! local particles
     ! 1. the division of the cells is good enough
     ! print*, 'dynamic grid: ', up_bound, act_n_dyncell, max_n_dcell
     IF(loc_np <= maxPart) THEN
      ! we move to the lower level of the grid to the "not ending" subcell
      DO
!       IF(loc_downcell /= 0) &
!         print*, 'dyn_cell(down_cell)%up_cell - act_n_dyncell == no_dcells - 1', &
!         act_n_dyncell - dyn_cell(loc_downcell)%up_cell 
       IF(dyn_cell(act_n_dyncell)%down_cell /= 0) THEN
         IF(act_n_dyncell - dyn_cell(loc_downcell)%up_cell == no_dcells - 1) &
          THEN
          ! we will move one level lower
          !print*, 'moving to the lower level...'
          act_n_dyncell = dyn_cell(act_n_dyncell)%down_cell
          loc_upcell = dyn_cell(act_n_dyncell)%up_cell
          loc_downcell = dyn_cell(act_n_dyncell)%down_cell
         ! we can go to the next subcell in the given level
         ELSE
         ! print*, 'moving to the next cell...'
          act_n_dyncell = act_n_dyncell + 1
          EXIT
         END IF
       ELSE
        EXIT
       END IF
      END DO
     ! we have to create additional cells otherwise
     ELSE
      ! we check if the width of new cells will be large enough
      if (dyn_cell(act_n_dyncell)%width(1)/2.E0 < minwidth .OR. &
          dyn_cell(act_n_dyncell)%width(2)/2.E0 < minwidth .OR. &
          dyn_cell(act_n_dyncell)%width(3)/2.E0 < minwidth) then
          loc_np = 0
          print*, 'CELL WOULD BE TOO SMALL...MOVING TO THE NEXT CELL...'
          CYCLE
      end if
     ! is there some free space left in the field dyn_cell?
      up_bound = SIZE(dyn_cell(:))
      if(max_n_dcell + no_dcells >= up_bound) then
       print*, 'creating a larger array dyn_cell...'
       ! define a new upper bound
       newbound = 2 * up_bound
       ALLOCATE(pom2(up_bound))
       do I = 1, up_bound
        pom2(I) = dyn_cell(I)
       end do
       DEALLOCATE(dyn_cell)
       ALLOCATE(dyn_cell(newbound))
       do I = 1, up_bound
        dyn_cell(I) = pom2(I)
       end do
       DEALLOCATE(pom2)
       up_bound = newbound
      end if
      dyn_cell(act_n_dyncell)%up_cell = max_n_dcell + 1
      ! subroutine which creates subcells in the given cell
      CALL create_subcells(max_n_dcell, act_n_dyncell, no_dcells)
      ! we have to increase the variable max_n_dcell
      max_n_dcell = max_n_dcell + no_dcells
      act_n_dyncell = dyn_cell(act_n_dyncell)%up_cell
     END IF
     ! if every single cell is divided correctly the loop will stop
       IF(act_n_dyncell == n_dyncell) THEN
        print*, 'all dynamic cells for this basic cell were created...'
        EXIT
       END IF
    END DO
   CLOSE(15)
  END SUBROUTINE create_dynamical_grid_cells
