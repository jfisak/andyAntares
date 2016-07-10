! this subroutine will create a dynamical grid cell for the given basic cell
  SUBROUTINE create_dynamical_grid_cell(n_dyncell,max_n_dcell)

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
   TYPE(dyn_cell), ALLOCATABLE     ::  pom2(:)
   ! local variables
   DOUBLE PRECISION, DIMENSION(3)       :: loc_corner, loc_cell_width
   INTEGER                              :: loc_np
   ! we have to go through every basic propagation grid cell
   OPEN(15,FILE='nop.dat')
    corner = cell(I)%corner
    cell_width_2(1) = 2.E0 * xmax / nx_cell
    cell_width_2(2) = 2.E0 * ymax / nx_cell
    cell_width_2(3) = 2.E0 * zmax / nx_cell
    Npart = SIZE(virtual_particle(:))
    ! at first we have to know, how many particles are
    ! in the given cell
     np = 0
    DO J = 1, Npart
     ! is the virtual particle in this cell?
     IF((virtual_particle(J)%pos(1) .GE. corner(1)) .AND. &
        (virtual_particle(J)%pos(1) < (corner(1) + cell_width_2(1))) .AND. &
        (virtual_particle(J)%pos(2) .GE. corner(2)) .AND. &
        (virtual_particle(J)%pos(2) < (corner(2) + cell_width_2(2))) .AND. &
        (virtual_particle(J)%pos(3) .GE. corner(3)) .AND. &
        (virtual_particle(J)%pos(3) < (corner(3) + cell_width_2(3)))) THEN
      np = np + 1
      if(np == 0) ALLOCATE(local_particle(up_bound))
      ! if the field is full, we will have to increase its size
      if(np == up_bound) then
       ! define a new upper bound
       newbound = 2.E0 * up_bound
       ALLOCATE(pom(up_bound))
       do I = 1, up_bound
        pom(I) = local_particle(I)
       end do
       DEALLOCATE(local_particle(:))
       ALLOCATE(local_particle(newbound))
       do I = 1, up_bound
        local_particle(I) = pom(I)
       end do
       DEALLOCATE(pom(:))
       up_bound = newbound
      end if
      local_particle(np) = virtual_particle(J)
     END IF
    END DO
    ! we have found particles included in this basic cell
    ! now we have to generate brand new dynamical cell
    ! if the number of particles np is equal to one or two
    ! we don't have to allocate any new dynamical cells
    ! number of the dynamic grid cell
    ! firstly for the basic cell division
    !n_dg = 0
    up_bound = SIZE(dyn_cell(:))
    DO
     loc_np = 0
      loc_corner(1) = dyn_cell(n_dyncell)%corner(1)
      loc_corner(2) = dyn_cell(n_dyncell)%corner(2)
      loc_corner(3) = dyn_cell(n_dyncell)%corner(3)
      loc_cell_width(1) = dyn_cell(n_dyncell)%width(1)
      loc_cell_width(2) = dyn_cell(n_dyncell)%width(3)
      loc_cell_width(3) = dyn_cell(n_dyncell)%width(2)
     ! how many particles is there in this subcell
     DO J = 1, np
      ! is the virtual particle in this cell?
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
     IF(loc_np <= maxPart ) THEN
      DO
       IF(dyn_cell(dyn_cell(act_n_dyncell)%down_cell)%up_cell - act_n_dyncell == no_dcells - 1 &
        .AND. dyn_cell(act_n_dyncell)%down_cell /= 0) THEN
        ! we will move one level lower
        act_n_dyncell = dyn_cell(act_n_dyncell)%down_cell
       ELSE
        EXIT
       END IF
      END DO
     ! we have to create additional cells otherwise
     ELSE
     ! is there some free space left in the field dyn_cell?
      bound = 
      if(act_n_dyncell + 8 > up_bound) then
       ! define a new upper bound
       newbound = 2 * up_bound
       ALLOCATE(pom2(up_bound))
       do I = 1, up_bound
        pom(I) = local_particle(I)
       end do
       DEALLOCATE(local_particle(:))
       ALLOCATE(local_particle(newbound))
       do I = 1, up_bound
        local_particle(I) = pom(I)
       end do
       DEALLOCATE(pom(:))
       up_bound = newbound
      end if
      local_particle(np) = virtual_particle(J)
     END IF
      dyn_cell(act_n_dyncell)%up_cell = max_n_dcell + 1
      DO I = 1, no_dcells
       dyn_cell(max_n_dcell + I)%width(1) = loc_cell_width(1) / 2.E0
       dyn_cell(max_n_dcell + I)%width(2) = loc_cell_width(2) / 2.E0
       dyn_cell(max_n_dcell + I)%width(3) = loc_cell_width(3) / 2.E0
       dyn_cell(max_n_dcell + I)%down_cell = act_n_dyncell
      END DO
      ! the first cell
      dyn_cell(max_n_dcell + 1)%corner(1) = loc_corner(1)
      dyn_cell(max_n_dcell + 1)%corner(2) = loc_corner(2)
      dyn_cell(max_n_dcell + 1)%corner(3) = loc_corner(3)
      dyn_cell(max_n_dcell + 1)%down_cell = act_n_dyncell
      ! the second cell
      dyn_cell(max_n_dcell + 2)%corner(1) = loc_corner(1) + loc_cell_width(1) / 2.E0
      dyn_cell(max_n_dcell + 2)%corner(2) = loc_corner(2)
      dyn_cell(max_n_dcell + 2)%corner(3) = loc_corner(3)
      dyn_cell(max_n_dcell + 2)%down_cell = act_n_dyncell
      ! the third cell
      dyn_cell(max_n_dcell + 3)%corner(1) = loc_corner(1)
      dyn_cell(max_n_dcell + 3)%corner(2) = loc_corner(2) + loc_cell_width(2) / 2.E0
      dyn_cell(max_n_dcell + 3)%corner(3) = loc_corner(3)
      dyn_cell(max_n_dcell + 3)%down_cell = act_n_dyncell
      ! the forth cell
      dyn_cell(max_n_dcell + 4)%corner(1) = loc_corner(1) + loc_cell_width(1) / 2.E0
      dyn_cell(max_n_dcell + 4)%corner(2) = loc_corner(2) + loc_cell_width(2) / 2.E0
      dyn_cell(max_n_dcell + 4)%corner(3) = loc_corner(3)
      dyn_cell(max_n_dcell + 4)%down_cell = act_n_dyncell
      ! the fifth cell
      dyn_cell(max_n_dcell + 5)%corner(1) = loc_corner(1)
      dyn_cell(max_n_dcell + 5)%corner(2) = loc_corner(2)
      dyn_cell(max_n_dcell + 5)%corner(3) = loc_corner(3) + loc_cell_width(3) / 2.E0
      dyn_cell(max_n_dcell + 5)%down_cell = act_n_dyncell
      ! the sixth cell
      dyn_cell(max_n_dcell + 6)%corner(1) = loc_corner(1) + loc_cell_width(1) / 2.E0
      dyn_cell(max_n_dcell + 6)%corner(2) = loc_corner(2)
      dyn_cell(max_n_dcell + 6)%corner(3) = loc_corner(3) + loc_cell_width(3) / 2.E0
      dyn_cell(max_n_dcell + 6)%down_cell = act_n_dyncell
      ! the seventh cell
      dyn_cell(max_n_dcell + 7)%corner(1) = loc_corner(1)
      dyn_cell(max_n_dcell + 7)%corner(2) = loc_corner(2) + loc_cell_width(2) / 2.E0
      dyn_cell(max_n_dcell + 7)%corner(3) = loc_corner(3) + loc_cell_width(3) / 2.E0
      dyn_cell(max_n_dcell + 7)%down_cell = act_n_dyncell
      ! the eighth cell
      dyn_cell(max_n_dcell + no_dcells)%corner(1) = loc_corner(1) + loc_cell_width(1) / 2.E0
      dyn_cell(max_n_dcell + no_dcells)%corner(2) = loc_corner(2) + loc_cell_width(2) / 2.E0
      dyn_cell(max_n_dcell + no_dcells)%corner(3) = loc_corner(3) + loc_cell_width(3) / 2.E0
      dyn_cell(max_n_dcell + no_dcells)%down_cell = act_n_dyncell
      ! we have to increase the variable max_n_dcell
      max_n_dcell = max_n_dcell + no_dcells
     END IF
    END DO
   CLOSE(15)
  END SUBROUTINE create_dynamical_grid_cell
