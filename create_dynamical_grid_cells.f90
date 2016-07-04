! this subroutine will create a dynamical grid cell for the given basic cell
  SUBROUTINE create_dynamical_grid_cell()

   USE types

   IMPLICIT NONE
   
   INTEGER                       :: I,J
   INTEGER                       :: Npart, np
   DOUBLE PRECISION, DIMENSION(3):: corner, cell_width_2
   ! we have to go through every basic propagation grid cell
   DO I = 1, Ngrid
    corner = cell(I)%corner
    cell_width_2(1) = 2.D0 * xmax / nx_cell
    cell_width_2(2) = 2.D0 * ymax / nx_cell
    cell_width_2(3) = 2.D0 * zmax / nx_cell
    Npart = SIZE(virtual_particle(:))
    ! at first we have to know, how many particles are
    ! in the given cell
     np = 0
    DO J = 1, Npart
     IF((virtual_particle(J)%pos(1) .GE. corner(1)) .AND. &
        (virtual_particle(J)%pos(1) < (corner(1) + cell_width_2(1))) .AND. &
        (virtual_particle(J)%pos(2) .GE. corner(2)) .AND. &
        (virtual_particle(J)%pos(2) < (corner(2) + cell_width_2(2))) .AND. &
        (virtual_particle(J)%pos(3) .GE. corner(3)) .AND. &
        (virtual_particle(J)%pos(3) < (corner(3) + cell_width_2(3)))) THEN
      np = np + 1
     END IF
    END DO
   END DO
  END SUBROUTINE create_dynamical_grid_cell
