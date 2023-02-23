! this subroutine will find dynamical cell from the given position
  SUBROUTINE find_dyn_cell(pos,actual_cell)

   USE types
USE constants
   IMPLICIT NONE

   INTEGER                              :: actual_cell
   ! position
   DOUBLE PRECISION, DIMENSION(3)       :: pos, width
   LOGICAL                              :: xpl, ypl, zpl
   ! local variables
   DOUBLE PRECISION, DIMENSION(3)       :: loc_wid, loc_cor
   INTEGER                              :: ind_cell_numb
   INTEGER                              :: ind_x, ind_y, ind_z
   INTEGER                              :: act_ndgrid
   width = dyn_cell(1)%width
   ! firstly we have to know the basic cell of the photon
   ind_x = FLOOR(pos(1)/width(1) + DBLE(nx_cell)/2) + 1
   ind_y = FLOOR(pos(2)/width(2) + DBLE(ny_cell)/2) + 1
   ind_z = FLOOR(pos(3)/width(3) + DBLE(nz_cell)/2) + 1
   IF(ind_x <= 0 .OR. ind_x > nx_cell ) STOP 'out of basic grid'
   IF(ind_y <= 0 .OR. ind_x > ny_cell ) STOP 'out of basic grid'
   IF(ind_z <= 0 .OR. ind_x > nz_cell ) STOP 'out of basic grid'
   ind_cell_numb = (ind_x - 1) * ny_cell * nz_cell + (ind_y - 1) * nz_cell + ind_z
   !print*, 'inx_x = ', ind_x, 'ind_y = ', ind_y, 'ind_z = ', ind_z, 'ind_cell_numb = ', ind_cell_numb

   act_ndgrid = ind_cell_numb
   ! now we are looking for the corresponding dynamic grid in the basic cell
   DO
    !print*, 'actual cell: ', act_ndgrid
    loc_cor = dyn_cell(act_ndgrid)%corner
    loc_wid = dyn_cell(act_ndgrid)%width
    ! now we will test the position of the particle in every eighth of cell
    ! => we will find the corresponding subcell
    IF(dyn_cell(act_ndgrid)%up_cell /= 0) THEN
     if(pos(1) >= loc_cor(1) .AND. pos(1) < loc_cor(1) + loc_wid(1)/2.E0 ) then
      !print*, 'it is false...'
      xpl = .FALSE.
     else if(pos(1) >= loc_cor(1) + loc_wid(1)/2.E0 .AND. pos(1) < loc_cor(1) + loc_wid(1)) then
      !print*, 'it is true...'
      xpl = .TRUE.
     else
      STOP 'it is not in the cell...'

     end if
     if(pos(2) >= loc_cor(2) .AND. pos(2) < loc_cor(2) + loc_wid(2)/2.E0 ) then
      ypl = .FALSE.
     else if(pos(2) >= loc_cor(2) + loc_wid(2)/2.E0 .AND. pos(2) < loc_cor(2) + loc_wid(2)) then
      ypl = .TRUE.
     else
      STOP 'it is not in the cell...'
     end if
     if(pos(3) >= loc_cor(3) .AND. pos(3) < loc_cor(3) + loc_wid(3)/2.E0 ) then
      zpl = .FALSE.
     else if(pos(3) >= loc_cor(3) + loc_wid(3)/2.E0 .AND. pos(3) < loc_cor(3) + loc_wid(3)) then
      zpl = .TRUE.
     else
      STOP 'it is not in the cell...'
     end if
     !print*, 'xpl = ', xpl, 'ypl = ', ypl, 'zpl = ', zpl
     IF(xpl .EQV. .FALSE.) THEN
      IF(ypl .EQV. .FALSE.) THEN
       IF(zpl .EQV. .FALSE.) THEN
        act_ndgrid = dyn_cell(act_ndgrid)%up_cell
        ! zpl .EQV. .TRUE.
       ELSE
        act_ndgrid = dyn_cell(act_ndgrid)%up_cell + 4
       END IF
        ! ypl .EQV. .TRUE.
      ELSE
       IF(zpl .EQV. .FALSE.) THEN
        act_ndgrid = dyn_cell(act_ndgrid)%up_cell + 2
        ! zpl .EQV. .TRUE.
       ELSE
        act_ndgrid = dyn_cell(act_ndgrid)%up_cell + 6
       END IF
      END IF
        ! xpl .EQV. .TRUE.
     ELSE
      IF(ypl .EQV. .FALSE.) THEN
       IF(zpl .EQV. .FALSE.) THEN
        act_ndgrid = dyn_cell(act_ndgrid)%up_cell + 1
        ! zpl .EQV. .TRUE.
       ELSE
        act_ndgrid = dyn_cell(act_ndgrid)%up_cell + 5
       END IF
        ! ypl .EQV. .TRUE.
      ELSE
       IF(zpl .EQV. .FALSE.) THEN
        act_ndgrid = dyn_cell(act_ndgrid)%up_cell + 3
        ! zpl .EQV. .TRUE.
       ELSE
        act_ndgrid = dyn_cell(act_ndgrid)%up_cell + 7
       END IF
      END IF
     END IF
    ELSE
     actual_cell = act_ndgrid
     EXIT
    END IF
   END DO

  END SUBROUTINE find_dyn_cell
