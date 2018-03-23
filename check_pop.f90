!_______________________ check_pop _________________________________
! this sbr checks if the population of the ionic levels corresponds
! to the total population of the given ion
SUBROUTINE check_pop()

USE types
IMPLICIT NONE

INTEGER                                 :: cur_mgi, indexe, indexi, indexl
INTEGER                                 :: n_ions, n_levels
DOUBLE PRECISION                        :: tot_pop, sum_pop, act_pop
DOUBLE PRECISION                        :: ratio

DO cur_mgi = 1, n_modelgrid
 IF (model_grid(cur_mgi)%assoc_cells == 0) CYCLE
 DO indexe = 1, n_elements
  n_ions = SIZE(elements(indexe)%ions)
  DO indexi = 1, n_ions
   tot_pop = model_grid(cur_mgi)%grid_comp(indexe)%grid_ion(indexi)%tot_pop
   ! write(*,*) 'check_pop: tot_pop = ', tot_pop
   n_levels = SIZE(elements(indexe)%ions(indexi)%levels)
   sum_pop = 0.D0
   DO indexl = 1, n_levels
    CALL populations(indexe, indexi, indexl, cur_mgi, act_pop)
    sum_pop = sum_pop + act_pop
   END DO
   ! write(*,*) 'check_pop: 2 tot_pop = ', tot_pop
   ratio = sum_pop / tot_pop - 1.D0
   IF(ABS(ratio) > 0.05) THEN
    write(*,*) 'check_pop: WARNING: el populations do not fit the total population'
    write(*,*) ' element index = ', indexe, ' ion = ', indexi
    write(*,*) ' sum_pop = ', sum_pop, ' tot_pop = ', tot_pop
   END IF
   IF(ratio > 1.D0) THEN
    write(*,*) 'check_pop: WARNING: el subpopulations are larger than the total population'
    write(*,*) ' element index = ', indexe, ' ion = ', indexi
    write(*,*) ' sum_pop = ', sum_pop, ' tot_pop = ', tot_pop
   END IF
  END DO
 END DO
END DO

END SUBROUTINE check_pop
