!___________________________________________________________________
!_______________________ check_pop _________________________________
!___________________________________________________________________
! this sbr checks if the population of the ionic levels corresponds
! to the total population of the given ion
SUBROUTINE check_pop()

USE types
USE constants
IMPLICIT NONE

INTEGER                                 :: cur_mgi, indexe, indexi, indexl
DOUBLE PRECISION                        :: cur_radius
INTEGER                                 :: n_ions, n_levels, at_index
DOUBLE PRECISION                        :: tot_pop_ion, sum_pop, act_pop
DOUBLE PRECISION                        :: ratio
DOUBLE PRECISION                        :: tot_pop_ele
DOUBLE PRECISION                        :: frac, tot_frac
DOUBLE PRECISION                        :: el_nd, temp
CHARACTER(LEN=60)                       :: filename

filename = trim(outputfolder)//'/ionFracs.dat'
OPEN(59,FILE=filename)

DO cur_mgi = 1, n_modelgrid
 IF (model_grid(cur_mgi)%assoc_cells == 0) CYCLE
 DO indexe = 1, n_elements
  n_ions = SIZE(elements(indexe)%ions)
  cur_radius = model_grid(cur_mgi)%rwind
  tot_pop_ele = 0.D0
  tot_pop_ion = 0.D0
  tot_frac = 0.D0
  DO indexi = 1, n_ions
   temp = model_grid(cur_mgi)%t
   el_nd = model_grid(cur_mgi)%e_dens
   CALL ionization_fraction(indexe, indexi, temp, el_nd, frac)
   at_index = elements(indexe)%atom_number
   write(59,*) cur_mgi, at_index, indexi, frac, cur_radius
   tot_frac = tot_frac + frac
   tot_pop_ion = model_grid(cur_mgi)%grid_comp(indexe)%grid_ion(indexi)%tot_pop
   tot_pop_ele = tot_pop_ele + tot_pop_ion
   ! write(*,*) 'check_pop: tot_pop_ion = ', tot_pop_ion
   n_levels = SIZE(elements(indexe)%ions(indexi)%levels)
   sum_pop = 0.D0
   DO indexl = 1, n_levels
    CALL populations(indexe, indexi, indexl, cur_mgi, act_pop)
    sum_pop = sum_pop + act_pop
   END DO
   ! write(*,*) 'check_pop: 2 tot_pop_ion= ', tot_pop_ion
   ratio = sum_pop / tot_pop_ion - 1.D0
   IF(ABS(ratio) > 0.05) THEN
    write(99,*) 'check_pop: WARNING: el populations do not fit the total population'
    write(99,*) ' element index = ', indexe, ' ion = ', indexi
    write(99,*) ' sum_pop = ', sum_pop, ' tot_pop_ion = ', tot_pop_ion
   END IF
   IF(ratio > 1.D0) THEN
    write(99,*) 'check_pop: WARNING: el subpopulations are larger than the total population'
    write(99,*) ' element index = ', indexe, ' ion = ', indexi
    write(99,*) ' sum_pop = ', sum_pop, ' tot_pop_ion = ', tot_pop_ion
   END IF
  END DO
  ! tot_pop_elea = model_grid(cur_mgi)%grid_comp(indexe)%abund * model_grid(cur_mgi)%rho / elements(indexe)%atom_mass
  ! ratio_el = tot_pop_ele / tot_pop_elea
  ! write(*,*) 'check_pop: indexe = ', indexe, ' ratio_el = ', ratio_el
  write(59,*) cur_mgi, at_index, ' TOT ', tot_frac, cur_radius
 END DO
END DO
CLOSE(59)

END SUBROUTINE check_pop
