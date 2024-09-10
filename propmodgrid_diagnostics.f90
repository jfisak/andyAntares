! this sbr tests the created propGrid if it covers the modGrid properly
!
! input: none
! output: none (yet)
!
!
SUBROUTINE propmodgrid_diagnostics()

USE types

IMPLICIT NONE

INTEGER                                         :: cur_mgi, n_associated_pgi

INTEGER                                         :: n_alone, n_assoc
INTEGER                                         :: n_plusone = -1, n_plustwo = -1, n_plusthree = -1

REAL                                            :: rel_coverage
REAL, PARAMETER                                 :: min_rel_coverage = 0.9

INTEGER                                         :: gridcell, cur_pgi
DOUBLE PRECISION                                :: loc_volume, tot_model_volume, tot_nonmodel_volume
INTEGER, ALLOCATABLE                            :: list_points(:)
INTEGER                                         :: np_wewant



! basic tests of the propGrid vers the modGrid
! test for the created
IF(xmax < R_inf .or. ymax < R_inf .or. zmax < R_inf) THEN
 write(*,*) 'propmodgrid_diagnostics: xmax, ymax, or zmax do not fit the modGrid'
 write(*,*) 'propmodgrid_diagnostics: xmax = ', xmax, ' ymax = ', ymax, ' zmax = ', zmax
 STOP 'propmodgrid_diagnostics'
END IF

write(*,*) 'main: xmax = ', xmax/R_inf, ' ymax = ', ymax/R_inf, ' zmax = ', zmax/R_inf

n_alone = 0
n_assoc = 0

DO cur_mgi = 1, n_modelgrid
 n_associated_pgi = model_grid(cur_mgi)%assoc_cells
 ! relative covering of the modGrid by the propGrid
 IF(n_associated_pgi == 0) THEN
  n_alone = n_alone + 1
 ELSE IF(n_associated_pgi > 0) THEN
  n_assoc = n_assoc + 1
 END IF

END DO

rel_coverage = real(n_assoc) / real(n_modelgrid)

! vacuum, photospheric and outer propGrid cells
n_plusone = model_grid(n_modelgrid + 1)%assoc_cells
IF(add_mg == 2) THEN
 n_plustwo = model_grid(n_modelgrid + add_mg)%assoc_cells
END IF
IF(add_mg == 3) THEN
 n_plustwo = model_grid(n_modelgrid + 2)%assoc_cells
 n_plusthree = model_grid(n_modelgrid + add_mg)%assoc_cells
END IF

write(*,*) 'propmodgrid_diagnostics: n+1 = ', n_plusone, ' n+2 = ', n_plustwo, ' n+3 = ', n_plusthree

IF(rel_coverage < min_rel_coverage) THEN
 write(*,*) '**********************************************************'
 write(*,*) 'propmodgrid_diagnostics: warning'
 write(*,*) 'propmodgrid_diagnostics: the relative coverage is too low!'
 write(*,*) 'propmodgrid_diagnostics: rel_coverage = ', rel_coverage * 100.0, ' %'
 write(*,*) '**********************************************************'
 write(99,*) '**********************************************************'
 write(99,*) 'propmodgrid_diagnostics: warning'
 write(99,*) 'propmodgrid_diagnostics: the relative coverage is too low!'
 write(99,*) 'propmodgrid_diagnostics: rel_coverage = ', rel_coverage * 100.0, ' %'
 write(99,*) '**********************************************************'
END IF

tot_model_volume = 0.D0
tot_nonmodel_volume = 0.D0
! volume of propGrid cells adjoncet to the modGrid cells
DO cur_pgi = 1, n_propgcells
 IF(dyn_cell(cur_pgi)%up_cell == 0) THEN
  gridcell = dyn_cell(cur_pgi)%model_index
  loc_volume = dyn_cell(cur_pgi)%width(1) * dyn_cell(cur_pgi)%width(2) * dyn_cell(cur_pgi)%width(3)
  ! write(*,*) 'propmodgrid_diagnostics: loc_volume = ', loc_volume
  IF(gridcell <= n_modelgrid) THEN
   tot_model_volume = tot_model_volume + loc_volume
  ELSE IF(gridcell > n_modelgrid) THEN
   tot_nonmodel_volume = tot_nonmodel_volume + loc_volume
  END IF
  model_grid(gridcell)%volume = model_grid(gridcell)%volume + loc_volume
  ! counters
  ! IF(dyn_cell(cur_pgi)%model_index /= 0) THEN
  !  count_pg_mcell = count_pg_mcell + 1
  ! ELSE IF(dyn_cell(cur_pgi)%model_index == n_modelgrid + 3) THEN
  !  count_pg_vacuum = count_pg_vacuum + 1
  ! END IF
 END IF
END DO

! obtaining a number of points we need to calculate a volume for each modGrid cell
IF(model_type == 1) THEN
 np_wewant = 2
ELSE IF(model_type == 2) THEN
 np_wewant = 4
ELSE IF(model_type == 3) THEN
 np_wewant = 6
END IF
ALLOCATE(list_points(np_wewant))

DO cur_mgi = 1, n_modelgrid
 ! find the closest neighbors
 IF(model_type == 1) THEN
 ELSE IF(model_type == 2) THEN
  CALL n_closest_points_2D(cur_mgi, np_wewant, list_points)
 ELSE IF(model_type == 3) THEN
 END IF
END DO

! DO cur_mgi = 1, n_modelgrid
!  write(*,*) 'propmodgrid_diagnostics: cur_mgi = ', cur_mgi, ' vol = ', model_grid(cur_mgi)%volume
! END DO

write(*,*) 'propmodgrid_diagnostics: tot_model_volume = ', tot_model_volume
write(*,*) 'propmodgrid_diagnostics: tot_nonmodel_volume = ', tot_nonmodel_volume

CALL save_output(12)



write(*,*) 'propmodgrid_diagnostics: n_assoc = ', n_assoc, ' n_alone = ', n_alone
STOP 'propmodgrid_diagnostics: testing'








END SUBROUTINE propmodgrid_diagnostics
