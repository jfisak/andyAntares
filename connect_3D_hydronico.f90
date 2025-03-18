SUBROUTINE connect_3D_hydronico()

USE types
IMPLICIT NONE

INTEGER                                 :: n_bas_pcell
INTEGER                         :: cur_bpgi
INTEGER                         :: cur_neighbour, last_index
DOUBLE PRECISION, DIMENSION(const_dimofspace)  :: cur_pos

DOUBLE PRECISION, DIMENSION(const_dimofspace) :: cur_center
INTEGER, ALLOCATABLE            :: list_index(:)
INTEGER                         :: new_index
INTEGER                         :: count_vacuum, count_out, count_in, count_ok
DOUBLE PRECISION, DIMENSION(const_dimofspace)  :: cur_corner, cur_width
INTEGER                         :: up_cell, cur_lowcell
INTEGER                        :: cur_mcell
INTEGER                         :: cur_n_mgi
INTEGER                         :: cur_pgcell, cur_pgi
INTEGER                         :: n_virtpoints
DOUBLE PRECISION               :: delta!, delta2
DOUBLE PRECISION               :: dist
INTEGER                         :: ind_J
DOUBLE PRECISION                :: r_0
LOGICAL                                         :: vacuum_found
DOUBLE PRECISION, PARAMETER     :: large_number=1.d90



INTEGER                         :: down_cell
INTEGER                         :: start_vp_index, end_vp_index
INTEGER                         :: cur_vp_nearest, cur_vp_index
INTEGER                         :: act_pgcell, cur_vp, cur_mgi

! create an array with saved indexes
! write(*,*) 'connection_prop_model_grid: inputmodel = ', inputmodel, ' dyngrid = ', dyngrid
count_ok = 0
count_in = 0
count_out = 0
count_vacuum = 0

IF(dyngrid > 0) THEN
 n_bas_pcell = nx_cell * ny_cell * nz_cell
 ALLOCATE(list_index(n_bas_pcell))
 ! initial setup
 DO cur_pgcell = 1, n_bas_pcell
  list_index(cur_pgcell) = 0
 END DO
 ! calculating of indeces
 last_index = 0
 n_virtpoints = SIZE(virtual_point)
 DO cur_pgcell = 1, n_virtpoints
  new_index = virtual_point(cur_pgcell)%ind_pcell
  IF(new_index /= last_index) THEN
   list_index(new_index) = cur_pgcell
   last_index = new_index
  END IF
 END DO

 ! we calculate associated cells for the rest of propGrid cells
 DO cur_pgcell = 1, n_propgcells
  up_cell = dyn_cell(cur_pgcell)%up_cell
  IF(up_cell == 0) THEN
   cur_corner = dyn_cell(cur_pgcell)%corner
   cur_width = dyn_cell(cur_pgcell)%width
   cur_center = cur_corner + cur_width/2.0
   cur_lowcell = dyn_cell(cur_pgcell)%down_cell
   r_0 = sqrt(cur_center(ind_x)**2.0 + cur_center(ind_y)**2.0 + &
    cur_center(ind_z)**2.0)
   ! add_mg = 1 r < R_star
   ! add_mg = 2 r > R_inf
   ! add_mg = 3 r > R_star && r < R_inf, vacuum cell
   IF(r_0 < R_star) THEN
    dyn_cell(cur_pgcell)%model_index = photosphere_index
   ELSE IF(r_0 > R_inf) THEN
    dyn_cell(cur_pgcell)%model_index = outerspace_index
   ELSE 
    ! we connect a modgrid from the current basic cell
    act_pgcell = cur_pgcell
    DO
     down_cell = dyn_cell(act_pgcell)%down_cell
     IF(down_cell == 0) EXIT
     act_pgcell = down_cell
    END DO
    cur_bpgi = act_pgcell

    start_vp_index = list_index(cur_bpgi)
    

    IF(start_vp_index == 0) THEN
     dyn_cell(cur_pgcell)%model_index = vacuum_index
    ELSE
     ! we must find an end index
     cur_vp_index = start_vp_index
     DO
      cur_vp_index = cur_vp_index + 1
      if(cur_vp_index == n_virtpoints) then
       end_vp_index = cur_vp_index
       EXIT
      end if
      if(virtual_point(cur_vp_index)%ind_pcell /= cur_bpgi) then
       end_vp_index = cur_vp_index - 1
       EXIT
      end if
     END DO

     delta = large_number
     cur_vp_nearest = 0
     DO cur_vp = start_vp_index, end_vp_index
      cur_pos = virtual_point(cur_vp)%pos
      dist = sqrt((cur_center(ind_x) - cur_pos(ind_x))**2.0 +&
       (cur_center(ind_y) - cur_pos(ind_y))**2.0 +&
       (cur_center(ind_z) - cur_pos(ind_z))**2.0)
      if(dist< delta) then
       delta = dist
       cur_vp_nearest = cur_vp
      end if
     END DO
     cur_mgi = virtual_point(cur_vp_nearest)%ind_mcell
     dyn_cell(cur_pgcell)%model_index = cur_mgi
     model_grid(cur_mgi)%assoc_cells = model_grid(cur_mgi)%assoc_cells + 1
      write(*,*) 'connection_prop_model_grid: cur_pgcell = ', cur_pgcell, ' modindex = ', virtual_point(cur_vp_nearest)%ind_mcell

    END IF
   END IF
  END IF ! up_cell == 0
 END DO 
ELSE IF(dyngrid == 0) THEN
 DO cur_pgcell = 1, n_modelgrid
  cur_pos = model_grid(cur_pgcell)%vec_pos
  CALL find_dyn_cell1(cur_pos, cur_pgi)
  IF(dyn_cell(cur_pgi)%model_index == 0) THEN
   dyn_cell(cur_pgi)%model_index = cur_pgcell
   model_grid(cur_pgcell)%assoc_cells = model_grid(cur_pgcell)%assoc_cells + 1
   count_ok = count_ok + 1
  END IF
  DO ind_J = 1,6
   cur_neighbour = dyn_cell(cur_pgi)%neighbor(ind_J)
   IF(cur_neighbour > 0) then
    cur_n_mgi = dyn_cell(cur_neighbour)%model_index
    IF(cur_n_mgi == 0) THEN
     dyn_cell(cur_neighbour)%model_index = cur_pgcell
     model_grid(cur_pgcell)%assoc_cells = model_grid(cur_pgcell)%assoc_cells + 1
     count_ok = count_ok + 1
    END IF
   END IF
  END DO
 END DO
 ! 
 vacuum_found = .true.
 DO WHILE(vacuum_found)
  vacuum_found = .false.
  DO cur_pgi = 1, n_propgcells
   cur_mcell = dyn_cell(cur_pgi)%model_index
   if(cur_mcell == 0) then
    vacuum_found = .true.
   else
    DO ind_J = 1,6
     IF(ind_J == 5) cycle
     cur_neighbour = dyn_cell(cur_pgi)%neighbor(ind_J)
     if(cur_neighbour > 0) then
      cur_n_mgi = dyn_cell(cur_neighbour)%model_index
      if(cur_n_mgi == 0) then
       dyn_cell(cur_neighbour)%model_index = cur_mcell
       model_grid(cur_mcell)%assoc_cells = model_grid(cur_mcell)%assoc_cells + 1
       ! write(*,*) 'connection_prop_model_grid: cur_pgcell = ', cur_pgcell, ' cur_mcell = ', cur_mcell
      end if ! cur_n_mgi == 0
     end if ! cur_neighbour > 0
    END DO 
   end if ! cur_mcell == 0
  END DO
 END DO
 ! other cells
 DO cur_pgcell = 1, n_propgcells
  cur_corner = dyn_cell(cur_pgcell)%corner/R_star
  cur_width = dyn_cell(cur_pgcell)%width/R_star
  cur_center = cur_corner + cur_width/2.0
  r_0 = sqrt(cur_center(ind_x)**2.0 + cur_center(ind_y)**2.0 + &
   cur_center(ind_z)**2.0)
  IF(r_0 < 1.0) THEN
   dyn_cell(cur_pgcell)%model_index = photosphere_index
   count_in = count_in + 1
  ELSE IF(r_0 > R_inf/R_star) THEN
   dyn_cell(cur_pgcell)%model_index = outerspace_index
   count_out = count_out + 1
  END IF
  IF(dyn_cell(cur_pgcell)%model_index == 0) THEN
   dyn_cell(cur_pgcell)%model_index = vacuum_index
   count_vacuum = count_vacuum + 1
   ! delta2 = large_number
   ! DO J = 1, n_modelgrid
   !  mod_pos = model_grid(J)%vec_pos
   !  delta = sqrt((cur_center(ind_x)-mod_pos(ind_x))**2.0 + (cur_center(2)-mod_pos(2))**2.0 + (cur_center(3)-mod_pos(3))**2.0)
   !  if(delta < delta2) THEN
   !   delta2 = delta
   !   best_index = J
   !  end if
   !  if(mod(cur_pgcell,1000)==0) write(*,*) 'connection_prop_model_grid: working on propGrid cell ', cur_pgcell, ' from ', max_n_dcell
   ! END DO ! loop over all modGrid cells
   ! dyn_cell(cur_pgcell)%model_index = M
  END IF ! if model_index == 0
 END DO ! loop over all propGrid cells
 ! write(*,*) 'connection_prop_model_grid: in = ', count_in, ' out = ', count_out, ' vacuum = ', count_vacuum, &
 !  ' count_ok = ', count_ok
END IF ! dyncell > 0











END SUBROUTINE connect_3D_hydronico
