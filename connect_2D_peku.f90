SUBROUTINE connect_2D_peku()

USE types
IMPLICIT NONE

INTEGER                                 :: N_vgrid_r, N_vgrid_t, N_vgrid_cells
DOUBLE PRECISION                        :: w_vgrid_r, w_vgrid_t
INTEGER                                 :: cur_point, n_rt_A, n_rt_B
INTEGER                                 :: cur_n_r_A, cur_n_t_A, cur_n_r_B, cur_n_t_B
INTEGER                                 :: cur_ind_A, cur_ind_B, cur_vpg_cell

DOUBLE PRECISION                        :: cur_r, cur_t
DOUBLE PRECISION                        :: rmax, rmin, tmax, tmin

INTEGER, ALLOCATABLE                    :: n_points_A(:), n_points_B(:)
INTEGER, ALLOCATABLE                    :: dummy_counter_A(:), dummy_counter_B(:)
INTEGER, ALLOCATABLE                    :: indices_A(:), indices_B(:)


! virtual grid definition
! division of virGrid A and B
N_vgrid_r = 10
N_vgrid_t = 10
N_vgrid_cells = N_vgrid_r * N_vgrid_t
ALLOCATE(n_points_A(N_vgrid_cells), n_points_B(N_vgrid_cells))
ALLOCATE(indices_A(N_vgrid_cells), indices_B(N_vgrid_cells))
n_points_A(:) = 0
n_points_B(:) = 0

rmax = MAXVAL(model_grid(:)%rwind)+1e1
rmin = MINVAL(model_grid(:)%rwind)-1e1
tmax = MAXVAL(model_grid(:)%angle)+1e-2
tmin = MINVAL(model_grid(:)%angle)-1e-2

w_vgrid_r = abs(rmax - rmin)/N_vgrid_r
w_vgrid_t = abs(tmax - tmin)/N_vgrid_t



! calculation of the virGrid index
DO cur_point = 1, n_modelgrid
 cur_r = model_grid(cur_point)%rwind
 cur_t = model_grid(cur_point)%angle
 cur_n_r_A = floor(cur_r/w_vgrid_r) + 1
 cur_n_t_A = floor(cur_t/w_vgrid_t) + 1
 cur_n_r_B = floor(cur_r/w_vgrid_r + 1/2) + 1
 cur_n_t_B = floor(cur_t/w_vgrid_t + 1/2) + 1

 n_rt_A = cur_n_r_A + N_vgrid_r * (cur_n_t_A - 1)
 ! write(*,*) 'connect_2D_peku: cur_n_t_A = ', cur_n_t_A, ' cur_n_r_A = ', cur_n_r_A, ' n_rt_A = ', n_rt_A
 n_points_A(n_rt_A) = n_points_A(n_rt_A) + 1
 n_rt_B = cur_n_r_B + N_vgrid_r * (cur_n_t_B - 1)
 n_points_B(n_rt_B) = n_points_B(n_rt_B) + 1
END DO

! index array
cur_ind_A = 0
cur_ind_B = 0

! create arrays with indeces pointing to an ordered list of modCell grids indeces
DO cur_vpg_cell = 1, N_vgrid_cells
 indices_A(cur_vpg_cell) = cur_ind_A + 1
 indices_B(cur_vpg_cell) = cur_ind_B + 1

 cur_ind_A = cur_ind_A + n_points_A(cur_vpg_cell)
 cur_ind_B = cur_ind_B + n_points_B(cur_vpg_cell)

END DO

! 
dummy_counter_A = n_points_A
dummy_counter_B = n_points_B





STOP 'connect_2D_peku: testing'




! loop over every propGrid cell to calculate associated modGrid cells



END SUBROUTINE connect_2D_peku
