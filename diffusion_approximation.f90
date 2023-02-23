! this sbr should sign every model cell if a diffusion approximation
SUBROUTINE diffusion_approximation(cur_mgi)

USE types
USE constants
IMPLICIT NONE

INTEGER                                 :: cur_mgi

DOUBLE PRECISION                        :: electron_density, chi_cont

LOGICAL                                 :: is_diff

DOUBLE PRECISION                        :: cur_lambda

DOUBLE PRECISION, PARAMETER             :: chi_min = 5e-16, lambda_min = 0.30


! write(*,*) 'diffusion_approximation: start'
! hydronico model, diffaaprox is decidead based on the parameter lambda
IF(model_type == 3 .and. inputmodel == 1) THEN
 cur_lambda = model_grid(cur_mgi)%diff_param
 IF(cur_lambda > lambda_min) THEN
  is_diff = .true.
  model_grid(cur_mgi)%is_difapp = is_diff
 ELSE
  is_diff = .false.
  model_grid(cur_mgi)%is_difapp = is_diff
 END IF
 ! write(*,*) 'diffusion_approximation: cur_mgi = ', cur_mgi, ' cur_lambda = ', cur_lambda, ' is_diff? = ', is_diff
ELSE
 electron_density = model_grid(cur_mgi)%e_dens
 chi_cont = sigma_e * electron_density
 
 IF(chi_cont > chi_min) THEN
  is_diff = .true.
  model_grid(cur_mgi)%is_difapp = is_diff
 ELSE 
  is_diff = .false.
  model_grid(cur_mgi)%is_difapp = is_diff
 END IF
 ! write(*,*) 'diffusion_approximation: cur_mgi = ', cur_mgi, ' chi_cont = ', chi_cont, ' is_diff? ', is_diff
 
END IF







END SUBROUTINE diffusion_approximation
