! this sbr should sign every model cell if a diffusion approximation
SUBROUTINE diffusion_approximation(cur_mgi)

USE types
IMPLICIT NONE

INTEGER                                 :: cur_mgi

DOUBLE PRECISION                        :: electron_density, chi_cont

DOUBLE PRECISION, PARAMETER             :: chi_min = 1e-16


electron_density = model_grid(cur_mgi)%e_dens
chi_cont = sigma_e * electron_density
write(*,*) 'diffusion_approximation: cur_mgi = ', cur_mgi, ' chi_cont = ', chi_cont

IF(chi_cont > chi_min) THEN
 model_grid(cur_mgi)%is_difapp = .true.
END IF








END SUBROUTINE diffusion_approximation
