SUBROUTINE bound_free_rates(approximation, pack_index, act_level, rate)
USE types
USE constants
IMPLICIT NONE

! input variables
INTEGER                         :: approximation
INTEGER                         :: pack_index, act_level
! ion informations
INTEGER                         :: indexe, indexi, freq
! grid informations
INTEGER                         :: current_mgi
DOUBLE PRECISION                :: el_dens, temp, gl_pop_ip1e, x
! output variables
DOUBLE PRECISION                :: rate

SELECT CASE(approximation)

!_______________________________________________________________
! hydrogenic approximation
! cross section is proportional to (f_ijk/f)^3
CASE (1)
 ! element index
 indexe = linelist(act_level)%indexe
 ! ion index
 indexi = linelist(act_level)%indexi
 ! actual model grid index
 current_mgi = get_package_model_index(pack_index)
 ! electron density
 el_dens = model_grid(current_mgi)%e_dens
 ! temperature
 temp = model_grid(current_mgi)%T
 ! number density of a ground state of ion indexi + 1, indexe
 gl_pop_ip1e = model_grid(current_mgi)%grid_comp(indexe)%grid_ion(indexi + 1)%gl_pop
 ! frequency
 freq = (element(indexe)%ion(indexi + 1)%level(0) - &
        element(indexe)%ion(indexi)%level(act_level)) / h
 ! argument of E_1(x)
 x = (h * freq) / (const_kB * temp)
 

CASE DEFAULT
 STOP 'bound_free_rates: non valid approximation was chosen'
END SELECT





END SUBROUTINE bound_free_rates
