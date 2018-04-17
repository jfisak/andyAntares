SUBROUTINE k_freq_ff(pack_index, ran_freq)
USE types
IMPLICIT NONE

! input
INTEGER                         :: pack_index, indexe, indexi
! output
DOUBLE PRECISION                :: ran_freq
! random number
DOUBLE PRECISION                :: ran_z
DOUBLE PRECISION                       :: ran2
! properties of model grid
INTEGER                         :: cur_mgi
INTEGER                         :: get_package_model_index
DOUBLE PRECISION                :: e_dens, temp
DOUBLE PRECISION, PARAMETER     :: ff_const = 3.6955e8
! generate a random number
ran_z = ran2(idum)
! calculate variables of free-free emission coefficient
cur_mgi = get_package_model_index(pack_index)
e_dens = model_grid(cur_mgi)%e_dens
!Nij = model_grid(cur_mgi)%grid_comp(indexe)%grid_ion(indexi)%tot_pop
temp = model_grid(cur_mgi)%t

ran_freq = - BOLK * temp / h * log(ran_z)


END SUBROUTINE k_freq_ff
