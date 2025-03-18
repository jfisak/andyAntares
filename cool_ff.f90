! calculating free-free cooling rates
! 
! INPUT: cur_mgi(INT): current modGrid index
! OUTPUT: Zcool(DBLE): total free-free collision index
!
SUBROUTINE cool_ff(cur_mgi, Zcool)
USE types
USE constants
USE rates_k
IMPLICIT NONE
! constant
DOUBLE PRECISION, PARAMETER                             :: C0 = 1.426D-27 ! Osterbrock 1974
! element index
INTEGER                                                 :: indexe
! ion index and other informations
INTEGER                                                 :: indexi, n_ions
! model grid information
DOUBLE PRECISION                                        :: ion_charge
INTEGER                                                 :: cur_mgi
DOUBLE PRECISION                                        :: cur_temp, e_dens
INTEGER                                                 :: act_cooling
DOUBLE PRECISION                                        :: act_pop, actVal
! output variables
DOUBLE PRECISION                                        :: Zcool
                                                       
cur_temp = model_grid(cur_mgi)%t
e_dens = model_grid(cur_mgi)%e_dens
!write(*,*) 'cool_ff: pack_index = ', pack_index, ' cur_mgi = ', cur_mgi, 'e_dens = ', e_dens

act_cooling = 0
Zcool = 0.D0
DO indexe = 1, n_elements
 n_ions = SIZE(elements(indexe)%ions)
 DO indexi = 2, n_ions
  act_cooling = act_cooling + 1
  ion_charge = DBLE(indexi - 1)
  act_pop = model_grid(cur_mgi)%grid_comp(indexe)%grid_ion(indexi)%tot_pop
  actVal = C0 * ion_charge**2.0 * cur_temp**(1.0/2.0) &
                          * act_pop * e_dens
  Zcool = Zcool + actVal
  ! write(*,*) 'cool_ff: ion_charge = ', ion_charge, ' cur_temp = ', cur_temp, &
  !  ' act_pop = ', act_pop, ' e_dens = ', e_dens, ' cooling rate = ', actVal
 END DO
END DO



END SUBROUTINE cool_ff
