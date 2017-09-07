! calculating free-free cooling rates
SUBROUTINE cool_ff(pack_index, Zcool)
USE types
USE rates
IMPLICIT NONE
! input variables
INTEGER                                                 :: pack_index
! constant
DOUBLE PRECISION, PARAMETER                             :: C0 = 1.426D-27 ! Osterbrock 1974
! element index
INTEGER                                                 :: indexe
! ion index and other informations
INTEGER                                                 :: indexi, ion_charge, n_ions
! model grid information
INTEGER                                                 :: cur_mgi, get_package_model_index
DOUBLE PRECISION                                        :: cur_temp, e_dens
INTEGER                                                 :: act_pop, act_cooling
INTEGER                                                 :: n_coll
! output variables
DOUBLE PRECISION                                        :: Zcool
                                                       
! number of possible rates
n_coll = 0
DO indexe = 1, n_elements
 n_ions = SIZE(elements(indexe)%ions)
 DO indexi = 1, n_ions
  n_coll = n_coll + 1
 END DO
END DO
!write(*,*) 'cool_ff: n_coll = ', n_coll
! we expect that number of included ions does not change in the stellar wind
! so we do not reallocate existing array
IF(.NOT. ASSOCIATED(Lcool_ff)) THEN
 ALLOCATE(Lcool_ff(n_coll))
END IF
cur_mgi = get_package_model_index(pack_index)
cur_temp = model_grid(cur_mgi)%t
e_dens = model_grid(cur_mgi)%e_dens
!write(*,*) 'cool_ff: pack_index = ', pack_index, ' cur_mgi = ', cur_mgi, 'e_dens = ', e_dens

act_cooling = 0
Zcool = 0.D0
DO indexe = 1, n_elements
 n_ions = SIZE(elements(indexe)%ions)
 DO indexi = 1, n_ions - 1
  act_cooling = act_cooling + 1
  ion_charge = indexi - 1
  !CALL populations(I, J, 1, cur_mgi, act_pop)
  act_pop = model_grid(cur_mgi)%grid_comp(indexe)%grid_ion(indexi)%tot_pop
  Lcool_ff(act_cooling) = C0 * ion_charge * cur_temp**(1.0/2.0) &
                          * act_pop * e_dens
  Zcool = Zcool + Lcool_ff(act_cooling)
!  write(*,*) 'cool_ff: ion_charge = ', ion_charge, ' cur_temp = ', cur_temp, &
!   ' act_pop = ', act_pop, ' e_dens = ', e_dens, ' cooling rate = ', Lcool_ff(act_cooling)
 END DO
END DO



END SUBROUTINE cool_ff
