MODULE counters

IMPLICIT NONE

! r-packets destroying
INTEGER                         :: count_des_phot = 0
INTEGER                         :: count_des_inte = 0
INTEGER                         :: count_des_esca = 0
INTEGER                         :: count_des_resd = 0
! INTEGER                         :: count_des_ipack = 0



! interaction
! r-packets counters
INTEGER                         :: count_r_line = 0
INTEGER                         :: count_r_thom = 0
INTEGER                         :: count_r_ph_k = 0
INTEGER                         :: count_r_ph_i = 0
INTEGER                         :: count_r_ff   = 0

! k-packets counters
INTEGER                         :: count_cool_ex = 0
INTEGER                         :: count_cool_ff = 0
INTEGER                         :: count_cool_io = 0
INTEGER                         :: count_cool_fb = 0

! i-packets counters
INTEGER                         :: count_i_int_down = 0
INTEGER                         :: count_i_rad_dxrs = 0
INTEGER                         :: count_i_rad_dxfl = 0
INTEGER                         :: count_i_rad_deex = 0
INTEGER                         :: count_i_int_upwa = 0
INTEGER                         :: count_i_col_deex = 0
INTEGER                         :: count_i_int_phot = 0
INTEGER                         :: count_i_int_reco = 0
INTEGER                         :: count_i_rad_reco = 0
INTEGER                         :: count_i_col_reco = 0

! d-packets counters
INTEGER                         :: count_d_change_cell = 0
INTEGER                         :: count_d_new_choice = 0
INTEGER                         :: count_d_radiative = 0
INTEGER                         :: count_d_rad_end = 0

! progrid cell statistics
INTEGER                         :: count_pg_vacuum = 0
INTEGER                         :: count_pg_mcell = 0

INTEGER                         :: count_too_dist = 0

END MODULE counters
