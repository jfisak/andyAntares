MODULE counters

IMPLICIT NONE

! r-packets destroying
INTEGER                         :: count_des_phot = 0
INTEGER                         :: count_des_inte = 0
INTEGER                         :: count_des_esca = 0



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

END MODULE counters
