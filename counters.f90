MODULE counters

IMPLICIT NONE

! r-packets counters
INTEGER                         :: count_r_line = 0
INTEGER                         :: count_r_thom = 0
INTEGER                         :: count_r_phot = 0
INTEGER                         :: count_r_ff   = 0

! k-packets counters
INTEGER                         :: count_cool_ex = 0
INTEGER                         :: count_cool_ff = 0
INTEGER                         :: count_cool_io = 0
INTEGER                         :: count_cool_fb = 0

END MODULE counters
