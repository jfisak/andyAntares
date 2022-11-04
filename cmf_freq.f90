! this sbr returns cmf frequency for the input position and the input rf frequency
SUBROUTINE cmf_freq(cur_pos, cur_freq_rf, cur_freq_cmf)

USE types
IMPLICIT NONE

DOUBLE PRECISION, DIMENSION(3)          :: cur_pos
DOUBLE PRECISION                        :: cur_freq_rf

DOUBLE PRECISION                        :: doppler

CALL doppler_factor2(cur_pos, cur_mgi, doppler)

cur_freq_cmf = cur_freq_rf * doppler






END SUBROUTINE cmf_freq

