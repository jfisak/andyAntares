! this sbr returns cmf frequency for the input position and the input rf frequency
SUBROUTINE cmf_freq(pack_index, cur_pos, cur_freq_rf, cur_freq_cmf)

USE constants
USE types
IMPLICIT NONE

INTEGER                                 :: pack_index
DOUBLE PRECISION, DIMENSION(3)          :: cur_pos
DOUBLE PRECISION                        :: cur_freq_rf, cur_freq_cmf

DOUBLE PRECISION                        :: doppler

INTEGER                                 :: cur_mgi

CALL doppler_factor2(pack_index, cur_pos, cur_mgi, doppler)

cur_freq_cmf = cur_freq_rf * doppler






END SUBROUTINE cmf_freq

