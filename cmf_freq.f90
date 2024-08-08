! this sbr returns cmf frequency for the input position and the input rf frequency
SUBROUTINE cmf_freq(pack_index,cur_freq_rf, cur_freq_cmf)

USE constants
USE types
USE dummypacket

IMPLICIT NONE

INTEGER                                 :: pack_index
DOUBLE PRECISION                        :: cur_freq_rf, cur_freq_cmf

DOUBLE PRECISION                        :: doppler

INTEGER                                 :: cur_mgi

CALL doppler_factor(pack_index, doppler)

cur_freq_cmf = cur_freq_rf * doppler

! write(*,*) 'cmf_freq: doppler = ', doppler






END SUBROUTINE cmf_freq

