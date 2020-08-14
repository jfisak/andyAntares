! this subroutine computes approximative value of the gamma
! function for the given x
! \Gamma(x) = max(\tilde{g}, 0.26 * exp(x) * E_1(x))
SUBROUTINE gamma_function(x, transition, gf)
USE types
IMPLICIT NONE

! input variables
DOUBLE PRECISION                        :: x
INTEGER                                 :: transition
DOUBLE PRECISION, PARAMETER             :: gf_const = 2.76D0
! exponential integal function
DOUBLE PRECISION                        :: eif
! computed values
DOUBLE PRECISION                        :: value_1, value_2
! output variables
DOUBLE PRECISION                        :: gf
INTEGER                                 :: indexe, indexi
CHARACTER(LEN=15)                       :: el_conf_lower, el_conf_upper
INTEGER                                 :: conf_len_lower, conf_len_upper
CHARACTER(LEN=2)                        :: sconf_l, sconf_u

! THE FIRST VALUE
! firstly we will calculate \tilde{g}
! \tilde{g} is equal to 0.7 for transitions <n, l> -> <n, l'>
!                       0.2 for transitions <n, l> -> <n', l'>
! now we have to extract which type of transition occurs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! now we choose only one value
! WORKS ONLY FOR THE OPACITY PROJECT DATA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
indexe = linelist(transition)%indexe
indexi = linelist(transition)%indexi
el_conf_lower = elements(indexe)%ions(indexi)%levels(linelist(transition)%lower)%elconf
el_conf_upper = elements(indexe)%ions(indexi)%levels(linelist(transition)%upper)%elconf
! lengths of these arrays without spaces
conf_len_lower = LEN_TRIM(el_conf_lower)
conf_len_upper = LEN_TRIM(el_conf_upper)

! write(*,*) 'gamma_function: testing configurations'
! write(*,*) 'elclower: conf_len_lower = ', conf_len_lower
! write(*,*) 'el_conf_lower = ', el_conf_lower, ' el_conf_upper = ', el_conf_upper
! write(*,*) 'elclower: ', el_conf_lower(conf_len_lower - 2:conf_len_lower - 1)
! write(*,*) 'elclower: ', el_conf_upper(conf_len_upper - 2:conf_len_upper - 1)

IF(orbitals_nl .EQV. .TRUE.) THEN
 IF(simpleTrans) THEN
  sconf_l = el_conf_lower(1:conf_len_lower - 1)
  sconf_u = el_conf_upper(1:conf_len_upper - 1)
 ELSE
  sconf_l = el_conf_lower(conf_len_lower - 2:)
  sconf_u = el_conf_upper(conf_len_upper - 2:conf_len_upper - 1)
 END IF
 
 IF(sconf_l == sconf_u) THEN
  value_1 = 7.D-1
 ELSE
  value_1 = 2.D-1
 END IF
ELSE
 value_1 = 7.D-1
END IF
!write(*,*) 'gamma_function: value_1 = ', value_1
!STOP 'gamma_function: only for testing'

!print*, 'el_conf_lower = ', el_conf_lower, ' el_conf_upper = ', el_conf_upper
!STOP
! THE SECOND VALUE
! the main thing in this part is to calculate an exponential integral function
CALL exp_int_func(1, x, eif)

value_2 = gf_const * exp(x) * eif

gf = MAX(value_1, value_2)

END SUBROUTINE gamma_function
