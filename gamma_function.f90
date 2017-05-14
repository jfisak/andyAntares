! this subroutine computes approximative value of the gamma
! function for the given x
! \Gamma(x) = max(\tilde{g}, 0.26 * exp(x) * E_1(x))
SUBROUTINE gamma_function(x, transition, gf)
USE types
IMPLICIT NONE

! input variables
DOUBLE PRECISION                        :: x, transition
DOUBLE PRECISION, PARAMETER             :: gf_const = 2.76D0
! exponential integal function
DOUBLE PRECISION                        :: eif
! computed values
DOUBLE PRECISION                        :: value_1, value_2
! output variables
DOUBLE PRECISION                        :: gf

! THE FIRST VALUE
! firstly we will calculate \tilde{g}
! \tilde{g} is equal to 0.7 for transitions <n, l> -> <n, l'>
!                       0.2 for transitions <n, l> -> <n', l'>
! configurations of initial and final state
!iconf = linelist(transition)%iconf
!jconf = linelist(transition)%jconf
! now we have to extract which type of transition occurs
!conf_leni = LEN_TRIM(iconf)
!conf_lenj = LEN_TRIM(iconf)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! now we choose only one value
! IT MUST BE DONE LATER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
value_1 = 0.7

! THE SECOND VALUE
! the main thing in this part is to calculate an exponential integral function
CALL exp_int_func(1, x, eif)

value_2 = gf_const * exp(x) * eif

gf = MAX(value_1, value_2)

END SUBROUTINE gamma_function
