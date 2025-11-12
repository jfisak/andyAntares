! this subroutine calculates exponential integral function for the given x
!
! INPUT: approx(INt): selected approximation
!        var_x(DBL): value x
! OUTPUT: eif(DBL): exponential integral function result
!
SUBROUTINE exp_int_func(approx, var_x, eif)
USE types
USE constants
IMPLICIT NONE
! input variable
INTEGER                                :: approx
DOUBLE PRECISION                        :: var_x
! output variable
DOUBLE PRECISION                        :: eif
! indexes
INTEGER                                 :: ind_I
! number of intervals
INTEGER                                 :: nintervals
! parameters describing the interval division
DOUBLE PRECISION                        :: al, bl
! 
DOUBLE PRECISION                        :: var_y, x1, x2
DOUBLE PRECISION                        :: loc_sum, summ
! approximation constants
DOUBLE PRECISION, PARAMETER             :: ap_1 = -0.5772156649
! calculated values of expintfunc
DOUBLE PRECISION, ALLOCATABLE           :: expintfunc(:)

SELECT CASE(approx)
CASE(1)
 eif = ap_1 - log10(var_x)
! not recommended to use

CASE(2)
 nintervals = 1e4
 ALLOCATE(expintfunc(nintervals))
 ! coefficients of linear points describing the interval division
 al = (const_pi / 2 - atan(var_x))/(nintervals - 1)
 bl = (nintervals * atan(var_x) - const_pi / 2) / (nintervals - 1)
 
 DO ind_I = 1, nintervals
  var_y = al * FLOAT(ind_I) + bl
  expintfunc(ind_I) = exp(-var_y)/atan(var_y) * 1/(1 + (var_y * var_y))
 END DO
 
 ! the trapezoid rule
 summ = 0
 DO ind_I = 1, nintervals - 1
  x1 = al * FLOAT(ind_I) + bl
  x2 = al * FLOAT(ind_I + 1) + bl
  loc_sum = (x2 - x1) * (expintfunc(ind_I + 1) - expintfunc(ind_I))
  summ = summ + loc_sum
 END DO
 
 eif = summ
CASE DEFAULT
 STOP 'exp_int_func: wrong approximation was chosen'
END SELECT

END SUBROUTINE exp_int_func
