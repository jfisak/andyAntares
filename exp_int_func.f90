! this subroutine calculates exponential integral function for the given x
SUBROUTINE exp_int_func(approx, x, eif)
USE types
USE constants
IMPLICIT NONE
! input variable
INTEGER                                :: approx
DOUBLE PRECISION                        :: x
! output variable
DOUBLE PRECISION                        :: eif
! indexes
INTEGER                                 :: I
! number of intervals
INTEGER                                 :: nintervals
! parameters describing the interval division
DOUBLE PRECISION                        :: al, bl
! 
DOUBLE PRECISION                        :: y, x1, x2
DOUBLE PRECISION                        :: loc_sum, summ
! approximation constants
DOUBLE PRECISION, PARAMETER             :: ap_1 = -0.5772156649
! calculated values of expintfunc
DOUBLE PRECISION, ALLOCATABLE           :: expintfunc(:)

SELECT CASE(approx)
CASE(1)
 eif = ap_1 - log10(x)
! not recommended to use

CASE(2)
 nintervals = 1e4
 ALLOCATE(expintfunc(nintervals))
 ! coefficients of linear points describing the interval division
 al = (pi / 2 - atan(x))/(nintervals - 1)
 bl = (nintervals * atan(x) - pi / 2) / (nintervals - 1)
 
 DO I = 1, nintervals
  y = al * FLOAT(I) + bl
  expintfunc(I) = exp(-y)/atan(y) * 1/(1 + (y * y))
 END DO
 
 ! the trapezoid rule
 summ = 0
 DO I = 1, nintervals - 1
  x1 = al * FLOAT(I) + bl
  x2 = al * FLOAT(I + 1) + bl
  loc_sum = (x2 - x1) * (expintfunc(I + 1) - expintfunc(I))
  summ = summ + loc_sum
 END DO
 
 eif = summ
CASE DEFAULT
 STOP 'exp_int_func: wrong approximation was chosen'
END SELECT

END SUBROUTINE exp_int_func
