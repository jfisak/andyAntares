!*****************************************************
!* Sorts an array ARR of length N in ascending order *
!* by straight insertion.                            *
!* ------------------------------------------------- *
!* INPUTS:                                           *
!*	    N	  size of table ARR                  *
!*          ARR	  table to be sorted                 *
!* OUTPUT:                                           *
!*	    ARR   table sorted in ascending order    *
!*                                                   *
!* NOTE: Straight insertion is a N² routine and      *
!*       should only be used for relatively small    *
!*       arrays (N<100).                             *
!*****************************************************         
SUBROUTINE sorting(N, ARR)

  USE types

  IMPLICIT NONE   

  INTEGER                           :: N, I, J
  DOUBLE PRECISION                  :: A
  DOUBLE PRECISION , DIMENSION(N)   :: ARR  

  DO J = 2, N
    A=ARR(J)
    DO I = J-1, 1, -1
       IF (ARR(I) .LE. A) GOTO 10
       ARR(I+1)=ARR(I)
    END DO
    I=0
10  ARR(I+1)=ARR(J)
  END DO
  RETURN

END SUBROUTINE sorting

