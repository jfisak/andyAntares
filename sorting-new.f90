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
!SUBROUTINE sorting_new(N, ARR)

!  USE types

!  IMPLICIT NONE   

!  INTEGER                                           :: N, I, J
!  DOUBLE PRECISION                                  :: A  
!  DOUBLE PRECISION, DIMENSION(N)                    :: ARR

!  DO J = 2, N
!     I = J - 1
!     A = ARR(J)
!     DO WHILE (I .GE. 1 .AND. ARR(I) .GT. A)
!       ARR(I+1) = ARR(I) 
!       I = I - 1
!     END DO
!     ARR(I+1) = A
!  END DO

!END SUBROUTINE sorting_new


SUBROUTINE sorting_new(n_points, ARR)

  USE types
USE constants

  IMPLICIT NONE   

  INTEGER                                           :: n_points, ind_I, ind_J
  DOUBLE PRECISION                                  :: var_A
  TYPE(line_list)                                   :: dummy
  TYPE(line_list), DIMENSION(n_points)                     :: ARR

  DO ind_J = 2, n_points
     ind_I = ind_J - 1
     
!     dummy = ARR(J) 
     var_A = ARR(ind_J)%freq
!     print*, a, arr(j)
!     DO WHILE (I .GE. 1 .AND. ARR(I)%freq .GT. A)
     DO WHILE (ind_I .GE. 1)
      IF(ARR(ind_I)%freq .LT. var_A) THEN
       dummy = ARR(ind_I+1)
       ARR(ind_I+1) = ARR(ind_I) 
       ARR(ind_I) = dummy
!       print*, arr(i)
!       print*, arr(i+1)
      END IF
       ind_I = ind_I - 1
     END DO
!     ARR(I+1) = dummy
  END DO

END SUBROUTINE sorting_new
