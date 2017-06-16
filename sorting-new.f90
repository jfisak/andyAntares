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


SUBROUTINE sorting_new(N, ARR)

  USE types

  IMPLICIT NONE   

  INTEGER                                           :: N, I, J
  DOUBLE PRECISION                                  :: A  
  TYPE(line_list)                                   :: dummy
  TYPE(line_list), DIMENSION(N)                     :: ARR

  DO J = 2, N
     I = J - 1
     
!     dummy = ARR(J) 
     A = ARR(J)%freq
!     print*, a, arr(j)
!     DO WHILE (I .GE. 1 .AND. ARR(I)%freq .GT. A)
     DO WHILE (I .GE. 1)
      IF(ARR(I)%freq .GT. A) THEN
       dummy = ARR(I+1)
       ARR(I+1) = ARR(I) 
       ARR(I) = dummy
!       print*, arr(i)
!       print*, arr(i+1)
      END IF
       I = I - 1
     END DO
!     ARR(I+1) = dummy
  END DO

END SUBROUTINE sorting_new
