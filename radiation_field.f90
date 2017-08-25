! this sbr computes a radiation field J(\nu) for the given cell
! there are (not now) several approximations how to do that
SUBROUTINE radiation_field(current_mgi)
USE types
IMPLICIT NONE

! input variables
INTEGER                         :: current_mgi

SELECT CASE(approx)

CASE(0)
 
CASE DEFAULT
 STOP 'radiation_field: this approximation is not known...'
END SELECT

END SUBROUTINE radiation_field
