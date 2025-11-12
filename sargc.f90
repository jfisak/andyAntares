! This function is equivalent to the corresponding C function. 
! It determines the number of arguments in a string. The string is
! parsed by the argp routine. The syntactic rules are also described
! there.
!
! INPUT: TEXT(CHAR): input text
!        Npoints(INT): number of points in the text
! OUTPUT: NONE
!
! 1x RETURN point
!
SUBROUTINE SARGC(TEXT,Npoints)

CHARACTER*(*) TEXT
INTEGER AS,AE,Npoints,ind_I

ind_I=0                               ! Do not look for any Argument
CALL SARGP(TEXT,Npoints,ind_I,AS,AE)

! RETURN point
RETURN

END
