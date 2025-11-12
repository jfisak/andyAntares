!  Returns the number of non-blank characters of string TEXT
!
! 1x RETURN point
FUNCTION IDX(TEXT)
CHARACTER*(*) TEXT
INTEGER IDX
DO IDX=LEN(TEXT),1,-1
   ! RETURN point
   IF(TEXT(IDX:IDX).NE.' ') RETURN
END DO
IDX=0
RETURN
END
