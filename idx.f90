
       FUNCTION IDX(TEXT)
!  Returns the number of non-blank characters of string TEXT
       CHARACTER*(*) TEXT
       INTEGER IDX
       DO IDX=LEN(TEXT),1,-1
          IF(TEXT(IDX:IDX).NE.' ') RETURN
       END DO
       IDX=0
       RETURN
       END
