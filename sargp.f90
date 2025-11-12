! The subroutine sargp parses a string.
! The following are determined:
! n : the number of arguments and, 
!         if i is in [1..n],
! as and ae : start and end index of the i-th argument.
! Parsing is performed according to the following rules:
! Spaces are generally ignored
! (see below for exceptions).
! Arguments are separated by spaces or commas or “=” or “:”.
! If an argument is separated by spaces and commas,
! this is considered a separation.
! If two commas follow without an argument, i.e. separated by spaces at most,
! this is considered an empty argument.
! Characters between two double quotation marks are considered
! Characters between two double quotation marks are considered
! to be one argument, even if they contain spaces or commas.
!  In this case, the quotation marks are
! not considered to be part of the argument.
! Interpretation of return values:
! n: always the number of arguments
! as=-1 -> No i-th argument was found (i not in [1..n])
! as=0  -> The i-th argument was empty (e.g. ‘,,’)
! otherwise : text(as:ae) = i-th argument
!
! INPUT: TEXT(CHR(*)): input text
!        Npoints(INT): number of arguments
!        AS(INT): position of the first argument
!        AE(INT): position of the last argument
! OUTPUT: ind_I(INT): desired argument 
! 
! 1x RETURN point
! 
SUBROUTINE SARGP(TEXT,Npoints,ind_I,AS,AE)

        IMPLICIT NONE

        CHARACTER*(*) TEXT
        INTEGER Npoints                  ! Out: Number of arguments
        INTEGER ind_I                  ! In : desired argument
        INTEGER AS,AE              ! first and last position of the
                                   ! argument in the string
        INTEGER TL,TI
        INTEGER STATE

        TL=LEN(TEXT)

        state=0
!	state=0 -> kein Argument aktiv
!	state=1 -> Argumentende gefunden, naechstes Komma
!			ist   k e i n   Leerargument
!	state=2 -> normales Argument aktiv
!	state=3 -> Argument in '"' aktiv


        Npoints=0
        AS=-1
        AE=TL

        DO 100 TI=1,TL

!** Go here looking for next start of an argument
        IF (STATE .EQ. 0) THEN
           IF (TEXT(TI:TI) .EQ. ' ') GOTO 100
           IF ((TEXT(TI:TI) .EQ. ',') &
           .OR.(TEXT(TI:TI) .EQ. '=') &
           .OR.(TEXT(TI:TI) .EQ. ':')) THEN
                            Npoints=Npoints+1
                            IF (Npoints .EQ. ind_I) THEN
                                AS=0
                            ENDIF
                            GOTO 100
           ENDIF
           IF (TEXT(TI:TI) .EQ. '"') THEN
                           Npoints=Npoints+1
                           IF (Npoints .EQ. ind_I) AS=TI+1
                           STATE=3
                           GOTO 100
           ENDIF
           STATE=2
           Npoints = Npoints + 1
           IF (Npoints .EQ. ind_I) AS=TI
           GOTO 100
        ELSEIF (STATE .EQ. 1) THEN
           IF (TEXT(TI:TI) .EQ. ' ') GOTO 100
           IF ((TEXT(TI:TI) .EQ. ',') &
           .OR.(TEXT(TI:TI) .EQ. '=') &
           .OR.(TEXT(TI:TI) .EQ. ':')) THEN
                           STATE=0
                           GOTO 100
           ENDIF
           IF (TEXT(TI:TI) .EQ. '"') THEN
                           Npoints=Npoints+1
                           IF (Npoints .EQ. ind_I) AS=TI+1
                           STATE=3
                           GOTO 100
           ENDIF
           STATE=2
           Npoints = Npoints + 1
           IF (Npoints .EQ. ind_I) AS=TI
           GOTO 100
        ELSEIF (STATE .EQ. 2) THEN
           IF (TEXT(TI:TI) .EQ. ' ') THEN
                          STATE=1
                          IF (Npoints .EQ. ind_I) AE=TI
                          GOTO 100
           ENDIF
           IF ((TEXT(TI:TI) .EQ. ',') &
           .OR.(TEXT(TI:TI) .EQ. '=') &
           .OR.(TEXT(TI:TI) .EQ. ':')) THEN
                          STATE=0
                          IF (Npoints .EQ. ind_I) AE=TI-1
                          GOTO 100
           ENDIF
           GOTO 100
        ELSE 
!***  ! IF (STATE .EQ. 3)
           IF (TEXT(TI:TI) .EQ. '"') THEN
                          STATE=1
                          IF (Npoints .EQ. ind_I) AE=TI-1
                          GOTO 100
           ENDIF
        ENDIF

100     CONTINUE

        ! RETURN point
        RETURN
        END
