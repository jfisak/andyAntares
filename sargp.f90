      SUBROUTINE SARGP(TEXT,N,I,AS,AE)

!**	Die Subroutine sargp zerlegt einen String.
!**	Ermittelt werden:
!**	n : die Anzahl der Argumente und, 
!**     	falls i in [1..n],
!**	as und ae : Start- und Endindex des i-ten Arguments.
!**	Geparsed wird nach folgenden Regel:
!**	Leerzeichen werden grunsaetzlich nicht beachtet
!**	(Ausnahmen siehe unten).
!**	Argumente werden durch Leerzeichen oder Komma oder '=' oder ':'
!**     getrennt.
!**	Wird ein Argument durch Leerzeichen und Komma getrennt,
!**     so gilt dies als eine Trennung.
!**	Folgen zwei Kommata ohne Argument, also hoechstens durch
!**	Leerzeichen getrennt, gilt dies als Leerargument.
!**	Zeichen zwischen zwei doppelten Anfuehrungszeichen gelten
!**	als ein Argument, auch wenn Leerzeichen oder Kommata
!**	enthalten sind. In diesem Fall werden die Anfuehrungszeichen
!**	als nicht zum Argument gehoerig betrachtet.
!**	Interpretaion der Rueckgabewerte:
!**	n: immer die Anzahl der Argumente
!**	as=-1 -> Es wurde kein i-tes Argument gefunden (i nicht in [1..n])
!**	as=0  -> Das i-te Argument war leer (z.B. ",,")
!**	sonst : text(as:ae) = i-tes Argument

        CHARACTER*(*) TEXT
        INTEGER N                  ! Out: Anzahl der Argumente
        INTEGER I                  ! In : gesuchtes Argument
        INTEGER AS,AE              ! Out: erste u. letzte Position des 
                                   !         Arguments im String

        INTEGER TL,TI
        INTEGER STATE

        TL=LEN(TEXT)

        state=0
!	state=0 -> kein Argument aktiv
!	state=1 -> Argumentende gefunden, naechstes Komma
!			ist   k e i n   Leerargument
!	state=2 -> normales Argument aktiv
!	state=3 -> Argument in '"' aktiv


        N=0
        AS=-1
        AE=TL

        DO 100 TI=1,TL

!** Go here looking for next start of an argument
        IF (STATE .EQ. 0) THEN
           IF (TEXT(TI:TI) .EQ. ' ') GOTO 100
           IF ((TEXT(TI:TI) .EQ. ',') &
           .OR.(TEXT(TI:TI) .EQ. '=') &
           .OR.(TEXT(TI:TI) .EQ. ':')) THEN
                            N=N+1
                            IF (N .EQ. I) THEN
                                AS=0
                            ENDIF
                            GOTO 100
           ENDIF
           IF (TEXT(TI:TI) .EQ. '"') THEN
                           N=N+1
                           IF (N .EQ. I) AS=TI+1
                           STATE=3
                           GOTO 100
           ENDIF
           STATE=2
           N=N+1
           IF (N .EQ. I) AS=TI
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
                           N=N+1
                           IF (N .EQ. I) AS=TI+1
                           STATE=3
                           GOTO 100
           ENDIF
           STATE=2
           N=N+1
           IF (N .EQ. I) AS=TI
           GOTO 100
        ELSEIF (STATE .EQ. 2) THEN
           IF (TEXT(TI:TI) .EQ. ' ') THEN
                          STATE=1
                          IF (N .EQ. I) AE=TI
                          GOTO 100
           ENDIF
           IF ((TEXT(TI:TI) .EQ. ',') &
           .OR.(TEXT(TI:TI) .EQ. '=') &
           .OR.(TEXT(TI:TI) .EQ. ':')) THEN
                          STATE=0
                          IF (N .EQ. I) AE=TI-1
                          GOTO 100
           ENDIF
           GOTO 100
        ELSE 
!***  ! IF (STATE .EQ. 3)
           IF (TEXT(TI:TI) .EQ. '"') THEN
                          STATE=1
                          IF (N .EQ. I) AE=TI-1
                          GOTO 100
           ENDIF
        ENDIF

100     CONTINUE

        RETURN
        END
