      SUBROUTINE SARGV(TEXT,I,ARGTEXT)

!**   Diese Funktion ist fast aequivalent zur entsprechenden C-Funktion.
!**   Sie ermittelt das i-te Argument in einem String.
!**   Die Regeln zur Argumenttrennung sind in sargp beschrieben.
!**   Ist das i-te Argument nicht vorhanden, wird argtext nicht veraendert.

      CHARACTER*(*) TEXT,ARGTEXT
      INTEGER I
      INTEGER N,AS,AE

      CALL SARGP(TEXT,N,I,AS,AE)

      IF (AS .EQ. -1) GOTO 10
      IF (AS .EQ. 0) THEN
             ARGTEXT=' '
      ELSE
             ARGTEXT=TEXT(AS:AE)
      ENDIF

10    RETURN

      END
