      SUBROUTINE SARGC(TEXT,N)

!   Diese Funktion ist aequivalent zur entsprechenden C-Funktion.
!   Sie ermittelt die Anzahl der Argumente in einem String.
!   Das Parsen des Strings uebernimmt die Routine sargp.
!   Dort sind auch die syntaktischen Regeln beschrieben.

      CHARACTER*(*) TEXT
      INTEGER AS,AE,N,I

      I=0 				! Do not look for any Argument
      CALL SARGP(TEXT,N,I,AS,AE)

      RETURN

      END
