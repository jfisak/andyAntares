SUBROUTINE warning(text)
IMPLICIT NONE

CHARACTER(LEN=*)              :: text

write(*,*) '****************************************************************'
write(*,*) text
write(*,*) '****************************************************************'

END SUBROUTINE
