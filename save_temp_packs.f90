!______________________________________________________
! save_temp_packs
! this SBR will save actually calculated packages which
! were caulculated in previous runs
!
!
SUBROUTINE save_temp_packs(pack_index)

USE types
USE constants
IMPLICIT NONE

INTEGER                         :: pack_index
INTEGER                         :: n_file
CHARACTER(100)                   :: temp_file_name
INTEGER                         :: I

n_file = INT(FLOAT(pack_index) / FLOAT(n_pack_save))

write(temp_file_name,"(A, A1, A, I3.3, A1, I6.6, A4)") TRIM(outputfolder), '/', TRIM(temp_filename), my_rank, "_", n_file, ".dat"
 
#if mpi!=1
 my_rank = 0
#endif
 ! write(*,*) "temp_packet.", my_rank, ".", n_file, temp_file_name

OPEN(73, form='unformatted', FILE=temp_file_name)
! OPEN(73, FILE=temp_file_name)
 ! write(*,*) 'save_temp_packs: I0 = ', n_pack_save * (n_file - 1) + 1,&
 !  ' I1 = ', n_pack_save * n_file
 DO I = n_pack_save * (n_file - 1) + 1, n_pack_save * n_file
  WRITE(73) package(I)
 END DO
CLOSE(73)

END SUBROUTINE save_temp_packs
