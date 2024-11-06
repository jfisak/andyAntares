SUBROUTINE find_unfinished_run()

USE types
USE constants
IMPLICIT NONE

CHARACTER(filename_length)          :: temp_file_name
INTEGER                         :: n_file, cur_pack
LOGICAL                         :: file_exists
INTEGER                         :: reading_packets
TYPE(photon)                    :: cur_package



n_file = 0
tot_saved_packets = 0
cur_pack = 0

DO
 n_file = n_file + 1
 write(temp_file_name,"(A, A1, A, I3.3, A1, I6.6, A4)") TRIM(outputfolder), '/', TRIM(temp_filename), my_rank, "_", n_file, ".dat"
 !write(temp_file_name,"(A, I3.3, A1, I6.6, A4)") TRIM(temp_filename), my_rank, "_", n_file, ".dat"
 inquire(FILE=temp_file_name, EXIST=file_exists)
 ! write(*,*) 'find_unfinished_run: my_rank = ', my_rank, ' n_file = ', n_file
 ! write(*,*) 'find_unfinished_run: filename = ', temp_file_name, ' E? = ', file_exists
 IF(file_exists) THEN
  OPEN(21, FORM="unformatted", FILE=temp_file_name)
  ! OPEN(21, FILE=temp_file_name)
   DO
    READ(21, iostat=reading_packets) cur_package
    ! write(*,*) cur_package
    IF(reading_packets /= 0) EXIT
    cur_pack = cur_pack + 1
    IF(cur_pack > SIZE(package) - 1) STOP 'number of loaded packets is larger than number of packets for a computation'
    package(cur_pack) = cur_package
    tot_saved_packets = tot_saved_packets + 1
   END DO
   ! write(*,*) 'find_unfinished_run: tot_saved_packets = ', tot_saved_packets
  CLOSE(21)
 END IF
 IF(.NOT. file_exists) EXIT
END DO
! write(*,*) 'find_unfinished_run: tot_saved_packets = ', tot_saved_packets


END SUBROUTINE find_unfinished_run

