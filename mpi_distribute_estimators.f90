! distribute estimators among processes
!
! INPUT: NONE
! OUTPUT: NONE
!
SUBROUTINE mpi_distribute_estimators()

USE MPI
USE types
USE constants
IMPLICIT NONE

INTEGER                                         :: ind_I
DOUBLE PRECISION, DIMENSION(n_modelgrid)        :: Jarray, recJarray
#if mpi==1
DO ind_I = 1, n_modelgrid
 Jarray(ind_I) = model_grid(ind_I)%J
 ! write(99,*) 'mpi_distribute_estimators: Jarray(ind_I) = ', Jarray(ind_I)
END DO
! radiation field in the cells
CALL MPI_REDUCE(Jarray, recJarray, n_modelgrid, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

IF(my_rank == 0) recJarray = recJarray / n_tasks

CALL MPI_BCAST(recJarray, n_modelgrid, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

DO ind_I = 1, n_modelgrid
 model_grid(ind_I)%J = recJarray(ind_I)
 ! write(99,*) 'mpi_distribute_estimators: J(ind_I) = ', model_grid(ind_I)%J
END DO
#endif

END SUBROUTINE mpi_distribute_estimators
