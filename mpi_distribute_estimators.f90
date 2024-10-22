SUBROUTINE mpi_distribute_estimators()

USE MPI
USE types
USE constants
IMPLICIT NONE

INTEGER                                         :: I
DOUBLE PRECISION, DIMENSION(n_modelgrid)        :: Jarray, recJarray
#if mpi==1
DO I = 1, n_modelgrid
 Jarray(I) = model_grid(I)%J
 ! write(99,*) 'mpi_distribute_estimators: Jarray(I) = ', Jarray(I)
END DO
! radiation field in the cells
CALL MPI_REDUCE(Jarray, recJarray, n_modelgrid, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

IF(my_rank == 0) recJarray = recJarray / n_tasks

CALL MPI_BCAST(recJarray, n_modelgrid, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

DO I = 1, n_modelgrid
 model_grid(I)%J = recJarray(I)
 ! write(99,*) 'mpi_distribute_estimators: J(I) = ', model_grid(I)%J
END DO
#endif

END SUBROUTINE mpi_distribute_estimators
