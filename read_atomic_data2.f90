! reads atomic data for the given element (from the file compose_adata.dat)
! in the format given by the parameter 
SUBROUTINE read_atomic_data2(element,lowerion,upperion,levels_type,filename)

  USE types

  IMPLICIT NONE    

 ! input parameters
 INTEGER                        :: element, lowerion, upperion, levels_type
 CHARACTER (LEN=20)             :: filename
 ! loop variables
 INTEGER                        :: I, J, K
 ! reading from file variables
 INTEGER                        :: ios, read_levels
 INTEGER                        :: current_element, current_ion, ions
 INTEGER                        :: n_levels, junk, l_numb
 CHARACTER (LEN=200)            :: line
 DOUBLE PRECISION               :: i_pot, l_energy, ionoffset, ionstage, s_weight
 ! basic setting of variables
 ionoffset = 0
 ions = 0
OPEN(8,status='old',FILE=filename)
 ! now we will choose the reading file using the given levels type
 SELECT CASE(levels_type)
 ! standard atomic data (26.03.2016)
 ! **0** the first input of the atomic data for the 3D wind code
 CASE(0)
  DO
   ! reading from the file filename
   ! we ignore rows starting *
   ! firstly we read the whole line
   READ(8,'(A)',IOSTAT=ios) line
   IF( INDEX(line, '*') /= 0) CYCLE
   IF (ios /= 0) THEN
    print*, 'the subroutine read_atomic data:'
    STOP 'ERROR: NO VALID ATOMIC DATA...'
   END IF
   ! if everything is OK, we will read from the variable line variables
   READ(line,*) current_element, current_ion, n_levels, i_pot
   print*, 'current element: ', current_element, 'current_ion: ', current_ion, &
        'number of levels: ', n_levels, 'i_pot: ', i_pot
   ! do we read the right file?
   IF((current_element /= element).OR. (current_ion < lowerion) .OR. (current_ion > upperion)) STOP 'WRONG ATOMIC DATA...'
   ! computation and allocation of important variables
   ionstage = elements(current_element)%ions(current_ion)%ion_stage
   elements(current_element)%ions(current_ion)%ion_potential = i_pot * e_v
   elements(current_element)%ions(current_ion)%nlevels = n_levels
   ALLOCATE(elements(current_element)%ions(current_ion)%levels(n_levels))
   ! now we will read every single atomic levels
   DO J=1, n_levels
    READ(8,*,IOSTAT=read_levels) l_numb, l_energy, s_weight, junk
    ! did we red anything?
    IF(read_levels /= 0) STOP 'WRONG NUMBER OF LEVELS IN THE FILE...'
    print*, 'from the level file: ', l_numb, l_energy, s_weight, junk
    ! increase the level energy l_e (multiply with e_v) by ionoffset i. e.
    ! for the ion. potential of the ground level
    elements(current_element)%ions(current_ion)%levels(J)%exci_energy = l_energy * e_v + ionoffset
    elements(current_element)%ions(current_ion)%levels(J)%stat_waight = s_weight
   END DO
    ! we can increase number of read ions :-)
    ions = ions + 1
    ! are every single ions already read?
    IF (ions == (upperion - lowerion + 1)) EXIT
  END DO
 ! **DEFAULT** default case
 CASE DEFAULT
  print*, 'choice has not been found'
  print*, 'reading only hydrogen data inluded in the code...'
  STOP 'with no valid atomic data...'
 END SELECT
CLOSE(8)
END SUBROUTINE read_atomic_data2

