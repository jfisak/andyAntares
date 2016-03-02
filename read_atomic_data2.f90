! simple subroutine for the read atomic data from the multiple input files
SUBROUTINE read_atomic_data()

  USE types

  IMPLICIT NONE    

  INTEGER          :: I, J, L, K, n_ions, atomicnumber, flag
  INTEGER          :: a_numb, i_stat, n_lev, l_numb, n_tran, remember_last_element, n_levels
  DOUBLE PRECISION :: i_pot, l_e, s_waight, ionstage, ionoffset
  CHARACTER(20)    :: elfile

  ntransitions = 0
  remember_last_element = 0
! every element will be red from the given file
 DO I = 1,n_elements
  ! name of the file containing energy level informations
  elfile=elements(I)%levelfile
  PRINT*, 'levelfile ', elfile
  OPEN(UNIT=8,FILE=elfile)
   READ(8,*) a_numb, i_stat, n_lev, i_pot
   atomicnumber = elements(I)%atom_number
   n_ions = elements(I)%nions
   ! loop over all ions
        DO J = 1, n_ions                   
          ionstage = elements(I)%ions(J)%ion_stage 
          ! If in the file atomic_data.dat we found particular element in the given atomic stage 
          ! wish we have in the compose_adata.dat list, this condition will be done
          IF ((a_numb .EQ. atomicnumber) .AND. (i_stat .EQ. ionstage)) THEN 
              ! Associate the pot.energy i_pot (multiple with e_v) from compose_adata.dat
              ! to the elements(I)%ions(J)%ion_potential
              elements(I)%ions(J)%ion_potential = i_pot * e_v
              ! Associate the numb.lavel n_lev from compose_adata.dat
              ! to the elements(I)%ions(J)%nlevels
              elements(I)%ions(J)%nlevels = n_lev
!              print*, I, J, n_lev
              ! Allocate the memory; define how large should be elements(I)%ions(J)%levels (n_lev)
              ALLOCATE (elements(I)%ions(J)%levels(n_lev))
              ! Loop over all levels 
              DO K = 1, n_lev
                 READ(8,*) l_numb, l_e, s_waight, n_tran
!                PRINT*, l_numb, l_e, s_waight, n_tran
                 ! Increase the level energy l_e (multiply with e_v) by ionoffset i.e. 
                 ! for the ion.pot. of the ground level
                 elements(I)%ions(J)%levels(K)%exci_energy = l_e * e_v + ionoffset
!                IF (K .EQ. 2) PRINT*,  l_e * e_v + ionoffset
                 elements(I)%ions(J)%levels(K)%stat_waight = s_waight
                 ! Calculate total number of transitions of the given level
                 ntransitions =  ntransitions + n_tran
              END DO
              ! When the all levels of the given element in the particular ionization stage
              ! is used, then put flag to 1. It means that we will continue to read file
              ! atomic_data.dat from the next element
              flag = 1
           END IF
        END DO
  CLOSE(8)
 END DO

10 CONTINUE

  ! Real number of transitions (how many transitions are involved, using only 
  ! transition to the one direction, upward or downward)
  ntransitions = ntransitions/2.D0
  PRINT*, 'transition =', ntransitions

END SUBROUTINE read_atomic_data
