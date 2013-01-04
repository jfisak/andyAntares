SUBROUTINE read_atomic_data()

! Read the file compose_data.DAT and remember
! elements(I)%ions(J)%ion_potential
! elements(I)%ions(J)%levels
! elements(I)%ions(J)%levels(n_Lev)
! elements(I)%ions(J)%levels(K)%exci_energy
! elements(I)%ions(J)%levels(K)%stat_waight
! ntransitions

  USE types

  IMPLICIT NONE    

  INTEGER          :: I, J, L, K, n_ions, atomicnumber, flag
  INTEGER          :: a_numb, i_stat, n_lev, l_numb, n_tran, remember_last_element, n_levels
  DOUBLE PRECISION :: i_pot, l_e, s_waight, ionstage, ionoffset


  OPEN (UNIT=8, FILE='atomic_data.dat')

  ntransitions = 0

  remember_last_element = 0

  ! Read data from atomic_data.dat until the end of the file and use only data
  ! of the elements and their ionization stages given in compose_adata.dat
  DO 
     ! Read only the first line of every element in particular ionization stage
     READ(8,*, END=10) a_numb, i_stat, n_lev, i_pot
     PRINT*, 'atomic data =', a_numb, i_stat, n_lev, i_pot
     ! In the case that the a_nub of the next reading line is not equal to the atomic number of 
     ! the last read element the last read element (donated with remember_last_element) then
     ! ionization potential will not change (or increase), e.i. ion.pot. increment ionoffset = 0.
     IF (a_numb .NE. remember_last_element) ionoffset = 0.
     ! Flag for the deciding for which elements of given ionization stage we will use data 
     ! For flag = 0. we will only read data but we will not use them, otherwise we will use them
     flag = 0
     ! Loop over all elements (which we already read from compose_adata.dat)
     DO I = 1, n_elements
        atomicnumber = elements(I)%atom_number
        n_ions = elements(I)%nions
!        PRINT*, 'n_ions=', n_ions
        ! Loop over all ions
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
      END DO

      ! In case we did not find the given element in the particular ionization stage
      ! we will just read data without using them. Print statement is only for testing
      ! if this case work properly when it should work
      IF (flag .EQ. 0) THEN
          DO K = 1, n_lev
             READ(8,*) l_numb, l_e, s_waight, n_tran
!             PRINT*, l_numb, l_e, s_waight, n_tran
          END DO
      END IF

      ! Increase the ionoffset by the pot.ene. i_pot of particular ion.stage
      ! For the neutral elements ionoffset = 0
      ionoffset = ionoffset + i_pot * e_v 
      ! Remember last element we pass
      remember_last_element = a_numb
  END DO
10 CONTINUE

  ! Real number of transitions (how many transitions are involved, using only 
  ! transition to the one direction, upward or downward)
  ntransitions = ntransitions/2.D0
  PRINT*, 'transition =', ntransitions

! Testing
  PRINT*, 'TESTING'
  DO I = 1, n_elements
     n_ions = elements(I)%nions
     PRINT*, 'element index is =', I, 'atomic number =',  elements(I)%atom_number, 'ionn stages =', elements(I)%nions
     DO J = 1, n_ions        
        PRINT*, 'ionisation index =', J,  'ionisation stage  =', elements(I)%ions(J)%ion_stage, &
                'ionisation pot. =', elements(I)%ions(J)%ion_potential/e_v,                     &
                'number of levels =', elements(I)%ions(J)%nlevels
        n_levels = elements(I)%ions(J)%nlevels
	DO K = 1, n_levels
	   PRINT*, K, elements(I)%ions(J)%levels(K)%exci_energy, elements(I)%ions(J)%levels(K)%stat_waight
	END DO
     END DO
  END DO

END SUBROUTINE read_atomic_data

