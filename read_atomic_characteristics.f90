SUBROUTINE read_atomic_characteristics()

  USE types
! Read the file compose_data.DAT and remember
! elements(I)%ions(J)%ion_potential
! elements(I)%ions(J)%levels
! elements(I)%ions(J)%levels(n_Lev)
! elements(I)%ions(J)%levels(K)%exci_energy
! elements(I)%ions(J)%levels(K)%stat_waight
! ntransitions


  IMPLICIT NONE    

  INTEGER          :: I, J, L, K, n_ions, atomicnumber, flag
  INTEGER          :: a_numb, i_stat, n_lev, l_numb, n_tran, remember_last_element, n_levels
  DOUBLE PRECISION :: i_pot, l_e, s_waight, ionstage, ionoffset
  CHARACTER(20)    :: elfile, junk
  ! for atomic transitions
  INTEGER          :: M, kindex, low_level, up_level
  DOUBLE PRECISION :: A, col_str, constant
  TYPE(line_list)  :: dammy


!  OPEN (UNIT=8, FILE='atomic_data.dat')

  ntransitions = 0

DO I = 1, n_elements
!________________________________________________________________
!___________ firstly we read atomic levels ______________________
!________________________________________________________________
 ionoffset = 0
 elfile=elements(I)%levelfile
 PRINT*, 'levelfile ', elfile
 OPEN(UNIT=8,status='old',FILE=elfile)
   atomicnumber = elements(I)%atom_number
   n_ions = elements(I)%nions
   ! Loop over all ions
   DO J = 1, n_ions
    READ(8,*) a_numb, i_stat, n_lev, i_pot
    PRINT*, 'i_pot=', i_pot, 'n_ions = ', n_ions
     ionstage = elements(I)%ions(J)%ion_stage 
     ! If in the file atomic_data.dat we found particular element in the given atomic stage 
     ! wish we have in the compose_adata.dat list, this condition will be done
     !IF ((a_numb .EQ. atomicnumber) .AND. (i_stat .EQ. ionstage)) THEN 
         ! Associate the pot.energy i_pot (multiple with e_v) from compose_adata.dat
         ! to the elements(I)%ions(J)%ion_potential
         elements(I)%ions(J)%ion_potential = i_pot * e_v
         ! Associate the numb.lavel n_lev from compose_adata.dat
         ! to the elements(I)%ions(J)%nlevels
         elements(I)%ions(J)%nlevels = n_lev
         print*, 'I = ', I, 'J = ', J, 'n_lev = ', n_lev
         ! Allocate the memory; define how large should be elements(I)%ions(J)%levels (n_lev)
         ALLOCATE (elements(I)%ions(J)%levels(n_lev))
         ! Loop over all levels 
         DO K = 1, n_lev
            READ(8,*) l_numb, l_e, s_waight, n_tran
            PRINT*, l_numb, l_e, s_waight, n_tran
            ! Increase the level energy l_e (multiply with e_v) by ionoffset i.e. 
            ! for the ion.pot. of the ground level
            elements(I)%ions(J)%levels(K)%exci_energy = l_e * e_v + ionoffset
         !   IF (K .EQ. 2) PRINT*,  l_e * e_v + ionoffset
            elements(I)%ions(J)%levels(K)%stat_waight = s_waight
            ! Calculate total number of transitions of the given level
            ntransitions =  ntransitions + n_tran
         END DO
         ! When the all levels of the given element in the particular ionization stage
         ! is used, then put flag to 1. It means that we will continue to read file
         ! atomic_data.dat from the next element
         flag = 1
      !END IF
! this part will be soon erased, I hope
   END DO
   READ(8,*) junk
   print*, 'number of transitions: ', ntransitions/2
   IF(I .EQ. n_elements) THEN
    ntransitions = ntransitions/2.D0
    ALLOCATE (linelist(ntransitions))
    L=1
   END IF
!________________________________________________________________
!_______ now we read atomic lines _______________________________
!________________________________________________________________
   !     atomicnumber = elements(I)%atom_number
   !    elfile=elements(I)%transitionfile
   !     OPEN (UNIT=9, status='old', FILE=elfile)
   !     n_ions = elements(I)%nions
        !PRINT*, I
        ! Loop over all ions
        DO J = 1, n_ions                   
         READ(8,*) a_numb, i_stat, n_tran
   !       ionstage = elements(I)%ions(J)%ion_stage 
          !PRINT*, '  ',J
   !       IF ((a_numb .EQ. atomicnumber) .AND. (i_stat .EQ. ionstage)) THEN 
              DO K = 1, n_tran
                 READ(8,*) kindex, low_level, up_level, A, col_str
                 PRINT*, kindex, low_level, up_level, A, col_str
!                 linelist(L)%atom_number = a_numb
!                 linelist(L)%ion = ionstage
                 linelist(L)%indexe = I
                 linelist(L)%indexi = J
                 linelist(L)%lower = low_level
                 linelist(L)%upper = up_level
                 linelist(L)%freq = (elements(I)%ions(J)%levels(up_level)%exci_energy - &
                                     elements(I)%ions(J)%levels(low_level)%exci_energy) / h
                 linelist(L)%A_ul = A
                 linelist(L)%f_ul = constant * (elements(I)%ions(J)%levels(up_level)%stat_waight / &
                                    elements(I)%ions(J)%levels(low_level)%stat_waight) *           &
                                    (linelist(L)%A_ul / linelist(L)%freq ** 2)
                 L = L+1  
              END DO
!           END IF

        END DO
!       CLOSE(9)
     CLOSE(8)
     print*, 'uzavirani souboru ', elfile
END DO

  PRINT*, 'check lenght of line list', L-1, ntransitions
  print*, 'n_transitions = ', ntransitions
  CALL sorting_new(ntransitions, linelist)

!      IF (flag .EQ. 0) THEN
!          DO K = 1, n_tran
!                 READ(9,*) kindex, low_level, up_level, A, col_str
!!                PRINT*, kindex, low_level, up_level, A, col_str
!          END DO
!      END IF


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
!      remember_last_element = a_numb
!      II=II+1
!END DO
!  END DO
!10 CONTINUE

  ! Real number of transitions (how many transitions are involved, using only 
  ! transition to the one direction, upward or downward)
!  ntransitions = ntransitions/2.D0
!  PRINT*, 'transition =', ntransitions
!
!! Testing
!  PRINT*, 'ATOMIC DATA: TESTING'
!  DO I = 1, n_elements
!     print*, elements(I)%levelfile
!     n_ions = elements(I)%nions
!     PRINT*, 'element index is =', I, 'atomic number =',  elements(I)%atom_number, 'ionn stages =', elements(I)%nions
!     DO J = 1, n_ions        
!        PRINT*, 'ionisation index =', J,  'ionisation stage  =', elements(I)%ions(J)%ion_stage, &
!                'ionisation pot. =', elements(I)%ions(J)%ion_potential/e_v,                     &
!                'number of levels =', elements(I)%ions(J)%nlevels
!        n_levels = elements(I)%ions(J)%nlevels
!      DO K = 1, n_levels
!         PRINT*, K, elements(I)%ions(J)%levels(K)%exci_energy, elements(I)%ions(J)%levels(K)%stat_waight
!      END DO
!     END DO
!  END DO

END SUBROUTINE read_atomic_characteristics
