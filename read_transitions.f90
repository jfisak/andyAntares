SUBROUTINE read_transitions()

! Read the filetransitions.dat and remember
! linelist(L)%atom_number
! linelist(L)%ion
! linelist(L)%lower
! linelist(L)%upper
! linelist(L)%freq
! linelist(L)%A_ul
! linelist(L)%f_ul

  USE types

  IMPLICIT NONE    

  INTEGER                                       :: I, J, L, K, M, kindex, low_level, up_level, flag, atomicnumber
  INTEGER                                       :: a_numb, i_stat, n_tran, n_ions, ionstage
  DOUBLE PRECISION                              :: A, col_str, constant
  TYPE(line_list)                               :: dammy

!  OPEN (UNIT=9, FILE='transitions.dat')
  OPEN (UNIT=10, FILE='frequencies.dat')

  ! Calculate the constant for the oscilatior strength calculation
  constant = ( me_g * light_speed ** 3)/(8.D0 * pi ** 2 * e_charge**2)

  ALLOCATE (linelist(ntransitions))

  L = 1 

!  DO 
     PRINT*, 'transition data =', a_numb, i_stat, n_tran
     flag = 0

     DO I = 1, n_elements
        atomicnumber = elements(I)%atom_number
        OPEN (UNIT=9, status='old', FILE='transitions.dat')
        READ(9,*) a_numb, i_stat, n_tran
        n_ions = elements(I)%nions
        !PRINT*, I
        ! Loop over all ions
        DO J = 1, n_ions                   
          ionstage = elements(I)%ions(J)%ion_stage 
          !PRINT*, '  ',J
          IF ((a_numb .EQ. atomicnumber) .AND. (i_stat .EQ. ionstage)) THEN 
              DO K = 1, n_tran
                 READ(9,*) kindex, low_level, up_level, A, col_str
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
                 print*, 'L = ', L, ' linelist: ', linelist(L)
                 L = L+1  
              END DO
              flag = 1
           END IF

        END DO
       CLOSE(9)
      END DO


      IF (flag .EQ. 0) THEN
          DO K = 1, n_tran
                 READ(9,*) kindex, low_level, up_level, A, col_str
!                PRINT*, kindex, low_level, up_level, A, col_str
          END DO
      END IF
  
!  END DO
10 CONTINUE

  PRINT*, 'check lenght of line list', L-1, ntransitions
!  dammy = linelist(ntransitions)
!  print*, dammy
  
!  DO I = 1, ntransitions
!     PRINT*, I, linelist(I)%freq, linelist(I)%atom_number, linelist(I)%ion,linelist(I)%lower,linelist(I)%upper, linelist(I)%A_ul
!  END DO

  print*, 'n_transitions = ', ntransitions
  CALL sorting_new(ntransitions, linelist)
!  linelist(1) = dammy

!  DO I = 1, ntransitions
!     PRINT*, I, linelist(I)%freq, linelist(I)%atom_number, linelist(I)%ion,linelist(I)%lower,linelist(I)%upper, linelist(I)%A_ul
!  END DO
     
END SUBROUTINE read_transitions


