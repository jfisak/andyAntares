SUBROUTINE part_fun(indexe, indexi, temp, U)

  ! Calculate partition function of element indexe in th egiven ionization
  ! stage indexi at the given temperature temp
  USE types

  IMPLICIT NONE    

  INTEGER             :: indexe, indexi, indexl, e_gl, nlevels
  DOUBLE PRECISION    :: U, temp, g_level, e_level
 
  
!  print*, 'partition func. called for:', indexe, indexi, temp

  IF(.NOT. ALLOCATED(elements(indexe)%ions(indexi)%levels)) THEN
   write(*,*) 'part_fun: indexe = ', indexe, ' indexi = ', indexi
   write(*,*) 'levels are not allocated'
   STOP
  END IF
  U = elements(indexe)%ions(indexi)%levels(1)%stat_waight
  e_gl = elements(indexe)%ions(indexi)%levels(1)%exci_energy
!  print*, '  Part.func. initialisation:', U, e_gl

  ! Number ov levels for the given element indexe in ionisation stage indexi
  nlevels = SIZE(elements(indexe)%ions(indexi)%levels)

  DO indexl = 2, nlevels
     ! Statistical weight of th egrpund level
     g_level = elements(indexe)%ions(indexi)%levels(indexl)%stat_waight   
     ! Excitation energy of the excited level
     e_level = elements(indexe)%ions(indexi)%levels(indexl)%exci_energy
     ! Partition function
     IF(temp == 0 ) STOP 'part_fun: temperature = 0...'
     U = U + g_level * EXP(-(e_level - e_gl) / BOLK / temp)  
!     print*, '   part.func. calculation:', indexl, g_level, e_level/e_v
  END DO

!  print*, '   part.func. calculation done:', U

!  Only for testing
!  IF (!finite(U)) 
!    {
!      printout("element %d ion %d\n",element,ion);
!      printout("modelgridindex %d\n",modelgridindex);
!      printout("level %d, nlevels %d\n",level,nlevels);
!      printout("sw %g\n",stat_weight(element,ion,0));
!      printout("T_exc %g \n",T_exc);
!      abort();
!    }
 
END SUBROUTINE part_fun

