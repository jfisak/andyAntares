SUBROUTINE do_rpackage_event(pack_index, event)

  USE types

  IMPLICIT NONE    

  INTEGER                           :: pack_index, event
  DOUBLE PRECISION                  :: dist, rand_numb, ran2, tau, tau_rand
  DOUBLE PRECISION, DIMENSION(3)    :: direction


  IF (event .EQ. rpkt_eventtype_lineinteraction) THEN
     ! In this case the package interacts with a line. In the general
     ! case we now would make it an i-package and assign the i-package
     ! to the upper level of the corresponding bound-bound transition
     
     ! However, for now we restrict ourselves to the pure scattering
     ! case, i.e. we immediately re-emit the package isotropically in
     ! the cmf and conserve the cmf frequency print*,
     ! freq_line,package(pack_index)%freq_cmf
!     print*, 'package: ', pack_index, 'line interaction...'
     CALL emit_rpackage(pack_index)
  ELSE IF (event .EQ. rpkt_eventtype_continuum) THEN
     ! In this case the package undergoes a continuum event. In the
     ! general case we need to decide now if this was a
     ! electron-scattering, free-free or bound-free absorption (sample
     ! the different opacities) and act accordingly. In case of e/s
     ! the packet stays and an rpkt and is re-emitted isotropically in
     ! the cmf. For ff it becomes a kpkt, In the case of bf we have to
     ! check further if it will go to a kpkt or ipkt (bf contribute to
     ! both the thermal kinetic and internal energy pools).
     
     ! For now we set the electron number density to zero, so no
     !continuum opacity should be there print*, 'Continuum event
     !occurred. This should not happen for now!'
     
     ! As an easy next step make sure that a proper e/s opacity is
     ! calcualted, than the following lines should work and give
     ! isotropic re-emission in the cmf
!     print*, 'package: ', pack_index, 'continuum interaction...'
     CALL emit_rpackage(pack_index)
  ELSE
     STOP 'ERROR in do_rpackage event'
  END IF
  
END SUBROUTINE do_rpackage_event


!!!!!! This is for test analitic expresion !!!!!!!!
!!    print*, rand_numb
!    rand_numb = 0.1D0
!!    rand_numb = ran2(idum)     
!    IF (rand_numb .LT. 0.5D0) THEN
!        package(pack_index)%active = 0
!!!        print*, 'DO ABSORPTION'
!    ELSE
!        CALL emit_rpackage(pack_index)
!        tau = 0.D0
!10      rand_numb = ran2(idum)     
!        IF (rand_numb .EQ. 0.D0) GOTO 10    
!        tau_rand = -LOG(rand_numb)
!!!        print*, 'DO SCATTERING'
!    END IF
    


! For absorption only  

!     package(pack_index)%active = 0

!END
