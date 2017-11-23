SUBROUTINE do_rpackage_event(pack_index, event, Lcont)

  USE types

  IMPLICIT NONE    

  INTEGER                           :: pack_index, event
  DOUBLE PRECISION                  :: dist, rand_numb, ran2, tau, tau_rand
  DOUBLE PRECISION, DIMENSION(3)    :: direction
  ! loop variables
  INTEGER                               :: I
  DOUBLE PRECISION                      :: summ, rand, ZcontTot
  DOUBLE PRECISION, POINTER               :: Lcont(:)


  IF (event .EQ. rpkt_eventtype_lineinteraction) THEN
     ! In this case the package interacts with a line. In the general
     ! case we now would make it an i-package and assign the i-package
     ! to the upper level of the corresponding bound-bound transition
     
     ! However, for now we restrict ourselves to the pure scattering
     ! case, i.e. we immediately re-emit the package isotropically in
     ! the cmf and conserve the cmf frequency print*,
     ! freq_line,package(pack_index)%freq_cmf
!     print*, 'photon ', pack_index, ' line interaction...'
     package(pack_index)%n_interactions = package(pack_index)%n_interactions + 1
     package(pack_index)%typ = type_ipkt
     !print*, 'photon ', pack_index, ' line interaction...'
!     CALL emit_rpackage(pack_index)
  ELSE IF (event .EQ. rpkt_eventtype_continuum) THEN
     ! In this case the package undergoes a continuum event. In the
     ! general case we need to decide now if this was a
     ! electron-scattering, free-free or bound-free absorption (sample
     ! the different opacities) and act accordingly. In case of e/s
     ! the packet stays and an rpkt and is re-emitted isotropically in
     ! the cmf. For ff it becomes a kpkt, In the case of bf we have to
     ! check further if it will go to a kpkt or ipkt (bf contribute to
     ! both the thermal kinetic and internal energy pools).
   ! total number of continuum rates
   ZcontTot = 0.D0
   package(pack_index)%n_interactions = package(pack_index)%n_interactions + 1
   DO I = 1, SIZE(Lcont)
    ZcontTot = ZcontTot + Lcont(I)
    !write(*,*) 'do_rpackage_event: I = ', I, ' Lcont = ', Lcont(I)
   END DO
   ! generating a random number
   rand = ran2(idum) * ZcontTot
   summ = 0.D0
   IF(rand >= summ .AND. rand <= Lcont(1) + summ) THEN
    ! electron scattering occures
    ! changes only a direction of propagation
    count_thomson = count_thomson + 1
    CALL emit_rpackage(pack_index)
    !write(*,*) 'do_rpackage_event: electron scattering'
    RETURN
   END IF
   summ = Lcont(1)
   ! photoionization
   DO I = 2, SIZE(Lcont)
    IF(rand >= summ .AND. rand <= Lcont(I) + summ) THEN
     ! temporary solution
     !IF(rand >= summ .AND. rand <= Lcont(I)/2.D0 + summ) THEN
      package(pack_index)%typ = type_ipkt
      !write(*,*) 'do_rpackage_event: packet = ', pack_index, ' b-f process'
     !ELSE
     !  package(pack_index)%typ = type_ipkt
     EXIT
    END IF
    summ = summ + Lcont(I)
   END DO

     
     
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
