SUBROUTINE do_rpackage_event(pack_index, event, actirrates)

  USE types
  USE rates_r

  IMPLICIT NONE    

  INTEGER                           :: pack_index, event
  DOUBLE PRECISION                  :: dist, rand_numb, tau, tau_rand
  DOUBLE PRECISION, DIMENSION(3)    :: direction
  ! loop variables
  INTEGER                               :: I
  DOUBLE PRECISION                      :: summ, rand
  TYPE(rrates)                          :: actirrates
  REAL(8)                               :: random
  DOUBLE PRECISION                      :: freq, freqt
  INTEGER                               :: indexe, indexi, indexl
  INTEGER                               :: actIndex
  INTEGER                               :: n_ions, n_levels
  INTEGER                               :: nline
  ! total rates for the given processes
  DOUBLE PRECISION                      :: Zthomson, Zphotion, Zff
  DOUBLE PRECISION                      :: ZcontTot
  LOGICAL                               :: procout = .FALSE.



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
     ! in this sbr we get only excited states from the upper states
     isUpperTransition = .TRUE.
     package(pack_index)%typ = type_ipkt
      if(procout) write(*,*) 'photon ', pack_index, ' line interaction...'
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
   package(pack_index)%n_interactions = package(pack_index)%n_interactions + 1
   ! total number of continuum rates
    if(procout) write(*,*) 'do_rpackage_event: photon ', pack_index, ' continuum interaction...'
   Zff = 0.D0
   Zthomson = actirrates%Lcont(4,1)
   Zphotion = 0.D0
   DO I = 2, n_photcrossect + 1
    Zphotion = Zphotion + actirrates%Lcont(4,I)
   END DO
   DO I = 2 + n_photcrossect , 1 + n_photcrossect + n_ff
    Zff = Zff + actirrates%Lcont(4,I)
   END DO
   ZcontTot = Zthomson + Zphotion + Zff
   
   ! just only for now
   DO I = 1, n_tot_cont ! three columns, we are interested in nof rows
    ZcontTot = ZcontTot + actirrates%Lcont(4, I)
    ! write(*,*) 'do_rpackage_event: I = ', I, ' actirrates%Lcont = ', actirrates%Lcont(4,I)
   END DO
   ! generating a random number
   rand = DBLE(random()) * ZcontTot
   ! write(*,*) 'do_rpackage_event: rand = ', rand, ' ZcontTot = ', ZcontTot
   !______________________________________________________________________
   !_________________ ELECTRON SCATTERING ________________________________
   !______________________________________________________________________
   I = 1
   summ = 0.D0
   IF(rand >= summ .AND. rand <= Zthomson) THEN
    ! electron scattering occures
    ! changes only a direction of propagation
    count_thomson = count_thomson + 1
      if(procout) write(*,*) 'do_rpackage_event: Thomson scattering'
    CALL emit_rpackage(pack_index)
    RETURN
   END IF
   summ = actirrates%Lcont(4, I)
   !______________________________________________________________________
   !_____________________ PHOTOIONIZATION ________________________________
   !______________________________________________________________________
   summ = Zthomson
   IF(rand > summ .AND. rand <= Zphotion + summ) THEN
    DO I = 2, n_photcrossect + 1! three columns, we are interested in nof rows
     IF(rand >= summ .AND. rand <= actirrates%Lcont(4, I) + summ) THEN
      ! we have to choose if the packet transofrms onto i or k packet
      ! we will get it from the treshold frequency for the given ion
      indexe = actirrates%Lcont(1,I)
      indexi = actirrates%Lcont(2,I)
      indexl = actirrates%Lcont(3,I)
      ! treshold frequency
      freqt = elements(indexe)%ions(indexi)%levels(indexl)%phfreq
      freq = package(pack_index)%freq_cmf
      rand = DBLE(random())
      IF(rand < freqt / freq) THEN
         if(procout) write(*,*) 'do_rpackage_event: photoionization -> i packet'
       package(pack_index)%typ = type_ipkt
       ! we have to find corresponding transition for the do_ipacket sbr
       DO nline = 1, ntransitions
        IF(indexe == linelist(nline)%indexe .AND. &
           indexi == linelist(nline)%indexi) THEN
         IF(indexl == linelist(nline)%upper) THEN
          package(pack_index)%last_line = nline
          isUpperTransition = .TRUE.
          EXIT
         ELSE IF(indexl == linelist(nline)%lower) THEN
          package(pack_index)%last_line = nline
          isUpperTransition = .FALSE.
          EXIT
         ! test for level number
         END IF
        ! test for indexe and indexi
        END IF
       ! loop over lines
       END DO
      ELSE
       if(procout) write(*,*) 'do_rpackage_event: photoionization -> k packet'
       package(pack_index)%typ = type_kpkt
      END IF
      EXIT
     END IF
     summ = summ + actirrates%Lcont(4,I)
    END DO
   END IF
   !______________________________________________________________________
   !_____________________ FREE-FREE PROCESS ______________________________
   !______________________________________________________________________
   summ = summ + Zphotion
   IF(rand >= summ .AND. rand < summ + Zff) THEN
    package(pack_index)%typ = type_kpkt
     if(procout) write(*,*) 'do_rpackage_event: free-free'
   END IF
   !______________________________________________________________________
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
