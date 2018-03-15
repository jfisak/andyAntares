SUBROUTINE do_rpackage_event(pack_index, event, actirrates)

  USE types
  USE rates_r
  USE counters

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
  DOUBLE PRECISION                      :: D



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
     !$OMP ATOMIC
     linelist(package(pack_index)%last_line)%n_int = &
      linelist(package(pack_index)%last_line)%n_int + 1
     !$OMP ATOMIC
     count_r_line = count_r_line + 1
     ! in this sbr we get only excited states from the upper states
     package(pack_index)%typ = type_ipkt
     ! CALL emit_rpackage(pack_index)
     ! CALL doppler_factor(pack_index, D)
     ! D = 1.D0
     ! package(pack_index)%freq_rf = package(pack_index)%freq_cmf / D
     package(pack_index)%l_ele = linelist(package(pack_index)%last_line)%indexe
     package(pack_index)%l_ion = linelist(package(pack_index)%last_line)%indexi
     package(pack_index)%l_lev = linelist(package(pack_index)%last_line)%upper
     
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
   ! if(procout) write(*,*) 'do_rpackage_event: photon ', pack_index, ' continuum interaction...'
   Zthomson = actirrates%Lcont(4,1)
   Zphotion = 0.D0
   DO I = 2, n_photcrossect + 1
    Zphotion = Zphotion + actirrates%Lcont(4,I)
   END DO
   Zff = actirrates%Lcont(4,n_photcrossect + 2)
   ! write(*,*) 'do_rpackage_event: Zff = ', Zff
   ZcontTot = Zthomson + Zphotion + Zff

   ! DO I = 1, SIZE(actirrates%Lcont(1,:))
   !  write(*,*) 'do_rpackage_event: Lcont = ', actirrates%Lcont(4,I)
   ! END DO
   
   ! generating a random number
   rand = DBLE(random()) * ZcontTot
   ! write(*,*) 'do_rpackage_event: rand = ', rand, ' ZcontTot = ', ZcontTot, &
   !  ' Zthomson = ', Zthomson, ' Zphotion = ', Zphotion, ' Zff = ', Zff
   !______________________________________________________________________
   !_________________ ELECTRON SCATTERING ________________________________
   !______________________________________________________________________
   I = 1
   summ = 0.D0
   IF(rand >= summ .AND. rand <= Zthomson) THEN
    ! electron scattering occures
    ! changes only a direction of propagation
    !$OMP ATOMIC
    count_r_thom = count_r_thom + 1
      if(procout) write(*,*) 'do_rpackage_event: packet = ', pack_index, ' Thomson scattering'
    CALL emit_rpackage(pack_index)
    CALL doppler_factor(pack_index, D)
    package(pack_index)%freq_rf = package(pack_index)%freq_cmf / D
    package(pack_index)%e_rf = package(pack_index)%e_cmf / D
    RETURN
   END IF
   summ = actirrates%Lcont(4, I)
   !______________________________________________________________________
   !_____________________ PHOTOIONIZATION ________________________________
   !______________________________________________________________________
   summ = Zthomson
   ! write(*,*) 'do_rpackage_event: summ = ', summ
   IF(rand > summ .AND. rand <= Zphotion + summ) THEN
    DO I = 2, n_photcrossect + 1
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
         if(procout) write(*,*) 'do_rpackage_event: package = ', pack_index, ' photoionization -> i packet'
       !$OMP ATOMIC
       count_r_ph_i = count_r_ph_i + 1
       package(pack_index)%typ = type_ipkt
       ! we have to find corresponding transition for the do_ipacket sbr
       package(pack_index)%last_line = no_line
       package(pack_index)%l_ele = indexe
       package(pack_index)%l_ion = indexi + 1
       package(pack_index)%l_lev = 1
      ELSE
       if(procout) write(*,*) 'do_rpackage_event: package = ', pack_index, ' photoionization -> k packet'
       !$OMP ATOMIC
       count_r_ph_k = count_r_ph_k + 1
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
   ! write(*,*) 'do_rpackage_event: summ + Zff= ', summ + Zff
   IF(rand >= summ .AND. rand < summ + Zff) THEN
    package(pack_index)%typ = type_kpkt
    if(procout) write(*,*) 'do_rpackage_event: package = ', pack_index, ' free-free'
     !$OMP ATOMIC
     count_r_ff = count_r_ff + 1
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
