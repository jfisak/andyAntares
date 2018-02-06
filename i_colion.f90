SUBROUTINE i_colion(approximation, indexe, indexi, act_level, pack_index, act_pop, Zion, &
        Zrecomb, Zintrecomb, actirates)
USE types
USE rates_i
IMPLICIT NONE

! input variables
INTEGER                         :: approximation
INTEGER                         :: pack_index, act_level
DOUBLE PRECISION                :: act_pop
INTEGER                         :: indexe, indexi
INTEGER                         :: nlevels
! grid informations
INTEGER                         :: current_mgi
DOUBLE PRECISION                :: el_dens, temp, gl_pop_ip1e, x
INTEGER                         :: nfreq
DOUBLE PRECISION                :: gindex
INTEGER                         :: get_package_model_index
INTEGER                         :: I, J
INTEGER                         :: npoints
DOUBLE PRECISION, PARAMETER     :: coll_const = 1.55D13
DOUBLE PRECISION, ALLOCATABLE   :: crossfreq(:)
DOUBLE PRECISION                :: eif, freq
! linear interpolation
DOUBLE PRECISION                :: ali, bli, freq1, freq2, func1, func2
INTEGER                         :: actPoint
DOUBLE PRECISION                :: cross_sect
DOUBLE PRECISION                :: exci_energy, gr_exci_energy, pop_number
DOUBLE PRECISION                :: stat_weight
! output variables
DOUBLE PRECISION                :: Zion, Zrecomb, Zintrecomb
! recombination
DOUBLE PRECISION, PARAMETER             :: times = 1e0
TYPE(irates)                   :: actirates



SELECT CASE(approximation)

!_______________________________________________________________
! cross section from file saved as a table
CASE (1)
 ! actual model grid index
 current_mgi = get_package_model_index(pack_index)
 ! electron density
 el_dens = model_grid(current_mgi)%e_dens
 ! temperature
 temp = model_grid(current_mgi)%T
 IF(indexi < elements(indexe)%atom_number) THEN
  nfreq = SIZE(elements(indexe)%ions(indexi)%levels(act_level)%photcros(1,:))
  IF(nfreq /= 0) THEN
   ALLOCATE(crossfreq(nfreq))
   ! number density of a ground state of ion indexi + 1, indexe
   ! frequency
   freq = (elements(indexe)%ions(indexi + 1)%levels(1)%exci_energy - &
           elements(indexe)%ions(indexi)%levels(act_level)%exci_energy) / h
   ! argument of E_1(x)
   x = (h * freq) / (BOLK * temp)
   ! exponential integral function (calculation of eif)
   CALL exp_int_func(1, x, eif)
   ! photoionization cross section
   !CALL bound_free_rates(pack_index, act_level, rad_rate)
   ! photoionization cross section for the given frequency freq
   crossfreq(:) = elements(indexe)%ions(indexi)%levels(act_level)%photcros(1,:)
   actPoint = 0
   DO I = 1, nfreq
    IF(freq > crossfreq(I)) THEN
     actPoint = I
     EXIT
    END IF
   END DO
   ! if the frequency is out of range of the frequency interval
   ! the total rate will be equal to zero
   IF(actPoint == 0 .OR. actPoint == 1) THEN
    Zion = 0
   ELSE
    ! now we have to do a linear interpolation between the points actPoint - 1 and actPoint
    freq1 = elements(indexe)%ions(indexi)%levels(act_level)%photcros(1, actPoint - 1)
    freq2 = elements(indexe)%ions(indexi)%levels(act_level)%photcros(1, actPoint)
    func1 = elements(indexe)%ions(indexi)%levels(act_level)%photcros(2, actPoint - 1)
    func2 = elements(indexe)%ions(indexi)%levels(act_level)%photcros(2, actPoint)
    ali = (func1 - func2) / (freq1 - freq2)
    bli = (func2 * freq1 - func1 * freq2) / (freq1 - freq2)
    cross_sect = ali * freq + bli
    !print*, 'photoionization: ali = ', ali, ' bli = ', bli, ' cross_sect = ', cross_sect
    ! gindex
    IF(indexi == 1) THEN
     gindex = 0.1
    ELSE IF(indexi == 2) THEN
     gindex = 0.2
    ELSE IF(indexi > 2) THEN
     gindex = 0.3
    END IF
    Zion = act_pop * el_dens *coll_const / temp**(1.0/2.0) * gindex * cross_sect * exp(-x) / x 
    Zion = times * Zion
    !print*, 'collion_rates: Zion = ', Zion
    ! write(*,*) 'i_colion: act_pop = ', act_pop, ' el_dens = ', el_dens, ' temp = ', temp, &
    !  ' gindex = ', gindex, ' cross_sect = ', cross_sect, ' exp(-x) = ', exp(-x)/x 
   END IF
   DEALLOCATE(crossfreq)
  ELSE ! npoints == 0
   Zion = 0.D0
  END IF
 ELSE ! indexi == atom number (atom is fully ionized)
  Zion = 0.D0
 END IF
 ! print*, 'collion_rates: Zion = ', Zion
!_________________________________________________________________________________________
! recombination
 SELECT CASE(nlte)
  CASE(0)
   Zrecomb = 0.D0
   Zintrecomb = 0.D0
   IF(indexi > 1) THEN
    nlevels = SIZE(actirates%Lma_int_reccol)
    !write(*,*) 'i_colion: number of points: ', nlevels
    DO I = 1, nlevels
     npoints = SIZE(elements(indexe)%ions(indexi - 1)%levels(I)%photcros(1,:))
     IF(npoints /= 0) THEN
      ALLOCATE(crossfreq(npoints))
      crossfreq(1:npoints) = elements(indexe)%ions(indexi - 1)%levels(I)%photcros(1, 1:npoints)
      ! frequency
      freq = (elements(indexe)%ions(indexi)%levels(act_level)%exci_energy - &
              elements(indexe)%ions(indexi - 1)%levels(I)%exci_energy) / h
      !write(*,*) 'i_colion: excienergy1 = ', &
      ! elements(indexe)%ions(indexi)%levels(act_level)%exci_energy, &
      ! ' excienergy2 =  ',  &
      ! elements(indexe)%ions(indexi - 1)%levels(I)%exci_energy, &
      ! 'freq = ', freq
      x = (h * freq) / (BOLK * temp)
      ! exponential integral function (calculation of eif)
      CALL exp_int_func(1, x, eif)
      actPoint = 0
      DO J = 1, npoints
       IF(freq < crossfreq(J)) THEN
        actPoint = J
        EXIT
       END IF
      END DO
      ! if the frequency is out of range of the frequency interval
      ! the total rate will be equal to zero
      !write(*,*) 'i_colion: actPoint = ', actPoint
      IF(actPoint == 0 .OR. actPoint == 1) THEN
       actirates%Lma_int_reccol(I) = 0.D0
       actirates%Lma_reccol(I) = 0.D0
       DEALLOCATE(crossfreq)
       CYCLE
      END IF
      ! now we have to do a linear interpolation between the points actPoint - 1 and actPoint
      freq1 = elements(indexe)%ions(indexi - 1)%levels(I)%photcros(1, actPoint - 1)
      freq2 = elements(indexe)%ions(indexi - 1)%levels(I)%photcros(1, actPoint)
      func1 = elements(indexe)%ions(indexi - 1)%levels(I)%photcros(2, actPoint - 1)
      func2 = elements(indexe)%ions(indexi - 1)%levels(I)%photcros(2, actPoint)
      ali = (func1 - func2) / (freq1 - freq2)
      bli = (func2 * freq1 - func1 * freq2) / (freq1 - freq2)
      cross_sect = ali * freq + bli
      ! gindex
      IF(indexi - 1 == 1) THEN
        gindex = 0.1
      ELSE IF(indexi - 1 == 2) THEN
        gindex = 0.2
      ELSE IF(indexi - 1 > 2) THEN
        gindex = 0.3
      END IF
      ! populations calculation
      CALL populations(indexe, indexi - 1, I, current_mgi, pop_number)
      exci_energy = elements(indexe)%ions(indexi - 1)%levels(I)%exci_energy
      gr_exci_energy = elements(indexe)%ions(indexi - 1)%levels(1)%exci_energy
      stat_weight = elements(indexe)%ions(indexi - 1)%levels(I)%stat_waight
      actirates%Lma_int_reccol(I) = pop_number * el_dens * coll_const / temp**(1.0/2.0) * gindex * cross_sect * &
        exp(-x) / x * exci_energy * stat_weight
      actirates%Lma_reccol(I) = pop_number * el_dens * coll_const / temp**(1.0/2.0) * gindex * cross_sect * &
        exp(-x) / x * (exci_energy - gr_exci_energy) * stat_weight
      ! for now it will be equal to zero
      !actirates%Lma_reccol(I) = 0.D0
      Zintrecomb = Zintrecomb + actirates%Lma_int_reccol(I)
      Zrecomb = Zrecomb + actirates%Lma_reccol(I)
      !print*, 'collion_rates: actirates%Lma_reccol = ', actirates%Lma_reccol(I)
      DEALLOCATE(crossfreq)
     END IF
     !write(*,*) 'i_colion: Zintrecomb = ', Zintrecomb, ' Zrecomb = ', Zrecomb
    END DO
   ELSE
    Zintrecomb = 0.D0
    Zrecomb = 0.D0
   END IF
   !print*, 'collion_rates: Zion = ', Zion, ' Zrecomb = ', Zrecomb
 END SELECT

CASE DEFAULT
 STOP 'bound_free_rates: non valid approximation was chosen'
END SELECT


END SUBROUTINE i_colion
