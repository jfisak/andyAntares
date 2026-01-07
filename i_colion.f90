! ionization collisional rates for i-packets 
! 
! INPUT: approximation(INT): obvious
!        indexe(INT): element index
!        indexi(INT): ion index
!        act_level(INT): an index of the current level
!        pack_index(INT): index of the packet
!        act_pop(DBLE): a population in the current level
! OUPTUT: Zion(DBLE): the total ionization rate
!         Zintrecomb(DBLE): the total rate of internal recombination
!         Zrecomb(DBLE): the total rate of recombination
!         actirates(irates): the rates for corresponding transitions
!
SUBROUTINE i_colion(approximation, indexe, indexi, act_level, pack_index, act_pop, Zion, &
        Zintrecomb, Zrecomb, actirates)
USE types
USE constants
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
DOUBLE PRECISION                :: el_dens, temp, var_x
INTEGER                         :: nfreq
DOUBLE PRECISION                :: gindex
INTEGER                         :: get_package_model_index
INTEGER                         :: ind_I, ind_J
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


Zion = 0.D0
Zrecomb = 0.D0
Zintrecomb = 0.D0
write(*,*) 'i_colion: Zion = ', Zion, ' Zrecomb = ', Zrecomb, ' Zintrecomb = ', Zintrecomb

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
 IF(indexi <= elements(indexe)%atom_number) THEN
  IF(ALLOCATED(elements(indexe)%ions(indexi)%levels(act_level)%photcros)) THEN
   nfreq = SIZE(elements(indexe)%ions(indexi)%levels(act_level)%photcros(1,:))
  ELSE
   nfreq = 0
  END IF
  IF(nfreq /= 0) THEN
   ALLOCATE(crossfreq(nfreq))
   ! number density of a ground state of ion indexi + 1, indexe
   ! frequency
   freq = (elements(indexe)%ions(indexi + 1)%levels(1)%exci_energy - &
           elements(indexe)%ions(indexi)%levels(act_level)%exci_energy) / const_h
   ! argument of E_1(var_x)
   var_x = (const_h * freq) / (const_kB * temp)
   ! exponential integral function (calculation of eif)
   CALL exp_int_func(1, var_x, eif)
   ! photoionization cross section
   !CALL bound_free_rates(pack_index, act_level, rad_rate)
   ! photoionization cross section for the given frequency freq
   crossfreq(:) = elements(indexe)%ions(indexi)%levels(act_level)%photcros(1,:)
   actPoint = 0
   DO ind_I = 1, nfreq
    IF(freq > crossfreq(ind_I)) THEN
     actPoint = ind_I
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
    Zion = act_pop * el_dens *coll_const / temp**(1.0/2.0) * gindex * cross_sect * exp(-var_x) / var_x 
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
  ! print*, 'i_colion: Zion = ', Zion
!_________________________________________________________________________________________
! recombination
 SELECT CASE(nlte)
  CASE(0)
   actirates%Lma_int_reccol(:) = 0.D0
   actirates%Lma_reccol(:) = 0.D0
   write(*,*) 'i_colion: indexi = ', indexi
   IF(indexi > 1) THEN
    nlevels = SIZE(actirates%Lma_int_reccol)
    !write(*,*) 'i_colion: number of points: ', nlevels
    DO ind_I = 1, nlevels
     IF(ALLOCATED(elements(indexe)%ions(indexi - 1)%levels(ind_I)%photcros)) THEN
      npoints = SIZE(elements(indexe)%ions(indexi - 1)%levels(ind_I)%photcros(1,:))
     ELSE
      npoints = 0
     END IF
     IF(npoints /= 0) THEN
      ALLOCATE(crossfreq(npoints))
      crossfreq(1:npoints) = elements(indexe)%ions(indexi - 1)%levels(ind_I)%photcros(1, 1:npoints)
      ! frequency
      freq = (MINVAL(elements(indexe)%ions(indexi)%levels(:)%exci_energy) - &
              elements(indexe)%ions(indexi - 1)%levels(ind_I)%exci_energy) / const_h
      !write(*,*) 'i_colion: excienergy1 = ', &
      ! elements(indexe)%ions(indexi)%levels(act_level)%exci_energy, &
      ! ' excienergy2 =  ',  &
      ! elements(indexe)%ions(indexi - 1)%levels(I)%exci_energy, &
      ! 'freq = ', freq
      var_x = (const_h * freq) / (const_kB * temp)
      ! exponential integral function (calculation of eif)
      CALL exp_int_func(1, var_x, eif)
      actPoint = 0
      DO ind_J = 1, npoints
       IF(freq < crossfreq(ind_J)) THEN
        actPoint = ind_J
        EXIT
       END IF
      END DO
      ! if the frequency is out of range of the frequency interval
      ! the total rate will be equal to zero
      !write(*,*) 'i_colion: actPoint = ', actPoint
      IF(actPoint == 0 .OR. actPoint == 1) THEN
       actirates%Lma_int_reccol(ind_I) = 0.D0
       actirates%Lma_reccol(ind_I) = 0.D0
       DEALLOCATE(crossfreq)
       CYCLE
      END IF
      ! now we have to do a linear interpolation between the points actPoint - 1 and actPoint
      freq1 = elements(indexe)%ions(indexi - 1)%levels(ind_I)%photcros(1, actPoint - 1)
      freq2 = elements(indexe)%ions(indexi - 1)%levels(ind_I)%photcros(1, actPoint)
      func1 = elements(indexe)%ions(indexi - 1)%levels(ind_I)%photcros(2, actPoint - 1)
      func2 = elements(indexe)%ions(indexi - 1)%levels(ind_I)%photcros(2, actPoint)
      ali = (func1 - func2) / (freq1 - freq2)
      bli = (func2 * freq1 - func1 * freq2) / (freq1 - freq2)
      cross_sect = ali * freq + bli
      ! gindex
      IF(indexi == 1) THEN
        gindex = 0.1
      ELSE IF(indexi == 2) THEN
        gindex = 0.2
      ELSE IF(indexi > 2) THEN
        gindex = 0.3
      END IF
      ! populations calculation
      CALL populations(indexe, indexi - 1, ind_I, current_mgi, pop_number)
      exci_energy = elements(indexe)%ions(indexi - 1)%levels(ind_I)%exci_energy
      gr_exci_energy = MINVAL(elements(indexe)%ions(indexi)%levels(:)%exci_energy)
      stat_weight = elements(indexe)%ions(indexi - 1)%levels(ind_I)%stat_waight
      actirates%Lma_int_reccol(ind_I) = pop_number * el_dens * coll_const / temp**(1.0/2.0) * gindex * cross_sect * &
        exp(-var_x) / var_x * exci_energy 
      actirates%Lma_reccol(ind_I) = pop_number * el_dens * coll_const / temp**(1.0/2.0) * gindex * cross_sect * &
        exp(-var_x) / var_x * (gr_exci_energy - exci_energy) 
      !write(*,*) 'i_colion: el = ', indexe, ' ion = ', indexi, ' e - e0 = ', (exci_energy - gr_exci_energy)
      ! for now it will be equal to zero
      !actirates%Lma_reccol(ind_I) = 0.D0
      Zintrecomb = Zintrecomb + actirates%Lma_int_reccol(ind_I)
      Zrecomb = Zrecomb + actirates%Lma_reccol(ind_I)
      !print*, 'collion_rates: actirates%Lma_reccol = ', actirates%Lma_reccol(ind_I)
      DEALLOCATE(crossfreq)
     END IF
     !write(*,*) 'i_colion: Zintrecomb = ', Zintrecomb, ' Zrecomb = ', Zrecomb
    END DO
   ELSE
    Zintrecomb = 0.D0
    Zrecomb = 0.D0
   END IF
   ! print*, 'i_colion: Zion = ', Zion, ' Zrecomb = ', Zrecomb
 END SELECT

CASE DEFAULT
 STOP 'bound_free_rates: non valid approximation was chosen'
END SELECT

write(*,*) 'i_colion: END Zion = ', Zion, ' Zrecomb = ', Zrecomb, ' Zintrecomb = ', Zintrecomb

END SUBROUTINE i_colion
