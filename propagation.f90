SUBROUTINE propagation(n_pack, opa_cell, lower_opa, delta_opa)
 
! Propagation of the photon in 3D grid

  USE types

  IMPLICIT NONE    

    INTEGER                           :: I, I_esc, pack_index, nc, next_cell, n_pack
    DOUBLE PRECISION                  :: tau, xi, ran2, tau_rand, cell_dist, event_dist, r,   &
                                         rho_cell, I_beta, opa_cell, lower_opa, delta_opa
    DOUBLE PRECISION, PARAMETER       :: rho=1.D0

!    DOUBLE PRECISION, PARAMETER       :: opa_cell=1.D0/20.D0, rho=1.D0
!    DOUBLE PRECISION, PARAMETER       :: opa_cell=5.D-3, rho=5.D-2

  OPEN (UNIT=3, FILE='position.dat')  
  OPEN (UNIT=4, FILE='escprob-photsphere-sphere-abs-10000-nn.dat')

! Initial value of the escaped photons from the sphere which can be fit within the grid 
! (it corresponds to the emergent intensity)
  I_esc = 0

! Loop over the photon packages   
  DO pack_index=1,n_pack
    ! IF (pack_index .EQ. 164) THEN 
    !     debug = 2
!     ELSE 
!         debug = 0 
!     END IF

!     IF ((pack_index .EQ. 108) .AND. (opa_cell .GT.  lower_opa)) THEN
!	debug =3 
!     ENDIF

!    Initial value for optical depth (which will be increment)
     tau=0.D0
!    Random choice of optical depth the photon will travel
10   xi = ran2(idum)     
     IF (xi .EQ. 0.D0) GOTO 10    
     tau_rand = -LOG(xi)

!    Number of the current cell
     nc = package(pack_index)%cell_numb
     
!     WRITE(3,*) pack_index, package(pack_index)%cell_numb, package(pack_index)%pos
     I=0
!    Loop as long as the photon is active 
     DO  WHILE (package(pack_index)%active .EQ. 1)
       IF (debug .EQ. 1) THEN 
           print*, pack_index, I, nc, tau, tau_rand, package(pack_index)%active, package(pack_index)%typ
       END IF


!      Radius of the sphere which can be fit within the grid 
       r =SQRT(((cell(package(pack_index)%cell_numb)%corner(1) + cell_width/2.D0))**2 &
       + ((cell(package(pack_index)%cell_numb)%corner(2) + cell_width/2.D0))**2       &
       + ((cell(package(pack_index)%cell_numb)%corner(3) + cell_width/2.D0))**2)
!       print*,  r
!      Density out of the sphere is zero     
       IF (r .GT. xmax) THEN
          rho_cell = 0.D0
       ELSE 
           rho_cell = rho
       END IF
!       WRITE(3,*) pack_index, package(pack_index)%cell_numb, package(pack_index)%pos
!
!      Calculate the distance to cell surfaces and return the shortest distance and the cell 
!      the photon will go 
       IF (debug .EQ. 3) THEN
         print*, pack_index, I 
       ENDIF
       CALL boundary(pack_index, cell_dist, next_cell)
       IF (cell_dist .LT. 0.D0) STOP 'cell_dis < 0'
       IF (debug .EQ. 1) THEN 
           print*, package(pack_index)%dir
       END IF

!      Calculate the distance to the point of event for simplifying case without lines
       event_dist = (tau_rand-tau)/(opa_cell * rho)

      CALL event_dist()




       IF (debug .EQ. 1) THEN 
          print*, cell_dist, next_cell, event_dist , cell(package(pack_index)%cell_numb)%indexc
       END IF
!
       IF (debug .EQ. 3) THEN
!         print*, package(pack_index)%pos, package(pack_index)%dir
       ENDIF
       IF (event_dist .LT. cell_dist) THEN   
!         Move photon package from the curent position for some distance
          CALL move_package(pack_index, event_dist)
          CALL do_event(pack_index, tau, tau_rand)
	  IF (debug .EQ. 1) THEN 
              print*, 'do event', opa_cell * rho_cell * cell_dist
          END IF
       ELSE     
!         Calculate(accumulate) the optical depth along the package inside the cell (cell distance)
          tau = tau + opa_cell * rho_cell * cell_dist
!         Move package from the curent position for the cell_dist
          CALL move_package(pack_index, cell_dist)
!         If package escaped the calculation volume (next_cell=-99) then it become no-active and 
!         package type is update to the type_escaped, else the cell number is updated
          CALL change_cell(pack_index, next_cell)
          IF (debug .EQ. 1) THEN 
              print*, 'propagate ', opa_cell * rho_cell * cell_dist
          END IF
       END IF
       IF (debug .EQ. 3) THEN
         print*, cell_dist, event_dist
         print*, package(pack_index)%pos, package(pack_index)%dir, package(pack_index)%active
       ENDIF
!      Increase the 
       I=I+1
!         IF (I .EQ. 500) stop
     END DO
       
!     write(4,*) 'packet dead', package(pack_index)%typ, pack_index
!     print*, pack_index, nc, tau, tau_rand, package(pack_index)%active, package(pack_index)%typ
!
!    Count the escaped packages (all packages with type type_escaped)
     IF (package(pack_index)%typ .EQ. type_escaped) I_esc = I_esc +1
  END DO

! Calculate the ratio between number of escaped (I_esc) and sent (n_pack) packages - this corresponds 
! to the e^(-tau). 
  I_beta = float(I_esc)/float(n_pack)
  WRITE(4,*) opa_cell*(xmax-R_star), I_esc, n_pack, I_beta
!  print*, opa_cell, xmax, R_star
!  WRITE(4,*) opa_cell*xmax, I_esc, n_pack, I_beta

END 
