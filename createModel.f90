!_______________________________________________________________________________________________
! this program creates a data file containg an analytical physical model defined with parameters
! by a user
! output file: 'synthetic_model.dat'
!_______________________________________________________________________________________________
PROGRAM createModel

IMPLICIT NONE

INTEGER                                         :: clump_type
DOUBLE PRECISION                                :: Rstar, Rinf, Teff, Vinf
DOUBLE PRECISION                                :: R1, R2, stred1, stred2
DOUBLE PRECISION                                :: delta, deltacl
DOUBLE PRECISION                                :: rho0, rhocl

DOUBLE PRECISION                                :: sigma
DOUBLE PRECISION                                :: rdist
DOUBLE PRECISION                                :: actrho, actvel
DOUBLE PRECISION                                :: rho, rhoclump
LOGICAL                                         :: clump=.false.
LOGICAL                                         :: turnclumpoff=.true.

DOUBLE PRECISION                                :: velocity
INTEGER                                         :: vel_approx
INTEGER                                         :: model_type

INTEGER                                         :: I, J, K
INTEGER                                         :: Nx, Ny, Nz
DOUBLE PRECISION                                :: xsour, ysour, zsour
DOUBLE PRECISION                                :: w_x, w_y, w_z
DOUBLE PRECISION                                :: xmax, ymax, zmax
DOUBLE PRECISION                                :: radial, cur_vel
DOUBLE PRECISION, DIMENSION(3)                  :: cur_vel_vec

DOUBLE PRECISION                                :: temperature
DOUBLE PRECISION                                :: Tstar, Tinf
DOUBLE PRECISION                                :: cur_temp

CHARACTER (LEN=200)                             :: output_filename

REAL                                            :: beta

CALL get_command_argument(1, output_filename)

IF(LEN_TRIM(output_filename) == 0) THEN
 output_filename = 'synthetic_model.dat'
END IF
write(*,*) 'createModel: output_filename = ', output_filename

beta = 1.3
! model_type
! 1 -- spherically symmetric model
! 3 -- full 3D model
model_type = 3
! vel_approx
! -1 -- zero velocity
! 0 -- homologous approximation
! 1 -- beta law
! 2 -- sin velocity
! 10 -- inverse homologous approximation
vel_approx = 1
clump_type = 1
! 3D grid informations
Nx = 15
Ny = 15
Nz = 15
! lower boundary condition
Teff = 1.473441441281968036e+04
Rstar = 5.616000000000000000e+12
! every distance is written as a multiply of R_star
R1 = 6.D0
R2 = 6.5D0
sigma = 0.7D0
stred1 = 6.05D0
stred2 = 6.45D0

rho0 = 1.D-8
rhocl = 1.D-12

delta = 0.02
deltacl = 0.01

Rinf = 20.0 * Rstar
rdist = 1.0
Vinf = 3000.D+5


xmax = Rinf/Rstar + 1.0/float(10*Nx)
ymax = Rinf/Rstar + 1.0/float(10*Ny) 
zmax = Rinf/Rstar + 1.0/float(10*Nz) 

Tstar = 90000
Tinf = 30000


I = 0

OPEN(1, FILE=output_filename)

write(1, *) Teff
write(1, *) Rstar
write(1, *) Rinf
write(1, *) Vinf

IF(model_type == 3) THEN
 ! write(1,*) xmax, ymax, zmax
 write(1,*) Nx, Ny, Nz
END IF

SELECT CASE(model_type)

CASE(1)
 DO 
  ! v clumpu
  IF(rdist >= R1 .AND. rdist < R2) THEN
   IF(turnclumpoff) THEN
    rdist = rdist + delta
   ELSE
    rdist = rdist + deltacl
   END IF
   clump=.true.
  ! mimo clump
  ELSE IF (rdist < R1 .OR. (rdist >= R2 .AND. rdist < Rinf/Rstar)) THEN
   rdist = rdist + delta
   clump=.false.
  ! mimo model
  ELSE IF (rdist >= Rinf/Rstar) THEN
   EXIT
  END IF
  actrho = rho(rdist, R1, R2, sigma, stred1, stred2, rho0, rhocl, clump, &
   turnclumpoff, clump_type)
  actvel = velocity(vel_approx, rdist, Rstar, Rinf, Vinf, beta)
  I = I + 1
  write(1, *) I, rdist , actvel , actrho, 15500
 END DO
 
 CLOSE(1)
CASE(3)
 w_x = 2.0 * xmax / DBLE(Nx - 1)
 w_y = 2.0 * ymax / DBLE(Ny - 1)
 w_z = 2.0 * zmax / DBLE(Nz - 1)

 DO I = 1, Nx
  DO J = 1, Ny
   DO K = 1, Nz
    xsour = -xmax + DBLE((I - 1)) * w_x
    ysour = -ymax + DBLE((J - 1)) * w_y
    zsour = -zmax + DBLE((K - 1)) * w_z
    radial = sqrt(xsour**2 + ysour**2 + zsour**2)

    cur_vel = velocity(vel_approx, radial, Rstar, Rinf, Vinf, beta)
    write(*,*) 'createModel: velocity = ', cur_vel
    IF(radial /= 0) THEN
     cur_vel_vec = cur_vel * (/xsour, ysour, zsour/)/radial
    ELSE
     cur_vel_vec = (/ 0.0, 0.0, 0.0 /)
    END IF

    cur_temp = temperature(0, radial, Tstar, Tinf, Rstar, Rinf)

    actrho = rho(radial, R1, R2, sigma, stred1, stred2, rho0, rhocl, clump, &
     turnclumpoff, clump_type)
    write(1,*) xsour, ysour, zsour, cur_vel_vec, actrho, 15500
   END DO
  END DO
 END DO


CASE DEFAULT
END SELECT
CLOSE(1)

END PROGRAM createModel

FUNCTION rho(r, R1, R2, sigma, stred1, stred2, rho0, rhocl, inClump, &
 turnclumpoff, clump_type)
DOUBLE PRECISION                :: rho
DOUBLE PRECISION                :: r
DOUBLE PRECISION                                :: R1, R2, stred1, stred2, sigma
DOUBLE PRECISION                                :: rho0, rhocl
LOGICAL                                         :: inClump, turnclumpoff
INTEGER                                         :: clump_type

SELECT CASE(clump_type)
 CASE(1)
  rhoclump = rhocl * exp(-(r - (R2 + R1)/2.0)**2/(R2 - R1))
 CASE(2)
  rhoclump = rhocl * exp(-(r - stred1)**2/sigma**2) + rhocl * exp(-(r - stred2)**2/sigma**2)
 CASE DEFAULT
  write(*,*) 'wrong clump type'
  STOP
END SELECT

IF (inClump .AND. .NOT. turnclumpoff) THEN
 rho = rho0 * r**(-2) + rhoclump
ELSE
 rho = rho0 * r**(-2)
END IF

IF(r == 0) rho = 0e0


END FUNCTION rho


FUNCTION velocity(approx, r, Rstar, Rinf, Vinf, beta)

INTEGER                         :: approx
DOUBLE PRECISION                :: velocity, r, v0
DOUBLE PRECISION                :: Rinf, Vinf, Rstar

DOUBLE PRECISION                :: aindex, bindex

SELECT CASE(approx)

! homologous approximation
CASE(-1)
 velocity = 0.0
CASE(0)
 IF(r*Rstar>Rinf) THEN
  velocity = Vinf
 ELSE IF(r*Rstar<Rstar) THEN
  velocity = 0.D0
 ELSE
  velocity = r * Vinf * Rstar/Rinf
 END IF
 RETURN
CASE(1)
 velocity = Vinf * (1.0 - 1.0/r)**beta
 write(*,*) 'Vinf = ', Vinf, ' Rstar = ', Rstar, ' r = ', r
CASE(2)
 velocity = 0.4 * Vinf * sin(4.0 * r*Rstar/(Rinf - Rstar)) + Vinf*0.5
! inv homologous approximation
CASE(10)
 v0 = Rstar/Rinf * Vinf
 aindex = - (Vinf - v0)/(Rinf - Rstar)
 bindex = (Vinf * Rinf - v0 * Rstar)/(Rinf - Rstar)
 velocity = aindex * r * Rstar + bindex
 RETURN
END SELECT


END FUNCTION


FUNCTION temperature(approx, r, Tstar, Tinf, Rstar, Rinf)

 INTEGER                                :: approx
 DOUBLE PRECISION                       :: r
 DOUBLE PRECISION                       :: temperature
 DOUBLE PRECISION                       :: Teff, Tstar, Rstar, Rinf, Tinf
 DOUBLE PRECISION                       :: a_ind, b_ind

 SELECT CASE(approx)
 CASE(0)
  temperature = 50000
 CASE(1)
  write(*,*) 'Rstar = ', Rstar, ' Rinf = ', Rinf

  a_ind = (Tinf - Tstar)/(Rinf/Rstar - 1)
  b_ind = (Tstar * Rinf/Rstar - Tinf)/(Rinf/Rstar - 1)
  write(*,*) 'a_ind = ', a_ind, ' b_ind = ', b_ind

  temperature = a_ind * r + b_ind
  write(*,*) 'temperature = ', temperature, ' r = ', r
  if(temperature < 0.0) temperature = Tinf
 END SELECT

END FUNCTION 
