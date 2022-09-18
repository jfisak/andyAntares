PROGRAM skipDensity

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
LOGICAL                                         :: turnclumpoff=.false.

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

model_type = 3
vel_approx = 0
clump_type = 1
! 3D grid informations
Nx = 5
Ny = 5
Nz = 5
xmax = 11
ymax = 11
zmax = 11
! star info
Teff = 1.473441441281968036e+04
Rstar = 5.616000000000000000e+14
! every distance is written as a multiply of R_star
R1 = 6.D0
R2 = 6.5D0
sigma = 0.7D0
stred1 = 6.05D0
stred2 = 6.45D0

rho0 = 1.D-13
rhocl = 1.D-12

delta = 0.25
deltacl = 0.01

Rinf = 10.0 * Rstar
rdist = 1.0
Vinf = 30000.D+5


I = 0

OPEN(1, FILE='skip_model.dat')

write(1, *) Teff
write(1, *) Rstar
write(1, *) Rinf
write(1, *) Vinf

IF(model_type == 3) THEN
 write(1,*) xmax, ymax, zmax
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
  actvel = velocity(vel_approx, rdist, Rstar, Rinf, Vinf)
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

    cur_vel = velocity(vel_approx, radial, Rstar, Rinf, Vinf)
    IF(radial /= 0) THEN
     cur_vel_vec = cur_vel * (/xsour, ysour, zsour/)/radial
    ELSE
     cur_vel_vec = (/ 0.0, 0.0, 0.0 /)
    END IF


    actrho = rho(radial, R1, R2, sigma, stred1, stred2, rho0, rhocl, clump, &
     turnclumpoff, clump_type)
    write(*,*) cur_vel_vec
    write(1,*) xsour, ysour, zsour, cur_vel_vec, actrho, 15500
   END DO
  END DO
 END DO


CASE DEFAULT
END SELECT

END PROGRAM skipDensity

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


END FUNCTION rho


FUNCTION velocity(approx, r, Rstar, Rinf, Vinf)

INTEGER                         :: approx
DOUBLE PRECISION                :: velocity, r, v0
DOUBLE PRECISION                :: Rinf, Vinf, Rstar

DOUBLE PRECISION                :: aindex, bindex

SELECT CASE(approx)

! homologous approximation
CASE(0)
 velocity = r*Rstar/Rinf * Vinf
 RETURN
! inv homologous approximation
CASE(1)
 v0 = Rstar/Rinf * Vinf
 aindex = - (Vinf - v0)/(Rinf - Rstar)
 bindex = (Vinf * Rinf - v0 * Rstar)/(Rinf - Rstar)
 velocity = aindex * r * Rstar + bindex
 RETURN
CASE(2)
 velocity = 0.4 * Vinf * sin(4.0 * r*Rstar/(Rinf - Rstar)) + Vinf*0.5
END SELECT


END FUNCTION
