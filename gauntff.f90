SUBROUTINE gauntff(freq, temp, gff)

USE types
USE constants
IMPLICIT NONE

! we calculate a FF Gaunt factor after Mihalas(1967) using the approximate
! formula (A2) which is valid for packets with wavelenghts from 100 A to
! 10000 A

! constants in the formula
DOUBLE PRECISION, PARAMETER                     :: const1 = 1.070192
DOUBLE PRECISION, PARAMETER                     :: const2 = 3.9999187E-3
DOUBLE PRECISION, PARAMETER                     :: const3 = -7.8622889E-5
DOUBLE PRECISION, PARAMETER                     :: const4 = 0.26061249
DOUBLE PRECISION, PARAMETER                     :: const5 = 6.4628601E-2
DOUBLE PRECISION, PARAMETER                     :: const6 = -6.1953813E-4
DOUBLE PRECISION, PARAMETER                     :: const7 = -0.57917786
DOUBLE PRECISION, PARAMETER                     :: const8 = -3.7542343E-2
DOUBLE PRECISION, PARAMETER                     :: const9 = -1.3983474E-5
DOUBLE PRECISION, PARAMETER                     :: const10 = 0.34169006
DOUBLE PRECISION, PARAMETER                     :: const11 = 1.1852264E-2

! theta constant
DOUBLE PRECISION, PARAMETER                     :: thetaConst = 5040.D0
! lambda const (we want lambda in microns)
DOUBLE PRECISION, PARAMETER                     :: laConst = 1.D6

! input variables
DOUBLE PRECISION                                :: freq, temp
! output variables
DOUBLE PRECISION                                :: gff
!
DOUBLE PRECISION                                :: theta, lam, x
DOUBLE PRECISION                                :: light_speed_m = 1.D-2 * light_speed

theta = thetaConst / temp
lam = light_speed_m / freq
x = 1.0 / (laConst * lam)

gff = (const1 + const2 / theta + const3 / theta**2) + &
          (const4 + const5 / theta + const6 / theta**2) / x + &
          (const7 + const8 / theta + const9 / theta**2) /x /x + &
          (const10 + const11 / theta) /x /x /x

! write(*,*) 'gauntff: gff = ', gff

END SUBROUTINE gauntff
