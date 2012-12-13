PROGRAM create_model_data

!  USE types

  IMPLICIT NONE 
   
  INTEGER                                   :: I, J
  INTEGER, PARAMETER                        :: n_modelgrid=200, n_elements=26, pi=3.1415926535897932D+00
  DOUBLE PRECISION, PARAMETER               :: m_sun=1.989D+33, r_sun=695990.D+05
  DOUBLE PRECISION                          :: delta_r, r, R_inf, R_star, T_eff, V_inf, M_dot, R_star_solarunit
  DOUBLE PRECISION, DIMENSION(n_modelgrid)  :: rwind, rho, vel, T
  DOUBLE PRECISION, DIMENSION(n_elements)   :: mass_frac
  CHARACTER(LEN=80)                         :: String


  OPEN (UNIT=2, FILE='model_data.dat')
 
  String = "(1I4, 1F10.5, 3E14.5, 26F6.2)"

  T_eff = 30000.D0  
  WRITE(2, *) T_eff
  R_star_solarunit = 9.9D0    ! in R_sun
  WRITE(2, *) R_star_solarunit
  R_star = 1.D0    ! in R_star
  R_inf  = 10.D0   ! in R_star
!  WRITE(2, *) R_inf
  V_inf = 3000.D0  ! km/s
!  WRITE(2, *) V_inf
  M_dot = 2.2D-7   ! m_sun/s
!  WRITE(2, *) M_dot
  WRITE(2, *) n_modelgrid
  delta_r = (R_inf - R_star)/n_modelgrid
   
  DO I = 1, n_modelgrid
     r = R_star + I * delta_r
     rwind(I) = r
     vel(I)   = (r / R_inf) * V_inf
     rho(I)   = (M_dot * m_sun/(3600.D0*24.D0*365.25D0)) / (4.D0 * pi * (rwind(I) * R_star_solarunit * r_sun)**2 * vel(I) * 1.D5)     
     T(I)     = T_eff * 3.D0 / 4.D0
     DO J = 1, n_elements
        IF (J .EQ. 1) THEN
            mass_frac(J) = 1.D0
        ELSE
            mass_frac(J) = 0.D0
        END IF
     END DO
!     massfraction(I)%H  = 1.D0!7.374E-01
!     massfraction(I)%He = 0.  !2.492E-01
!     massfraction(I)%Li = 0.  !5.698E-11
!     massfraction(I)%Be = 0.  !1.582E-10 
!     massfraction(I)%B  = 0.  !3.964E-09 
!     massfraction(I)%C  = 0.  !2.365E-03
!     massfraction(I)%N  = 0.  !6.928E-04 
!     massfraction(I)%O  = 0.  !5.733E-03
!     massfraction(I)%F  = 0.  !5.046E-07 
!     massfraction(I)%Ne = 0.  !1.257E-03
!     massfraction(I)%Na = 0.  !2.923E-05
!     massfraction(I)%Mg = 0.  !7.079E-04 
!     massfraction(I)%Al = 0.  !5.563E-05 
!     massfraction(I)%Si = 0.  !6.649E-04 
!     massfraction(I)%P  = 0.  !5.825E-06 
!     massfraction(I)%S  = 0.  !3.092E-04
!     massfraction(I)%Cl = 0.  !8.202E-06
!     massfraction(I)%Ar = 0.  !7.341E-05
!     massfraction(I)%K  = 0.  !3.065E-06
!     massfraction(I)%Ca = 0.  !6.415E-05
!     massfraction(I)%Sc = 0.  !4.646E-08
!     massfraction(I)%Ti = 0.  !3.121E-06
!     massfraction(I)%V  = 0.  !3.172E-07 
!     massfraction(I)%Cr = 0.  !1.660E-05 
!     massfraction(I)%Mn = 0.  !1.082E-05
!     massfraction(I)%Fe = 0.  !1.292E-03
     WRITE(2, String) I, rwind(I), vel(I), rho(I), T(I), mass_frac
!                massfraction(I)%H,  massfraction(I)%He, &
!                massfraction(I)%Li, massfraction(I)%Be, massfraction(I)%B,  massfraction(I)%C,  massfraction(I)%N, &
!                massfraction(I)%O,  massfraction(I)%F,  massfraction(I)%Ne, massfraction(I)%Na, massfraction(I)%Mg,&
!                massfraction(I)%Al, massfraction(I)%Si, massfraction(I)%P,  massfraction(I)%S,  massfraction(I)%Cl,&
!                massfraction(I)%Ar, massfraction(I)%K,  massfraction(I)%Ca, massfraction(I)%Sc, massfraction(I)%Ti,&
!                massfraction(I)%V,  massfraction(I)%Cr, massfraction(I)%Mn, massfraction(I)%Fe
  END DO

END PROGRAM create_model_data
