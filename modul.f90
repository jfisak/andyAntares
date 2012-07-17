MODULE types

IMPLICIT NONE

! Type Definitions

  TYPE grid_cell  
     INTEGER                         :: model_index
     INTEGER, DIMENSION(3)           :: indexc
     DOUBLE PRECISION, DIMENSION(3)  :: corner
!     REAL                :: deltax, deltay, deltaz 
  END TYPE grid_cell

  TYPE photon 
     INTEGER                         :: cell_numb, pack_numb, active, typ, last_cross
     DOUBLE PRECISION, DIMENSION(3)  :: pos, dir 
     DOUBLE PRECISION                :: e_cmf, e_rf, freq_cmf, freq_rf
  END TYPE photon



! Global variables
  DOUBLE PRECISION                   :: xmax, ymax, zmax, cell_width
  INTEGER                            :: nx_cell, ny_cell, nz_cell, Ngrid

  TYPE(grid_cell), ALLOCATABLE       :: cell(:) 
    
  TYPE(photon), ALLOCATABLE          :: package(:)

  INTEGER                            :: idum
  INTEGER                            :: debug



! Globally defined numerical constants 
!! Different packet types
  INTEGER, PARAMETER                 :: type_escaped=-99 
  INTEGER, PARAMETER                 :: type_rpkt=0 
  INTEGER, PARAMETER                 :: type_kpkt=1

!! Packet's memory of which cell surface it crossed last
  INTEGER, PARAMETER                 :: posx=1
  INTEGER, PARAMETER                 :: negx=2
  INTEGER, PARAMETER                 :: posy=3
  INTEGER, PARAMETER                 :: negy=4
  INTEGER, PARAMETER                 :: posz=5
  INTEGER, PARAMETER                 :: negz=6
  INTEGER, PARAMETER                 :: NONE = -99


!! Physical constants
  DOUBLE PRECISION, PARAMETER       :: pi=3.1415926535897932D+00,me_g=9.109534D-28,mp_g=1.6726485D-24,sigma_e=6.6516D-25,&
                                       h=6.626176D-27,light_speed=2.99792458D+10,e_charge=4.803242D-10,ftran=0.6407D+00, &       
                                       nio=4.5655967D+14,const=1.D-04,vel_ter=920.0D+05,r_sun=695990.D+05,beta=2.0D+00,  &   
                                       BOLK=1.380662D-16,m_sun=1.989D+33,  T=39000, sigma=5.6704D+05 !ergcm^(-2)s(-1)K(-4) !D. H. Cohen et al.2012

END MODULE types
