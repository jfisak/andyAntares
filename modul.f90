MODULE types

#ifdef MPI_ON
  include 'mpif.h'
#endif


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

  TYPE modelgrid 
     DOUBLE PRECISION                :: rho, vel, rwind
  END TYPE modelgrid

  TYPE spec_type
     DOUBLE PRECISION                :: freq, flux
     INTEGER                         :: esc
  END TYPE spec_type



! Global variables
  DOUBLE PRECISION                   :: xmax, ymax, zmax, cell_width
  DOUBLE PRECISION                   :: R_star, R_inf, V_inf, M_dot, T_eff
  INTEGER                            :: nx_cell, ny_cell, nz_cell, Ngrid, n_modelgrid, n_nubin
  INTEGER                            :: dummypackage

  TYPE(modelgrid), ALLOCATABLE       :: model_grid(:)
  TYPE(grid_cell), ALLOCATABLE       :: cell(:)   
  TYPE(photon), ALLOCATABLE          :: package(:)

  INTEGER                            :: idum
  INTEGER                            :: debug

#ifdef MPI_ON
  INTEGER                            :: n_tasks, my_rank
#endif


! Globally defined numerical constants 
!! Different packet types
  INTEGER, PARAMETER                 :: type_escaped=-99 
  INTEGER, PARAMETER                 :: type_rpkt=0 
  INTEGER, PARAMETER                 :: type_kpkt=1
  INTEGER, PARAMETER                 :: type_ipkt=2

  INTEGER, PARAMETER                 :: rpkt_eventtype_changecell=1
  INTEGER, PARAMETER                 :: rpkt_eventtype_lineinteraction=2
  INTEGER, PARAMETER                 :: rpkt_eventtype_continuum=3

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
                                       BOLK=1.380662D-16,m_sun=1.989D+33, sigma=5.6704D-05 !ergcm^(-2)s(-1)K(-4) !D. H. Cohen et al.2012
  DOUBLE PRECISION, PARAMETER       :: parsec=30.857D17

!! Parameters for testing
  DOUBLE PRECISION, PARAMETER        :: freq_line=light_speed/1216.D-8, osc_line=0.416D0 ! for Ly_alph line
!  DOUBLE PRECISION, PARAMETER        :: nu_min= 2.D15, nu_max=3.D15  !nu_min= 1.D14, nu_max=1.D17,
  DOUBLE PRECISION, PARAMETER         :: nu_min= 2.14286D15, nu_max=3.D15
!  DOUBLE PRECISION, PARAMETER        :: nu_min= 2.4D15, nu_max=2.542D15 


END MODULE types
