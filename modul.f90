MODULE types

!#ifdef MPI_ON
!  include 'mpif.h'
!#endif


  IMPLICIT NONE

! Type Definitions

  TYPE grid_cell  
     ! won't be used in dynamical grid
     INTEGER                         :: model_index
     INTEGER, DIMENSION(3)           :: indexc
     DOUBLE PRECISION, DIMENSION(3)  :: corner
!     REAL                :: deltax, deltay, deltaz 
  END TYPE grid_cell

  TYPE dyn_grid_cell
      INTEGER                        :: model_index
      INTEGER                        :: cell_index, up_cell, down_cell
      INTEGER, DIMENSION(6)          :: neighbour
      DOUBLE PRECISION, DIMENSION(3) :: corner, width
  END TYPE dyn_grid_cell

  TYPE photon 
     INTEGER                         :: cell_numb, pack_numb, active
     INTEGER                         :: typ, next_cross, last_line
     DOUBLE PRECISION, DIMENSION(3)  :: pos, dir 
     DOUBLE PRECISION                :: e_cmf, e_rf, freq_cmf, freq_rf, delta_s
  END TYPE photon


  TYPE grid_ion_t
     DOUBLE PRECISION                :: gl_pop, tot_pop
  END TYPE grid_ion_t


  TYPE grid_comp_t
     DOUBLE PRECISION                :: abund 
     TYPE(grid_ion_t), ALLOCATABLE   :: grid_ion(:)
  END TYPE grid_comp_t


  TYPE modelgrid 
     INTEGER                         :: assoc_cells
     DOUBLE PRECISION                :: zwind, rwind
     DOUBLE PRECISION                :: T, J, rho, vel, e_dens
     TYPE(grid_comp_t), ALLOCATABLE  :: grid_comp(:)
  END TYPE modelgrid


  TYPE spec_type
     DOUBLE PRECISION                :: freq, flux
     INTEGER                         :: esc
  END TYPE spec_type


  TYPE line_list
     INTEGER                         :: indexe, indexi, lower, upper
     DOUBLE PRECISION                :: freq, A_ul, f_ul
  END TYPE line_list

  TYPE ion_levels 
     INTEGER                         :: nuptrans, ndowtrans
     INTEGER                         :: l_index
     DOUBLE PRECISION                :: exci_energy, stat_waight
     CHARACTER(LEN=3)                :: elconf
  END TYPE ion_levels

  TYPE element_ions 
     INTEGER                         :: nlevels, ion_stage
     DOUBLE PRECISION                :: ion_potential
     TYPE(ion_levels), ALLOCATABLE   :: levels(:)
  END TYPE element_ions

  TYPE atom_elements 
     INTEGER                         :: indexe, atom_number, nions
     DOUBLE PRECISION                :: atom_mass
     CHARACTER(20)                   :: levelfile, transitionfile
     TYPE(element_ions), ALLOCATABLE :: ions(:)
  END TYPE atom_elements

  TYPE virt_particle
     DOUBLE PRECISION, DIMENSION(3)  :: pos
     DOUBLE PRECISION                :: weight
  END TYPE virt_particle

! Global variables
  DOUBLE PRECISION                   :: xmax, ymax, zmax
  DOUBLE PRECISION, DIMENSION(3)     :: basic_cell_width
  INTEGER                            :: dyngrid
  DOUBLE PRECISION                   :: R_star, R_inf, V_inf, M_dot, T_eff
  DOUBLE PRECISION, ALLOCATABLE      :: incomingflux(:,:)
  INTEGER                            :: nx_cell, ny_cell, nz_cell, Ngrid, destroyed_pack
  INTEGER                            :: dummypackage, n_nubin, n_modelgrid

  TYPE(modelgrid), ALLOCATABLE       :: model_grid(:)
  TYPE(grid_cell), ALLOCATABLE       :: cell(:)   
  TYPE(dyn_grid_cell), ALLOCATABLE   :: dyn_cell(:)   
  TYPE(photon), ALLOCATABLE          :: package(:)

  TYPE(line_list), ALLOCATABLE       :: linelist(:)
  TYPE(atom_elements), ALLOCATABLE   :: elements(:)
  TYPE(virt_particle), ALLOCATABLE   :: virtual_particle(:)

  INTEGER                            :: idum
  INTEGER                            :: debug
! flux from existing input file
  INTEGER                            :: inputflux, inputmodel
!#ifdef MPI_ON
!  INTEGER                            :: n_tasks, my_rank
!#endif


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
  INTEGER, PARAMETER                 :: no_line = -99

!! Atomic data
  ! Total number of chemical elements in the simulation
  INTEGER                            :: n_elements
  ! Total number of line transitions in the simulation
  INTEGER                            :: ntransitions
  ! Specify backgraund model type (1-D, 2-D, 3-D)
  INTEGER                            :: model_type
  ! Specify minimal size of dynamic cell
  DOUBLE PRECISION, PARAMETER          :: minwidth = 1E8



!! Physical constants
  DOUBLE PRECISION, PARAMETER        :: pi=3.1415926535897932D+00,me_g=9.109534D-28,mp_g=1.6726485D-24, sigma_e=6.6516D-25,&
                                       h=6.626176D-27,light_speed=2.99792458D+10,e_charge=4.803242D-10,ftran=0.6407D+00, &       
                                       nio=4.5655967D+14,const=1.D-04,vel_ter=920.0D+05,r_sun=695990.D+05,beta=2.0D+00,  &   
                                       BOLK=1.380662D-16,m_sun=1.989D+33, sigma =5.6704D-05 !ergcm^(-2)s(-1)K(-4) !D. H. Cohen et al.2012
  DOUBLE PRECISION, PARAMETER        :: parsec=30.857D17, e_v = 1.60217646D-12, saha_const=2.0706839D-16, b = 1.D0

!! Parameters for testing
  DOUBLE PRECISION                   :: freq_line

!  DOUBLE PRECISION, PARAMETER       :: osc_line=0.416D0 ! for Ly_alph line
!  DOUBLE PRECISION, PARAMETER       :: nu_min= 2.D15, nu_max=3.D15  !nu_min= 1.D14, nu_max=1.D17,
!  DOUBLE PRECISION, PARAMETER       :: nu_min= 2.14286D15, nu_max=3.D15

! Define the min and max wavelenght range in cm for the synthetic spectrum calculation 1A = 1.D-8 cm
   DOUBLE PRECISION, PARAMETER        :: nu_min = 3.D14, nu_max = 3.7D15 ! in cm (800 - 10000 A)
  !DOUBLE PRECISION, PARAMETER       :: nu_min = 2.4D15, nu_max = 2.5D15 ! in cm (1150 - 1250 A)

END MODULE types
