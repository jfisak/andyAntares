MODULE types

IMPLICIT NONE
SAVE
#if mpi==1
 include 'mpif.h'
#endif


 INTEGER, PARAMETER                 :: file_length = 180

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Type Definitions

  TYPE dyn_grid_cell
      INTEGER                        :: model_index=0, n_virt=0
      INTEGER                        :: up_cell, down_cell
      INTEGER, DIMENSION(6)          :: neighbor
      DOUBLE PRECISION, DIMENSION(3) :: corner, width
      DOUBLE PRECISION, DIMENSION(3) :: vec_vel
      INTEGER, DIMENSION(3)          :: n_sbgr
  END TYPE dyn_grid_cell

  TYPE photon 
     INTEGER                         :: cell_numb, active
     DOUBLE PRECISION                :: e_cmf, e_rf, freq_cmf, freq_rf, delta_s
     INTEGER                         :: typ, next_cross, last_line
     INTEGER                         :: n_interactions
     DOUBLE PRECISION, DIMENSION(3)  :: pos, dir 
     INTEGER                         :: l_ele, l_ion, l_lev, n_int = 0
     LOGICAL                         :: redShift, virtual
  END TYPE photon

  TYPE virtual_packet
   DOUBLE PRECISION, DIMENSION(3)       :: pos, dir
   INTEGER                              :: active, cell_index
   DOUBLE PRECISION                     :: freq_rf
  END TYPE virtual_packet


  TYPE grid_ion_t
     DOUBLE PRECISION                :: gl_pop, tot_pop
  END TYPE grid_ion_t


  TYPE grid_comp_t
     DOUBLE PRECISION                :: abund 
     TYPE(grid_ion_t), ALLOCATABLE   :: grid_ion(:)
  END TYPE grid_comp_t


  TYPE modelgrid 
     INTEGER                         :: assoc_cells
     DOUBLE PRECISION                :: width
     DOUBLE PRECISION                :: volume = 0.D0
     DOUBLE PRECISION                :: diff_param
     DOUBLE PRECISION                :: T = 0.D0, J = 0.D0, rho = 0.D0, vel = 0.D0, rwind, e_dens = 0.D0
     DOUBLE PRECISION                :: zwind, velang, angle
     DOUBLE PRECISION, DIMENSION(3)  :: vec_vel, vec_pos
     TYPE(grid_comp_t), ALLOCATABLE  :: grid_comp(:)
     ! is the diffusion approximation recommended?
     LOGICAL                         :: is_difapp = .false.
  END TYPE modelgrid


  TYPE spec_type
     DOUBLE PRECISION                :: freq, flux
     INTEGER                         :: esc
  END TYPE spec_type


  TYPE line_list
     INTEGER                         :: indexe, indexi, lower, upper
     DOUBLE PRECISION                :: freq, A_ul, f_lu
     INTEGER                         :: n_deexc, n_exc
     INTEGER(KIND=4)                 :: n_int
     LOGICAL                         :: counted
  END TYPE line_list

  TYPE ion_levels 
     INTEGER, ALLOCATABLE    :: linetransitions(:), lineuptransitions(:)
     DOUBLE PRECISION                :: exci_energy, stat_waight
     CHARACTER(LEN=30)               :: elconf
     LOGICAL                         :: phcrossform
     DOUBLE PRECISION, ALLOCATABLE   :: photcros(:,:), phcrosscoeff(:)
     DOUBLE PRECISION, ALLOCATABLE   :: population(:)
     DOUBLE PRECISION                :: phfreq
     INTEGER                         :: phfreqi
     INTEGER                         :: levelindex
     INTEGER                         :: vsplit
  END TYPE ion_levels

  TYPE element_ions 
     INTEGER                         :: ion_stage
     DOUBLE PRECISION                :: ion_potential
     TYPE(ion_levels), ALLOCATABLE   :: levels(:)
  END TYPE element_ions

  TYPE atom_elements 
     INTEGER                         :: indexe, atom_number, nions
     DOUBLE PRECISION                :: atom_mass
     DOUBLE PRECISION                :: abundance
     CHARACTER(LEN=file_length)                   :: levelfile='', transitionfile=''
     TYPE(element_ions), ALLOCATABLE :: ions(:)
  END TYPE atom_elements

  TYPE virt_point
     DOUBLE PRECISION, DIMENSION(3)  :: pos
     DOUBLE PRECISION                :: weight
     INTEGER                         :: ind_pcell, ind_mcell
  END TYPE virt_point

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Global variables
! propagation grid properties
  DOUBLE PRECISION                   :: xmax, ymax, zmax, xmin, ymin, zmin
  DOUBLE PRECISION, DIMENSION(3)     :: basic_cell_width
  INTEGER                            :: nx_cell, ny_cell, nz_cell, Ngrid, destroyed_pack
  INTEGER                            :: dyngrid
! information about saved propmod grid
  INTEGER                            :: saved_grid
! NLTE
  INTEGER                            :: nlte
  INTEGER                            :: abs_surface
! properties of a central star
  DOUBLE PRECISION                   :: R_star, R_inf, V_inf, V_0, M_dot, T_eff
  DOUBLE PRECISION                   :: Z_inf
! lower boundary condition
  DOUBLE PRECISION, ALLOCATABLE      :: incomingflux(:,:)
  ! number of points in a spectrum, nof points in modGrid, nop in propGrid
  INTEGER                            :: n_nubin, n_modelgrid, n_propgcells
  ! additional model grid variables
  INTEGER                            :: add_mg
  ! number of virtual point
  INTEGER                            :: Nvirtpoint
  ! velocity approximation
  INTEGER                            :: velApprox
  ! brtm activation
  LOGICAL                            :: calc_brtm

! fields for the given types
  TYPE(modelgrid), ALLOCATABLE       :: model_grid(:)
  TYPE(dyn_grid_cell), ALLOCATABLE, SAVE   :: dyn_cell(:)   
  TYPE(photon), ALLOCATABLE          :: package(:)
  TYPE(virtual_packet), ALLOCATABLE  :: vpackage(:)

  TYPE(line_list), ALLOCATABLE       :: linelist(:)
  TYPE(atom_elements), ALLOCATABLE   :: elements(:)
  TYPE(virt_point), ALLOCATABLE   :: virtual_point(:)
! variable for random number generation
  INTEGER                            :: idum
! is electron density values stored?
  INTEGER                            :: eldensfile
  CHARACTER(180)                      :: inputpopfile
! debug mode
  INTEGER                            :: debug
! flux from existing input file
  INTEGER                            :: inputflux, inputmodel
  CHARACTER(160)                         :: inputmodelFile, inputcomposition
! number of photoionization cross sections
  INTEGER                               :: n_photcrossect, n_tot_cont, n_ff = 0
! number of dummy packages
  INTEGER                            :: n_dummy_packs
! is the random seed initialized?
  LOGICAL, ALLOCATABLE               :: initrs(:)


! Globally defined numerical constants 
!! Different packet types
  INTEGER, PARAMETER                 :: type_escaped=-99 
  INTEGER, PARAMETER                 :: type_photosphere=99
  INTEGER, PARAMETER                 :: type_rpkt=0 
  INTEGER, PARAMETER                 :: type_kpkt=1
  INTEGER, PARAMETER                 :: type_ipkt=2
  INTEGER, PARAMETER                 :: type_dpkt=3
  INTEGER, PARAMETER                 :: type_vrpkt=4
  INTEGER, PARAMETER                 :: type_vkpkt=5
  INTEGER, PARAMETER                 :: type_vipkt=6
  INTEGER, PARAMETER                 :: type_vdpkt=7

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
  INTEGER, PARAMETER                 :: edxy=7
  INTEGER, PARAMETER                 :: edxz=8
  INTEGER, PARAMETER                 :: edyz=9

  INTEGER, PARAMETER                 :: ind_x = 1
  INTEGER, PARAMETER                 :: ind_y = 2
  INTEGER, PARAMETER                 :: ind_z = 3

  INTEGER, PARAMETER                 :: NONE = -99
  INTEGER, PARAMETER                 :: no_line = -99

  INTEGER                            :: my_rank
  INTEGER                            :: ierr
  INTEGER                            :: n_tasks
  CHARACTER(160)                      :: outputfolder=''
  CHARACTER(160)                      :: outputfile
!! Atomic data
 ! Total number of chemical elements in the simulation
  INTEGER                            :: n_elements
  ! Total number of line transitions in the simulation
  INTEGER                            :: ntransitions
  ! Specify backgraund model type (1-D, 2-D, 3-D)
  INTEGER                            :: model_type
  ! Specify minimal size of dynamic cell
  DOUBLE PRECISION, PARAMETER          :: minwidth = 1E8
  ! number of packets which will be saved into a file
  INTEGER                               :: n_add_pack
  INTEGER                               :: n_pack_save
  ! temporary file name
  CHARACTER(160)                     :: temp_filename = 'temp_packet'
  INTEGER                               :: tot_saved_packets
  ! for testing case
  LOGICAL                               :: simpleTrans, orbitals_nl
  INTEGER                               :: enable_diffusion
  INTEGER                               :: sobolev_approximation
  LOGICAL                               :: vel_propgrid = .false., vel_modgrid = .true.

  INTEGER, PARAMETER                    :: const_dimofspace = 3





END MODULE types
