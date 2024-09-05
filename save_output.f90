! saving output
! 
! #00 OUTPUT FOLDER
! #01 STANDARD OUTPUT
! #02 SPECTRAL LINES
! #03 COUNTERS
! #04 TEMPERATURE STRUCTURE AND IONIZATION BALANCE
! #06 PACKETS INFORMATION
! #07 IONIZATION FRACTIONS
! #10 VELOCITY FIELD IN THE PROPGRID CELLS
! #11 PROPMOD SLICES
SUBROUTINE save_output(otype)

USE types
USE constants
USE counters
IMPLICIT NONE


! type of output
INTEGER                                 :: otype
! folder variables
CHARACTER(LEN=file_length)                       :: lineOutput
! save informations about lines
INTEGER                                 :: I, J, K
DOUBLE PRECISION                        :: wavle
! occupation numbers
! levels index variables
INTEGER                                 :: act_elem, act_ion, act_lev
DOUBLE PRECISION                        :: act_pop
DOUBLE PRECISION                        :: eenergy
CHARACTER(LEN=file_length)                       :: fileTempStruct, fileOccNum
CHARACTER(LEN=file_length)                       :: fileHydrogenFrac, fileHeliumFrac
CHARACTER(LEN=file_length)                       :: fileGrid, filePart
! ionization fraction files
DOUBLE PRECISION                        :: frac, N_jk, totElPop
! DOUBLE PRECISION                        :: frac1, N_jk1, totElPop1
! DOUBLE PRECISION                        :: frac2, N_jk2, totElPop2
! DOUBLE PRECISION                        :: frac3, N_jk3, totElPop3
INTEGER                                 :: indexe, indexi
INTEGER                                 :: status
CHARACTER(LEN=file_length)                       :: filePackets
! testing PoWR ionization fractions
DOUBLE PRECISION                        :: ntot, nhi
INTEGER, PARAMETER                      :: indexH = 1, indexHI = 1, indexHII = 2
INTEGER, PARAMETER                      :: indexHe = 2, indexHeI = 1, indexHeII = 2, indexHeIII = 3
DOUBLE PRECISION                        :: abundance, density
CHARACTER(LEN=file_length)                       :: fileHI, fileHII, fileHeI, fileHeII, fileHeIII
CHARACTER(LEN=file_length)                       :: fileEldens, fileRho, temp_file_name
INTEGER                                 :: cell_index
DOUBLE PRECISION                        :: num_tot_pop
INTEGER                                 :: tot_n_ions, cur_ion, n_ions
DOUBLE PRECISION, ALLOCATABLE           :: part_functions(:)
DOUBLE PRECISION                        :: U, temperature
INTEGER                                 :: n_adgrids

! chemical composition
INTEGER                                 :: cur_indexe, cur_element, cur_Z, cur_nions
DOUBLE PRECISION                        :: cur_atom_mass, cur_abundance
CHARACTER(LEN=file_length)              :: cur_levelfile, cur_transfile

DOUBLE PRECISION, DIMENSION(3)          :: cur_vel, cur_centre
INTEGER                                 :: cur_index_I

!________________________________________________________________________________
! #11
INTEGER                                         :: Nx_cov, Ny_cov, Nz_cov, ind_I, ind_J
DOUBLE PRECISION                                :: x_cov, z_cov, cur_physvar
DOUBLE PRECISION, DIMENSION(const_dimofspace)   :: cur_pos
REAL, ALLOCATABLE                               :: coverage_matrix_T(:,:), coverage_matrix_rho(:,:), &
                                                   coverage_matrix_v(:,:)
INTEGER                                         :: cur_mgi, cur_pgi, cur_index
DOUBLE PRECISION                                :: cur_rho, cur_temp
CHARACTER(LEN=file_length)                      :: temp_file_name_t, temp_file_name_rho, temp_file_name_v 
DOUBLE PRECISION, DIMENSION(const_dimofspace)   :: width

!________________________________________________________________________________
! #00 output folder
!
! realizes if the output folder exists
!
! * if the output folder exists do nothing
! * if the output folder does not exist create a new one
!________________________________________________________________________________
#if mpi==1
 IF(my_rank == 0) THEN
#endif
  IF(outputfolder(:) == '') THEN
   CALL GET_ENVIRONMENT_VARIABLE("outputfolder", outputfolder)
  END IF
#if mpi==1
  DO I = 1, n_tasks - 1
   CALL MPI_SEND(outputfolder, 80, MPI_CHAR, I, 6, MPI_COMM_WORLD, ierr)
  END DO
 ELSE IF (outputfolder == '') THEN
  CALL MPI_RECV(outputfolder, 80, MPI_CHAR, 0, 6, MPI_COMM_WORLD, status, ierr)
 END IF
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! now will save important variables
SELECT CASE(otype)

!____________________________________________________________
! #01 standard output of procedures
!
! saves standard output into the file 99
!____________________________________________________________
CASE(0)
 ! write(*,*) 'save_output: outpfol = ', trim(outputfolder), ' my_rank = ', my_rank
 write(outputfile,"(A, A7, I3.3, A4)") trim(outputfolder), "/output", my_rank, '.dat'
 ! write(*,*) 'save_output: outputfile = ', outputfile
 ! inquire(unit=99, opened=itsopen)
 ! write(*,*) 'save_output: itsopen = ', itsopen
 ! IF(itsopen) THEN
 !  CLOSE(99)
 ! ELSE
  OPEN(99, FILE=outputfile) 
  ! write(*,*) 'save_output: opening the file ', outputfile
  write(99,*) '___________________________________________________________'
  write(99,*) '___________________________________________________________'
  write(99,*) 'ANDY ANTARES CODE'
  write(99,*) 'MŇAU'
  write(99,*) '___________________________________________________________'
  write(99,*) '___________________________________________________________'
 ! END IF
!____________________________________________________________
! #02 spectral lines
!
! saves info about spectral lines into the file 98
!
! atomic number, ion index, wavelength, lifetime, number of absorption in line,
! number of deexcitations in line
!____________________________________________________________
CASE(1)
 ! line rates
 write(lineOutput,"(A, A9, I3.3, A4)") trim(outputfolder), '/linevar.', my_rank, '.dat'
 OPEN(11,FILE=lineOutput)
  DO I = 1, ntransitions
   ! wavelength is in Angstroms
   wavle = 1e8 * light_speed / linelist(I)%freq
   WRITE(11,*) I, elements(linelist(I)%indexe)%atom_number, linelist(I)%indexi, wavle,&
    linelist(I)%lower, linelist(I)%upper, &
    linelist(I)%f_lu, linelist(I)%n_int, linelist(I)%n_deexc, linelist(I)%counted
  END DO
 CLOSE(11)
!________________________________________________________________________________
! #03 rate counters
! 
! saves number of processes for r, i and k packets
!
! 
!________________________________________________________________________________
CASE(2)
 write(99,*) 'number of processes'
 write(99,*) 'r-packets deactivation:'
 write(99,*) 'count_des_phot = ', count_des_phot, ' count_des_inte = ', count_des_inte, &
  ' count_des_esca = ', count_des_esca, ' count_des_resd = ', count_des_resd
 write(99,*) 'r-packets'
 write(99,*) 'count_r_line = ', count_r_line, ' count_r_thom = ', count_r_thom, &
  ' count_r_ph_k = ', count_r_ph_k, ' count_r_ph_i = ', count_r_ph_i, &
  ' count_r_ff = ', count_r_ff
 write(99,*) 'k-packets'
 write(99,*) 'count_cool_ex = ', count_cool_ex, ' count_cool_ff = ', count_cool_ff, &
  ' count_cool_io = ', count_cool_io, ' count_cool_fb = ', count_cool_fb
 write(99,*) 'i-packets'
 write(99,*) 'count_i_int_down = ', count_i_int_down, ' count_i_rad_dxrs = ', count_i_rad_dxrs,&
  ' count_i_rad_deex = ', count_i_rad_deex, &
  ' count_i_rad_dxfl = ', count_i_rad_dxfl, ' count_i_int_upwa = ', count_i_int_upwa, &
  ' count_i_col_deex = ', count_i_col_deex, ' count_i_int_phot = ', count_i_int_phot, &
  ' count_i_int_reco = ', count_i_int_reco, ' count_i_rad_reco = ', count_i_rad_reco, &
  ' count_i_col_reco = ', count_i_col_reco
 write(99,*) 'd-packets'
 write(99,*) 'count_d_change_cell = ', count_d_change_cell, ' count_d_new_choice = ', count_d_new_choice, &
  ' count_d_radiative = ', count_d_radiative, ' count_d_rad_end = ', count_d_rad_end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! #04 temperature structure and ionization balance
!
! temperature structure into file 12
!
! radius, temperature
!
! ionization structure into file 13
!
!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NOT WORKING FOR MPI YET
CASE(3)
 fileTempStruct = trim(outputfolder)//'/tempStruct.dat'
 ! write(fileTempStruct,"(A6,I3.3)") trim(outputfolder), tem
 OPEN(12, FILE=fileTempStruct)
  DO I = 1, n_modelgrid
   IF(model_grid(I)%assoc_cells == 0) CYCLE
   WRITE(12, *) model_grid(I)%rwind, model_grid(I)%T, model_grid(I)%e_dens
  END DO
 CLOSE(12)
 ! saving population numbers
 fileOccNum = trim(outputfolder)//'/occNums.dat'
 OPEN(13, FILE = fileOccNum)
  ! writing basic informations about chemical composition
  ! 000 -- model cell information
  ! 001 -- chemical composition
  ! 002 -- level population
  ! 003 -- total ion population
  DO I = 1, n_modelgrid
   IF(model_grid(I)%assoc_cells == 0) CYCLE
   WRITE(13, *) '000', I, model_grid(I)%rwind, model_grid(I)%t, model_grid(I)%rho, model_grid(I)%e_dens
   DO J = 1, n_elements
     WRITE(13, *) '001', elements(J)%atom_number, elements(J)%abundance
   END DO
   ! write the occupation numbers
   DO act_elem = 1, n_elements
    DO act_ion = 1, SIZE(elements(act_elem)%ions)
     ! WRITE(13, *) model_grid(I)%grid_comp(act_elem)%grid_ion(act_ion)%tot_pop
     num_tot_pop = 0.D0
     DO act_lev = 1, SIZE(elements(act_elem)%ions(act_ion)%levels)
      eenergy = elements(act_elem)%ions(act_ion)%levels(act_lev)%exci_energy
      CALL populations(act_elem, act_ion, act_lev, I, act_pop)
      num_tot_pop = num_tot_pop + act_pop
      WRITE(13, *) '002', elements(act_elem)%atom_number, act_ion, act_lev, &
      elements(act_elem)%ions(act_ion)%levels(act_lev)%stat_waight, &
      eenergy / e_v, act_pop
     END DO
     WRITE(13, *) '003', elements(act_elem)%atom_number, act_ion, ' TOT ', &
     model_grid(I)%grid_comp(act_elem)%grid_ion(act_ion)%tot_pop / num_tot_pop, &
      model_grid(I)%grid_comp(act_elem)%grid_ion(act_ion)%tot_pop
    END DO
   END DO
  END DO
 CLOSE(13)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! #04 IONIZATION BALANCE
!
! saving ionization balance for hydrogen and helium
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE(4)
 fileHydrogenFrac = trim(outputfolder)//'/hydrogenFrac.dat'
 OPEN(14, FILE=fileHydrogenFrac)
  DO I = 1, n_modelgrid
   indexe = 1
   indexi = 1
   IF(model_grid(I)%assoc_cells == 0) CYCLE
   totElPop = model_grid(I)%rho * model_grid(I)%grid_comp(indexe)%abund / &
    elements(indexe)%atom_mass
   N_jk = model_grid(I)%grid_comp(indexe)%grid_ion(indexi)%tot_pop
   frac = N_jk / totElPop
   write(14,*) model_grid(I)%rwind / R_inf, frac
  END DO
 CLOSE(14)
 fileHeliumFrac = trim(outputfolder)//'/heliumFrac.dat'
 ! OPEN(15, FILE=fileHeliumFrac)
 !  DO I = 1, n_modelgrid
 !  IF(model_grid(I)%assoc_cells == 0) CYCLE
 !  indexe = 2
 !   totElPop1 = model_grid(I)%rho * model_grid(I)%grid_comp(indexe)%abund / &
 !    elements(indexe)%atom_mass
 !   N_jk1 = model_grid(I)%grid_comp(indexe)%grid_ion(1)%tot_pop
 !   frac1 = N_jk1 / totElPop1
 !   totElPop2 = model_grid(I)%rho * model_grid(I)%grid_comp(indexe)%abund / &
 !    elements(indexe)%atom_mass
 !   N_jk2 = model_grid(I)%grid_comp(indexe)%grid_ion(2)%tot_pop
 !   frac2 = N_jk2 / totElPop2
 !   totElPop3 = model_grid(I)%rho * model_grid(I)%grid_comp(indexe)%abund / &
 !    elements(indexe)%atom_mass
 !   N_jk3 = model_grid(I)%grid_comp(indexe)%grid_ion(3)%tot_pop
 !   frac3 = N_jk3 / totElPop3
 !  write(15,*) model_grid(I)%rwind / R_inf, frac1, frac2, frac3
 !  END DO
 ! CLOSE(15)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! #05 packets informations
! 
! saves rf frequency and energy for the spectrum generation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE(5)
 write(filePackets,"(A, A7, I3.3, A4)") trim(outputfolder),&
  "/packets", my_rank, '.dat'
 write(99,*) 'save_output: filePackets = ', filePackets
 OPEN(98, FILE=filePackets)
  DO I = 1, SIZE(package) - 1
   write(98,*) package(I)%typ, package(I)%freq_rf, package(I)%e_rf
  END DO
 CLOSE(98)
#if mpi==1
 CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
#endif
CALL do_spectrum(SIZE(package))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! #07 ionization fractions
!
! only for testing 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE(6)
 fileHI = trim(outputfolder)//'/HI.dat'
 ! write(*,*) 'save_output: fileHI'
 OPEN(40, FILE=fileHI)
  DO I = 1, n_modelgrid
   IF(model_grid(I)%assoc_cells == 0) CYCLE
   abundance = model_grid(I)%grid_comp(indexH)%abund
   density = model_grid(I)%rho
   ntot = abundance * density / elements(indexH)%atom_mass
   nhi = model_grid(I)%grid_comp(indexH)%grid_ion(indexHI)%gl_pop
   cell_index = n_modelgrid - I + 1
   write(40,*) I, log10(ntot), log10(nhi/ntot)
  END DO
 CLOSE(40)
 ! write(*,*) 'save_output: fileHII'
 fileHII = trim(outputfolder)//'/HII.dat'
 OPEN(40, FILE=fileHII)
  DO I = 1, n_modelgrid
   IF(model_grid(I)%assoc_cells == 0) CYCLE
   abundance = model_grid(I)%grid_comp(indexH)%abund
   density = model_grid(I)%rho
   ntot = abundance * density / elements(indexH)%atom_mass
   nhi = model_grid(I)%grid_comp(indexH)%grid_ion(indexHII)%gl_pop
   cell_index = n_modelgrid - I + 1
   write(40,*) I, log10(ntot), log10(nhi/ntot)
  END DO
 CLOSE(40)
 ! write(*,*) 'save_output: fileHeI'
 fileHeI = trim(outputfolder)//'/HeI.dat'
 OPEN(40, FILE=fileHeI)
  DO I = 1, n_modelgrid
   IF(model_grid(I)%assoc_cells == 0) CYCLE
   abundance = model_grid(I)%grid_comp(indexHe)%abund
   density = model_grid(I)%rho
   ntot = abundance * density / elements(indexHe)%atom_mass
   nhi = model_grid(I)%grid_comp(indexHe)%grid_ion(indexHeI)%gl_pop
   cell_index = n_modelgrid - I + 1
   write(40,*) I, log10(ntot), log10(nhi/ntot)
  END DO
 CLOSE(40)
 ! write(*,*) 'save_output: fileHeII'
 fileHeII = trim(outputfolder)//'/HeII.dat'
 OPEN(40, FILE=fileHeII)
  DO I = 1, n_modelgrid
   IF(model_grid(I)%assoc_cells == 0) CYCLE
   abundance = model_grid(I)%grid_comp(indexHe)%abund
   density = model_grid(I)%rho
   ntot = abundance * density / elements(indexHe)%atom_mass
   nhi = model_grid(I)%grid_comp(indexHe)%grid_ion(indexHeII)%gl_pop
   cell_index = n_modelgrid - I + 1
   write(40,*) I, log10(ntot), log10(nhi/ntot)
  END DO
 CLOSE(40)
 ! write(*,*) 'save_output: fileHeIII'
 fileHeIII = trim(outputfolder)//'/HeIII.dat'
 OPEN(40, FILE=fileHeIII)
  DO I = 1, n_modelgrid
   IF(model_grid(I)%assoc_cells == 0) CYCLE
   abundance = model_grid(I)%grid_comp(indexHe)%abund
   density = model_grid(I)%rho
   ntot = abundance * density / elements(indexHe)%atom_mass
   ! nhi = model_grid(I)%grid_comp(indexHe)%grid_ion(indexHeIII)%gl_pop
   nhi = model_grid(I)%grid_comp(indexHe)%grid_ion(indexHeIII)%tot_pop
   cell_index = n_modelgrid - I + 1
   write(40,*) I, log10(ntot), log10(nhi/ntot)
  END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! #07 electron density and mass density
!
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE(7)
 fileEldens = trim(outputfolder)//'/elDens.dat'
 OPEN(40, FILE=fileEldens)
  DO I = 1, n_modelgrid
   IF(model_grid(I)%assoc_cells == 0) CYCLE
   write(40,'(I3,d12.4,d12.4)') I, model_grid(I)%rwind, model_grid(I)%e_dens
  END DO
 CLOSE(40)
 fileRho = trim(outputfolder)//'/rho.dat'
 OPEN(40, FILE=fileRho)
  DO I = 1, n_modelgrid
   IF(model_grid(I)%assoc_cells == 0) CYCLE
   write(40,*) I, model_grid(I)%rwind/R_star, model_grid(I)%rho
  END DO
 CLOSE(40)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! #08 propagation grid information
!
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE(8)
 fileGrid=trim(outputfolder)//'/dyn_cells.dat'
 OPEN(16,FILE=fileGrid)
 write(*,*) 'setup_grid2: SAVING CELLS INTO A FILE dyn_cells.dat'
  DO I = 1, SIZE(dyn_cell)
   write(16,*) I, dyn_cell(I)%corner, dyn_cell(I)%width, dyn_cell(I)%up_cell, dyn_cell(I)%model_index
  END DO
 CLOSE(16)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! #09 partition function
!
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE(9)
 filePart=trim(outputfolder)//'/partitionFunction.dat'
 tot_n_ions = 0
 DO J = 1, n_elements
  n_ions = SIZE(elements(J)%ions)
  tot_n_ions = tot_n_ions + n_ions
 END DO
 ALLOCATE(part_functions(tot_n_ions))
 OPEN(17, FILE=filePart)
 DO I = 1, n_modelgrid
  temperature = model_grid(I)%T
  cur_ion = 1
  DO J = 1, n_elements
   n_ions = SIZE(elements(J)%ions)
   DO K = 1, n_ions
    IF(temperature /= 0.D0) THEN
     CALL part_fun(J, K, temperature, U)
     part_functions(cur_ion) = U
    ELSE
     part_functions(cur_ion) = 0.D0
    END IF
    cur_ion = cur_ion + 1
   END DO ! over ions
  END DO ! over elements
  write(17,*) I, part_functions(1:)
 END DO ! over model cells
 CLOSE(17)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! #09 adaptive grid parts
!
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE(10)
 write(temp_file_name,"(A, A12, I3.3, A4)") TRIM(outputfolder), '/temp_adgrid', my_rank, ".dat"
 write(*,*) 'save_output: temp_file_name = ', temp_file_name
 Ngrid = nx_cell * ny_cell * nz_cell
 n_adgrids = SIZE(dyn_cell)
 IF(Ngrid == n_adgrids) RETURN
 OPEN(72, form='unformatted', FILE=temp_file_name)
  DO I = Ngrid, n_adgrids
   WRITE(72) dyn_cell(I)
  END DO
 CLOSE(72)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! #10 velocity field in the propGrid cells
!
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE(11)
 write(temp_file_name,"(A, A11, I3.3, A4)") TRIM(outputfolder), '/velo_field', my_rank, ".dat"
 write(*,*) 'save_output: temp_file_name = ', temp_file_name
 OPEN(73, FILE=temp_file_name)
  DO cur_index_I = 1, n_propgcells
   cur_centre = dyn_cell(cur_index_I)%corner + dyn_cell(cur_index_I)%width/2.0
   cur_vel = dyn_cell(cur_index_I)%vec_vel
   write(73,*) cur_centre, cur_vel
  END DO
 CLOSE(73)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! #11 density and temperature slices
!
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE(12)
 
# if mpi == 1
 IF(my_rank == 0) THEN
# endif
! maps of a coverage
! creates a matrix with a specified physical quantity
IF(dyngrid == 0) THEN
 Nx_cov = nx_cell
 Ny_cov = ny_cell
 Nz_cov = nz_cell
ELSE IF(dyngrid > 0) THEN
 width(ind_x) = MINVAL(dyn_cell(:)%width(ind_x))
 width(ind_y) = MINVAL(dyn_cell(:)%width(ind_y))
 width(ind_z) = MINVAL(dyn_cell(:)%width(ind_z))
 Nx_cov = (xmax - xmin)/width(ind_x)
 Ny_cov = (ymax - ymin)/width(ind_y)
 Nz_cov = (zmax - zmin)/width(ind_z)
END IF
x_cov = 0.D0
z_cov = 0.D0

ALLOCATE(coverage_matrix_rho(Nx_cov, Ny_cov), coverage_matrix_t(Nx_cov, Ny_cov), &
  coverage_matrix_v(Nx_cov * Ny_cov, 2 * const_dimofspace))
coverage_matrix_rho(:,:) = -1.0
coverage_matrix_t(:,:) = -1.0
coverage_matrix_v(:,:) = -1.0

write(temp_file_name_t,"(A, A17)") TRIM(outputfolder), '/propmod_t_xy.dat'
write(temp_file_name_rho,"(A, A19)") TRIM(outputfolder), '/propmod_rho_xy.dat'
write(temp_file_name_v,"(A, A19)") TRIM(outputfolder), '/propmod_vel_xy.dat'
write(*,*) 'save_output: temp_file_name_t = ', temp_file_name_t
! probably a temporary solution
cur_pos(ind_z) = z_cov
! xy plane
cur_index = 0
DO ind_I = 1, Nx_cov
 cur_pos(ind_x) = ((xmax - xmin) * ind_I + (Nx_cov * xmin - xmax))/DBLE(Nx_cov - 1)
 DO ind_J = 1, Ny_cov
  cur_index = cur_index + 1
  cur_pos(ind_y) = ((ymax - ymin) * ind_J + (Ny_cov * ymin - ymax))/DBLE(Ny_cov - 1)
  ! looking for a current propGrid cell index
  CALL find_dyn_cell1(cur_pos, cur_pgi)
  IF(cur_pgi > 0) THEN
   cur_mgi = dyn_cell(cur_pgi)%model_index
   IF(cur_mgi > 0) THEN
    cur_temp = model_grid(cur_mgi)%T
    cur_rho = model_grid(cur_mgi)%rho
    coverage_matrix_T(ind_I, ind_J) = REAL(cur_temp)
    coverage_matrix_rho(ind_I, ind_J) = REAL(cur_rho)
    coverage_matrix_v(cur_index, 1:3) = REAL(cur_pos)
    CALL velo_vector(cur_pos, cur_mgi, cur_vel)
    coverage_matrix_v(cur_index, 4:6) = REAL(cur_vel)
   END IF
  END IF ! cur_pgi > 0
 END DO
END DO

OPEN(173, FILE=temp_file_name_t)
OPEN(174, FILE=temp_file_name_rho)
OPEN(175, FILE=temp_file_name_v)

DO ind_I = 1, Ny_cov
 write(173,*) coverage_matrix_T(:,ind_I)
 write(174,*) coverage_matrix_rho(:,ind_I)
END DO
DO ind_I = 1, Nx_cov * Ny_cov
 write(175,*) coverage_matrix_v(ind_I,:)
END DO

CLOSE(173)
CLOSE(174)
CLOSE(175)

coverage_matrix_rho(:,:) = -1.0
coverage_matrix_t(:,:) = -1.0
coverage_matrix_v(:,:) = -1.0

write(temp_file_name_t,"(A, A17)") TRIM(outputfolder), '/propmod_t_yz.dat'
write(temp_file_name_rho,"(A, A19)") TRIM(outputfolder), '/propmod_rho_yz.dat'
write(temp_file_name_v,"(A, A19)") TRIM(outputfolder), '/propmod_vel_yz.dat'

! probably a temporary solution
cur_pos(ind_x) = x_cov
! yz plane
cur_index = 0
DO ind_I = 1, Ny_cov
 cur_pos(ind_y) = ((ymax - ymin) * ind_I + (Ny_cov * ymin - ymax))/DBLE(Ny_cov - 1)
 DO ind_J = 1, Ny_cov
  cur_index = cur_index + 1
  cur_pos(ind_z) = ((zmax - zmin) * ind_J + (Nz_cov * zmin - zmax))/DBLE(Nz_cov - 1)
  ! looking for a current propGrid cell index
  CALL find_dyn_cell1(cur_pos, cur_pgi)
  IF(cur_pgi > 0) THEN
   cur_mgi = dyn_cell(cur_pgi)%model_index
   IF(cur_mgi > 0) THEN
    cur_temp = model_grid(cur_mgi)%T
    cur_rho = model_grid(cur_mgi)%rho
    coverage_matrix_T(ind_I, ind_J) = REAL(cur_temp)
    coverage_matrix_rho(ind_I, ind_J) = REAL(cur_rho)
    coverage_matrix_v(cur_index, 1:3) = REAL(cur_pos)
    CALL velo_vector(cur_pos, cur_mgi, cur_vel)
    coverage_matrix_v(cur_index, 4:6) = REAL(cur_vel)
   END IF
  END IF ! cur_pgi > 0
 END DO
END DO

OPEN(173, FILE=temp_file_name_t)
OPEN(174, FILE=temp_file_name_rho)
OPEN(175, FILE=temp_file_name_v)

DO ind_I = 1, Ny_cov
 write(173,*) coverage_matrix_T(:,ind_I)
 write(174,*) coverage_matrix_rho(:,ind_I)
END DO
DO ind_I = 1, Ny_cov * Nz_cov
 write(175,*) coverage_matrix_v(ind_I,:)
END DO

CLOSE(173)
CLOSE(174)
CLOSE(175)

# if mpi == 1
 END IF
# endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! #100 input file
! 
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE(100)
 write(99,*) '___________________________________________________________'
 write(99,*) 'input file'
 write(99,*) '___________________________________________________________'
 write(99,*) 'set up variables:'
 ! write(99,*) 'number of packets: ', n_pack, ' a temporary file saves ', n_pack_save, ' packets'
 write(99,*) ' a temporary file saves ', n_pack_save, ' packets'
 write(99,*) 'n_nubin = ', n_nubin
 write(99,*) 'propagation grid parameters:'
 write(99,*) 'using a previously saved grid = ', saved_grid
 write(99,*) 'nx_cell = ', nx_cell, ' ny_cell = ', ny_cell, ' nz_cell = ', nz_cell
 write(99,*) 'xmax = ', xmax, ' ymax = ', ymax, ' zmax = ', zmax
 write(99,*) 'adaptive grid type = ', dyngrid
 write(99,*) 'number of virtual points = ', Nvirtpoint
 write(99,*) 'model grid parameters:'
 write(99,*) 'input model file = ', inputmodelFile
 write(99,*) 'input composition = ', inputcomposition
 write(99,*) 'input population file = ', inputpopfile
 write(99,*) 'electron density file = ', eldensfile
 write(99,*) 'model_type = ', model_type, 'inputmodel = ', inputmodel
 write(99,*) 'NLTE = ', nlte, ' velocity approximation = ', velApprox
 write(99,*) 'absortive surface = ', abs_surface
 write(99,*) 'diffusion = ', enable_diffusion
 write(99,*) 'calculate BRTM = ', calc_brtm
 write(99,*) '___________________________________________________________'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! #101 composition file
! 
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE(101)
 write(99,*) '___________________________________________________________'
 write(99,*) 'composition file'
 write(99,*) '___________________________________________________________'
 write(99,*) '1.) chemical compositiion'
 write(99,*) 'indexe    Z       n_ions  atom_mass       tot_abundance   level file      transition file'

 DO cur_element = 1, n_elements
  cur_indexe = elements(cur_element)%indexe
  cur_Z = elements(cur_element)%atom_number
  cur_nions = elements(cur_element)%nions
  cur_atom_mass = elements(cur_element)%atom_mass
  cur_abundance = elements(cur_element)%abundance
  cur_levelfile = elements(cur_element)%levelfile
  cur_transfile = elements(cur_element)%transitionfile
  write(99,*) cur_indexe, cur_Z, cur_nions, cur_atom_mass, cur_abundance, cur_levelfile, cur_transfile
 END DO
 write(99,*) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! #102 model grid description
! 
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CASE(102)
 
 write(99,*) '__________________________________________________'
 write(99,*) '_________MODEL GRID DESCRIPTION___________________'
 write(99,*) '__________________________________________________'
 write(99,*) 'R_star = ', R_star/R_sun, 'R_inf = ', R_inf/R_sun

 write(99,*) '__________________________________________________'
 
CASE DEFAULT
 write(99,*) 'save_output: this case is not known'
END SELECT





END SUBROUTINE save_output
