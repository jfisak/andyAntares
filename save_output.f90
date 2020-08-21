! saving output
! 
! #00 OUTPUT FOLDER
! #01 STANDARD OUTPUT
! #02 SPECTRAL LINES
! #03 COUNTERS
! #04 TEMPERATURE STRUCTURE AND IONIZATION BALANCE
! #06 PACKETS INFORMATION
! #07 IONIZATION FRACTIONS
SUBROUTINE save_output(otype)

USE types
USE counters
IMPLICIT NONE


! type of output
INTEGER                                 :: otype
! date and time info
INTEGER, DIMENSION (9)            :: TT
! folder variables
! CHARACTER(LEN=30)                       :: outputfolder
LOGICAL                                 :: dirExists
CHARACTER(LEN=30)                       :: mkdirCMD
CHARACTER(LEN=60)                       :: lineOutput
! save informations about lines
INTEGER                                 :: I, J, K
DOUBLE PRECISION                        :: wavle
! occupation numbers
! levels index variables
INTEGER                                 :: act_elem, act_ion, act_lev
DOUBLE PRECISION                        :: act_pop
DOUBLE PRECISION                        :: eenergy
CHARACTER(LEN=60)                       :: fileTempStruct, fileOccNum, fileFreqs
CHARACTER(LEN=60)                       :: fileHydrogenFrac, fileHeliumFrac
CHARACTER(LEN=60)                       :: fileGrid, filePart
! ionization fraction files
DOUBLE PRECISION                        :: frac, N_jk, totElPop
! DOUBLE PRECISION                        :: frac1, N_jk1, totElPop1
! DOUBLE PRECISION                        :: frac2, N_jk2, totElPop2
! DOUBLE PRECISION                        :: frac3, N_jk3, totElPop3
INTEGER                                 :: indexe, indexi
INTEGER                                 :: status
CHARACTER(LEN=60)                       :: filePackets
! testing PoWR ionization fractions
DOUBLE PRECISION                        :: ntot, nhi
INTEGER, PARAMETER                      :: indexH = 1, indexHI = 1, indexHII = 2
INTEGER, PARAMETER                      :: indexHe = 2, indexHeI = 1, indexHeII = 2, indexHeIII = 3
DOUBLE PRECISION                        :: abundance, density
CHARACTER(LEN=60)                       :: fileHI, fileHII, fileHeI, fileHeII, fileHeIII
CHARACTER(LEN=60)                       :: fileEldens, fileRho
INTEGER                                 :: cell_index
DOUBLE PRECISION                        :: num_tot_pop
INTEGER                                 :: tot_n_ions, cur_ion, n_ions
DOUBLE PRECISION, ALLOCATABLE           :: part_functions(:)
DOUBLE PRECISION                        :: U, temperature
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
! creates a folder, where an output will be saved
! it reads a shell variable OUTPUTFO, if it does not
! exist, it will create (or not, if it already exists)
IF(outputfolder(:) == '') THEN
 CALL GET_ENVIRONMENT_VARIABLE("OUTPUTFO", outputfolder)
END IF
IF(outputfolder(:) == '') THEN
 CALL DATE_AND_TIME(VALUES = TT)
 ! write(*,*) 'save_output: ', TT(1), TT(2), TT(3), TT(4), TT(5), TT(6), TT(7)
 write(outputfolder, "(A6, I4.4, I2.2, I2.2, I2.2, I2.2, I2.2)") "3Dwind", &
  TT(1), TT(2), TT(3), TT(5), TT(6), TT(7)
 ! write(*,*) 'save_output: outputfolder = ', outputfolder
END IF

inquire( file=trim(outputfolder)//'/.', exist=dirExists )

IF(.NOT. dirExists) THEN
 mkdirCMD = 'mkdir '//TRIM(outputfolder)
 CALL SYSTEM(mkdirCMD)
 ! write(99,*) 'save_output: creating a folder: ', outputfolder
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
 write(98,*) lineOutput!trim(outputfolder), '/linevar.', my_rank, '.dat'
! CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
! STOP
 write(98,*) 'save_output: lineOutput = ', TRIM(lineOutput)
 OPEN(11,FILE=lineOutput)
  DO I = 1, ntransitions
   ! wavelength is in Angstroms
   wavle = 1e8 * light_speed / linelist(I)%freq
   WRITE(11,*) I, elements(linelist(I)%indexe)%atom_number, linelist(I)%indexi, wavle,&
    linelist(I)%A_ul, linelist(I)%n_int, linelist(I)%n_deexc
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
   WRITE(12, *) model_grid(I)%rwind / R_inf, model_grid(I)%T, model_grid(I)%e_dens
  END DO
 CLOSE(12)
 ! saving population numbers
 fileOccNum = trim(outputfolder)//'/occNums.dat'
 OPEN(13, FILE = fileOccNum)
  ! writing basic informations about chemical composition
  DO I = 1, n_modelgrid
   IF(model_grid(I)%assoc_cells == 0) CYCLE
   WRITE(13, *) '**modelCell**'
   WRITE(13, *) I, model_grid(I)%rwind, model_grid(I)%t, model_grid(I)%rho, model_grid(I)%e_dens
   WRITE(13, *) '**composition**'
   DO J = 1, n_elements
     WRITE(13, *) elements(J)%atom_number, elements(J)%abundance
   END DO
   ! write the occupation numbers
   WRITE(13, *) '**occunumbs**'
   DO act_elem = 1, n_elements
    DO act_ion = 1, SIZE(elements(act_elem)%ions)
     ! WRITE(13, *) model_grid(I)%grid_comp(act_elem)%grid_ion(act_ion)%tot_pop
     num_tot_pop = 0.D0
     DO act_lev = 1, SIZE(elements(act_elem)%ions(act_ion)%levels)
      eenergy = elements(act_elem)%ions(act_ion)%levels(act_lev)%exci_energy
      CALL populations(act_elem, act_ion, act_lev, I, act_pop)
      num_tot_pop = num_tot_pop + act_pop
      WRITE(13, *) elements(act_elem)%atom_number, act_ion, act_lev, &
      elements(act_elem)%ions(act_ion)%levels(act_lev)%stat_waight, &
      eenergy / e_v, act_pop
     END DO
     WRITE(13, *) elements(act_elem)%atom_number, act_ion, ' TOT ', &
     model_grid(I)%grid_comp(act_elem)%grid_ion(act_ion)%tot_pop / num_tot_pop, &
      model_grid(I)%grid_comp(act_elem)%grid_ion(act_ion)%tot_pop
    END DO
   END DO
  END DO
 CLOSE(13)
 fileFreqs = 'freqs.dat'
 OPEN(14, FILE = fileFreqs)
  DO I = 1, SIZE(package)
   IF(package(I)%typ == type_escaped) THEN
    WRITE(14, *) 1.D8 * light_speed / package(I)%freq_rf
   END IF
  END DO
 CLOSE(14)
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
 ! CALL SLEEP(120)
#endif
! CALL do_spectrum(SIZE(package))
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
    CALL part_fun(J, K, temperature, U)
    part_functions(cur_ion) = U
    cur_ion = cur_ion + 1
   END DO ! over ions
  END DO ! over elements
  write(17,*) I, part_functions(1:)
 END DO ! over model cells
 CLOSE(17)
CASE DEFAULT
 write(99,*) 'save_output: this case is not known'
END SELECT





END SUBROUTINE save_output
