SUBROUTINE save_output(otype)

USE types
USE counters
IMPLICIT NONE


! type of output
INTEGER                                 :: otype
! folder variables
CHARACTER(LEN=30)                       :: outputfolder
LOGICAL                                 :: dirExists
CHARACTER(LEN=30)                       :: mkdirCMD
CHARACTER(LEN=60)                       :: lineOutput
! save informations about lines
INTEGER                                 :: I, J
DOUBLE PRECISION                        :: wavle
! occupation numbers
! levels index variables
INTEGER                                 :: act_elem, act_ion, act_lev
DOUBLE PRECISION                        :: act_pop
DOUBLE PRECISION                        :: eenergy
CHARACTER(LEN=60)                       :: fileTempStruct, fileOccNum, fileFreqs
CHARACTER(LEN=60)                       :: filePackets

! creates a folder, where an output will be saved
! it reads a shell variable OUTPUTFO, if it does not
! exist, it will create (or not, if it already exists)
! a directory 3dwindmodel
CALL GET_ENVIRONMENT_VARIABLE("OUTPUTFO", outputfolder)
IF(outputfolder(:) == '') THEN
 outputfolder = '3dwindmodel'
END IF

inquire( file=trim(outputfolder)//'/.', exist=dirExists )

IF(.NOT. dirExists) THEN
 mkdirCMD = 'mkdir '//TRIM(outputfolder)
 CALL SYSTEM(mkdirCMD)
 write(*,*) 'save_output: creating a folder: ', outputfolder
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! now will save important variables
SELECT CASE(otype)
CASE(1)
 ! line rates
 write(lineOutput,"(A, A9, I3.3, A4)") trim(outputfolder), '/linevar.', my_rank, '.dat'
 write(*,*) lineOutput!trim(outputfolder), '/linevar.', my_rank, '.dat'
! CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
! STOP
 write(*,*) 'save_output: lineOutput = ', TRIM(lineOutput)
 OPEN(11,FILE=lineOutput)
  DO I = 1, ntransitions
   ! wavelength is in Angstroms
   wavle = 1e8 * light_speed / linelist(I)%freq
   WRITE(11,*) elements(linelist(I)%indexe)%atom_number, linelist(I)%indexi, wavle,&
    linelist(I)%A_ul, linelist(I)%n_int, linelist(I)%n_deexc
  END DO
 CLOSE(11)
!________________________________________________________________________________
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! RATE COUNTERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!________________________________________________________________________________
CASE(2)
 write(*,*) 'save_output'
 write(*,*) 'count_cool_ex = ', count_cool_ex, ' count_cool_ff = ', count_cool_ff, &
  ' count_cool_io = ', count_cool_io, ' count_cool_fb = ', count_cool_fb
 write(*,*) 'count_i_int_down = ', count_i_int_down, ' count_i_rad_dxrs = ', count_i_rad_dxrs,&
  ' count_i_rad_deex = ', count_i_rad_deex, &
  ' count_i_rad_dxfl = ', count_i_rad_dxfl, ' count_i_int_upwa = ', count_i_int_upwa, &
  ' count_i_col_deex = ', count_i_col_deex, ' count_i_int_phot = ', count_i_int_phot, &
  ' count_i_int_reco = ', count_i_int_reco, ' count_i_rad_reco = ', count_i_rad_reco, &
  ' count_i_col_reco = ', count_i_col_reco
 write(*,*) 'count_r_line = ', count_r_line, ' count_r_thom = ', count_r_thom, &
  ' count_r_ph_k = ', count_r_ph_k, ' count_r_ph_i = ', count_r_ph_i, &
  ' count_r_ff = ', count_r_ff
 write(*,*) 'count_des_phot = ', count_des_phot, ' count_des_inte = ', count_des_inte, &
  ' count_des_esca = ', count_des_esca
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! TEMPERATURE STRUCTURE AND IONIZATION BALANCE !!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NOT WORKING FOR MPI YET
CASE(3)
 fileTempStruct = trim(outputfolder)//'/tempStruct.dat'
 ! write(fileTempStruct,"(A6,I3.3)") trim(outputfolder), tem
 OPEN(12, FILE=fileTempStruct)
  DO I = 1, n_modelgrid
   WRITE(12, *) model_grid(I)%rwind, model_grid(I)%T
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
     WRITE(13, *) model_grid(I)%grid_comp(act_elem)%grid_ion(act_ion)%tot_pop
     DO act_lev = 1, SIZE(elements(act_elem)%ions(act_ion)%levels)
      eenergy = elements(act_elem)%ions(act_ion)%levels(act_lev)%exci_energy
      CALL populations(act_elem, act_ion, act_lev, I, act_pop)
      WRITE(13, *) elements(act_elem)%atom_number, act_ion, act_lev, eenergy, act_pop
     END DO
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
! PACKETS INFORMATIONS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE(4)
 write(filePackets,"(A, A7, I3.3, A4)") trim(outputfolder),&
  "/packets", my_rank, '.dat'
 write(*,*) 'save_output: filePackets = ', filePackets
 OPEN(98, FILE=filePackets)
  DO I = 1, SIZE(package) - 1
   write(98,*) package(I)%typ, package(I)%freq_rf, package(I)%e_rf
  END DO
 CLOSE(98)
#if mpi==1
 CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
 ! CALL SLEEP(120)
#endif
 CALL do_spectrum(SIZE(package), outputfolder)
CASE DEFAULT
 write(*,*) 'save_output: this case is not known'
END SELECT





END SUBROUTINE save_output
