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
CHARACTER(LEN=60)                       :: fileTempStruct, fileOccNum

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
 lineOutput = trim(outputfolder)//'/linevar.dat'
 write(*,*) 'save_output: lineOutput = ', TRIM(lineOutput)
 OPEN(11,FILE=lineOutput)
  DO I = 1, ntransitions
   ! wavelength is in Angstroms
   wavle = 1e8 * light_speed / linelist(I)%freq
   WRITE(11,*) elements(linelist(I)%indexe)%atom_number, linelist(I)%indexi, wavle,&
    linelist(I)%A_ul, linelist(I)%n_int
  END DO
 CLOSE(11)
! rate counters
CASE(2)
 write(*,*) 'save_output'
 write(*,*) 'count_cool_ex = ', count_cool_ex, ' count_cool_ff = ', count_cool_ff, &
  ' count_cool_io = ', count_cool_io, ' count_cool_fb = ', count_cool_fb
 write(*,*) 'count_i_int_down = ', count_i_int_down, ' count_i_rad_dxrs = ', count_i_rad_dxrs,&
  ' count_i_rad_dxfl = ', count_i_rad_dxfl, ' count_i_int_upwa = ', count_i_int_upwa, &
  ' count_i_col_deex = ', count_i_col_deex, ' count_i_int_phot = ', count_i_int_phot, &
  ' count_i_int_reco = ', count_i_int_reco, ' count_i_rad_reco = ', count_i_rad_reco, &
  ' count_i_col_reco = ', count_i_col_reco
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! TEMPERATURE STRUCTURE AND IONIZATION BALANCE !!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE(3)
 fileTempStruct = trim(outputfolder)//'/tempStruct.dat'
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
   WRITE(13, *) I, model_grid(I)%rwind
   WRITE(13, *) '**composition**'
   DO J = 1, n_elements
     WRITE(13, *) elements(J)%atom_number, elements(J)%abundance
   END DO
   ! write the occupation numbers
   WRITE(13, *) '**occunumbs**'
   DO act_elem = 1, n_elements
    DO act_ion = 1, SIZE(elements(act_elem)%ions)
     DO act_lev = 1, SIZE(elements(act_elem)%ions(act_ion)%levels)
      eenergy = elements(act_elem)%ions(act_ion)%levels(act_lev)%exci_energy
      CALL populations(act_elem, act_ion, act_lev, I, act_pop)
      WRITE(13, *) elements(act_elem)%atom_number, act_ion, act_lev, eenergy, act_pop
     END DO
    END DO
   END DO
  END DO
 CLOSE(13)
CASE DEFAULT
 write(*,*) 'save_output: this case is not known'
END SELECT






END SUBROUTINE save_output
