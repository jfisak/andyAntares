 FUNCTION get_package_model_index(pack_index)

! For the pack_index calculate the model_index

USE types
USE constants
USE dummypacket

  IMPLICIT NONE    
    
    INTEGER                 :: pack_index, get_package_model_index
    INTEGER                     :: cur_cell_numb, dummy_pack_index

   IF(pack_index <= SIZE(package)) THEN
    cur_cell_numb = package(pack_index)%cell_numb
   ELSE IF(pack_index > SIZE(package)) THEN
    dummy_pack_index = pack_index - SIZE(package)
    cur_cell_numb = dummypackage(dummy_pack_index)%cell_numb
   END IF
   IF(cur_cell_numb < 0) THEN
    ! write(*,*) 'get_package_model_index: pack_index = ', pack_index, 'cur_cell_numb = ', cur_cell_numb
    get_package_model_index = -99
    RETURN
    ! STOP 'get_package_model_index: cur_cell_numb < 0'
   END IF
   get_package_model_index = dyn_cell(cur_cell_numb)%model_index
   ! write(*,*) 'get_package_model_index: cur_cell_numb = ', cur_cell_numb
   ! write(*,*) 'get_package_model_index: get_package_model_index = ', get_package_model_index
   IF(dyn_cell(cur_cell_numb)%model_index == 0) THEN
    write(*,*) 'get_package_model_index: cur_cell_numb = ', cur_cell_numb, ' n_propgridcells = ', n_propgcells
    write(*,*) 'get_package_model_index: up_cell = ', dyn_cell(cur_cell_numb)%up_cell
    STOP 'get_package_model_index: model_index = 0'
   END IF

  RETURN 

 END FUNCTION
