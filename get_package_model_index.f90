 FUNCTION get_package_model_index(pack_index)

! For the pack_index calculate the model_index

  USE types

  IMPLICIT NONE    
    
    INTEGER                 :: pack_index, get_package_model_index

   get_package_model_index = cell(package(pack_index)%cell_numb)%model_index

  RETURN 

 END FUNCTION
