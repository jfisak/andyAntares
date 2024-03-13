 FUNCTION get_package_model_index(pack_index)

! For the pack_index calculate the model_index

  USE types
USE constants

  IMPLICIT NONE    
    
    INTEGER                 :: pack_index, get_package_model_index
    INTEGER                     :: cur_cell_numb

   cur_cell_numb = package(pack_index)%cell_numb
   get_package_model_index = dyn_cell(cur_cell_numb)%model_index

  RETURN 

 END FUNCTION
