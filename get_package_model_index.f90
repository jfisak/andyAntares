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
   get_package_model_index = dyn_cell(cur_cell_numb)%model_index

  RETURN 

 END FUNCTION
