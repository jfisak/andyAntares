SUBROUTINE read_3D_model()

USE types
IMPLICIT NONE




SELECT CASE(inputmodel)


! testing pesudo 3D model
CASE(0)

 CALL read_3D_pseudo3D()

CASE(1)
 CALL read_3D_nico()
CASE DEFAULT
 write(*,*) 'read_3D_model: the choice: ', inputmodel, ' is not known'
 STOP
END SELECT

END SUBROUTINE read_3D_model
