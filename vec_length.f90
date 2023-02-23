
 FUNCTION vec_length(vector)

  USE types
USE constants

  IMPLICIT NONE    

    DOUBLE PRECISION                 :: vec_length
    DOUBLE PRECISION, DIMENSION(3)   :: vector 

    vec_length=SQRT(vector(1)**2 + vector(2)**2 + vector(3)**2)

  RETURN

 END FUNCTION
