FUNCTION get_element_name(input_Z)

INTEGER                         :: input_Z
CHARACTER(LEN=2)                :: get_element_name

SELECT CASE(input_Z)

CASE(1)
 get_element_name = 'H'
CASE(2)
 get_element_name = 'He'
CASE(3)
 get_element_name = 'Li'
CASE(4)
 get_element_name = 'Be'
CASE(5)
 get_element_name = 'B'
CASE(6)
 get_element_name = 'C'
CASE(7)
 get_element_name = 'N'
CASE(8)
 get_element_name = 'O'
CASE(9)
 get_element_name = 'F'
CASE(10)
 get_element_name = 'Ne'
CASE(11)
 get_element_name = 'Na'
CASE DEFAULT
 get_element_name = '??'
END SELECT


END FUNCTION
