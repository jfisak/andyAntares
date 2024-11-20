FUNCTION get_ion_number(input_n)

INTEGER                         :: input_n
CHARACTER(LEN=10)                :: get_ion_number

SELECT CASE(input_n)

CASE(0)
 get_ion_number = 'I'
CASE(1)
 get_ion_number = 'II'
CASE(2)
 get_ion_number = 'III'
CASE(3)
 get_ion_number = 'IV'
CASE(4)
 get_ion_number = 'V'
CASE(5)
 get_ion_number = 'VI'
CASE(6)
 get_ion_number = 'VII'
CASE(7)
 get_ion_number = 'VIII'
CASE(8)
 get_ion_number = 'IX'
CASE(9)
 get_ion_number = 'X'
CASE(10)
 get_ion_number = 'XI'
CASE(11)
 get_ion_number = 'XII'
CASE DEFAULT
 get_ion_number = '??'
END SELECT


END FUNCTION
