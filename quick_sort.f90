SUBROUTINE sort_quick (customer_list, nvals) 
!Sets up for the quick sort recursive method
	USE types
USE constants
	LOGICAL, EXTERNAL :: gt_person !funtion to tell which person goes first
	INTEGER :: i = 0
	INTEGER :: j = 0
	TYPE(person):: temp_person2
	TYPE(person), INTENT(INOUT), DIMENSION(*) :: customer_list !Define the customer list to handle whatever size is sent
	INTEGER, INTENT(IN) :: nvals !grab the number of values from the calling code

	! Here is where we will do the selection sort

	WRITE (*,*) 'Now doing Quick sort'
	CALL qsRecursive(0, nvals-1, customer_list) !kicks off the recursive process
	

END SUBROUTINE sort_quick

SUBROUTINE sort_quick (customer_list, nvals) 
!Sets up for the quick sort recursive method
	USE types
USE constants
	LOGICAL, EXTERNAL :: gt_person !funtion to tell which person goes first
	INTEGER :: i = 0
	INTEGER :: j = 0
	TYPE(person):: temp_person2
	TYPE(person), INTENT(INOUT), DIMENSION(*) :: customer_list !Define the customer list to handle whatever size is sent
	INTEGER, INTENT(IN) :: nvals !grab the number of values from the calling code

	! Here is where we will do the selection sort

	WRITE (*,*) 'Now doing Quick sort'
	CALL qsRecursive(0, nvals-1, customer_list) !kicks off the recursive process
	

END SUBROUTINE sort_quick

RECURSIVE SUBROUTINE qsRecursive (lo, hi, customer_list)
!This is the actualy recursive portion of the quicksort
	USE types
USE constants
	INTEGER :: pivotPoint
	INTEGER, INTENT(IN) :: lo
	INTEGER, INTENT(IN) :: hi
	TYPE(person), INTENT(INOUT), DIMENSION(*) :: customer_list
	pivotPoint = qsPartition(lo, hi, customer_list); !basically all we do is find the pivot point, adjust elements, then call it again
	IF (lo < pivotPoint) CALL qsRecursive(lo, pivotPoint -1, customer_list)
	IF (pivotPoint < hi) CALL qsRecursive(pivotPoint + 1, hi, customer_list)

END SUBROUTINE qsRecursive

FUNCTION qsPartition (loin, hiin, customer_list)
	!The partition portios of the Quick Sort is the must involved part
	USE types
USE constants
	LOGICAL, EXTERNAL :: gt_person !funtion to tell which person goes first
	TYPE(person), INTENT(INOUT), DIMENSION(*) :: customer_list
	INTEGER, INTENT(IN) :: loin
	INTEGER:: lo !variable so we can manipulate the hi and lo values without changing things elsewhere in the program by reference
	INTEGER, INTENT(IN) :: hiin
	INTEGER:: hi !variable so we can manipulate the hi and lo values without changing things elsewhere in the program by reference
	TYPE(person)::pivot !the temp location for the pivitoal element to which everything will be compaired
	hi = hiin
	lo = loin
	pivot = customer_list(lo)
	DO
		IF (lo >= hi) EXIT !exit the loop when done
		DO !move in from the right
			IF ((gt_person(pivot, customer_list(hi))) .OR. (lo >= hi)) EXIT
			hi = hi - 1
		END DO	
		IF (hi /= lo) then !move the entry indexed by hi to left side of partition
				customer_list(lo) = customer_list(hi) 
				lo = lo + 1
		END IF
		DO !move in from the left
			IF ((gt_person(customer_list(lo),pivot)) .OR. (lo >= hi)) EXIT
			lo = lo + 1
		END DO	
		IF (hi /= lo) then !move the entry indexed by hi to left side of partition
			customer_list(hi) = customer_list(lo) 
			hi = hi - 1
		END IF
	END DO
	customer_list(hi) = pivot !put the pivot element back when we're done
	qsPartition = hi !return the correct position of the pivot element
END FUNCTION qsPartition
