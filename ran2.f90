DOUBLE PRECISION FUNCTION  ran2(idum)

 USE ran2_class
 IMPLICIT NONE 

 ! input variables
 INTEGER                     :: idum
 INTEGER                     :: j,k
 INTEGER                     :: iy, idum2
 INTEGER, DIMENSION(NTAB)    :: iv
 DOUBLE PRECISION            :: out_ran

 IF(.NOT. ASSOCIATED(ranum)) THEN
  !write(*,*) 'ran2: allocating ranum'
  ALLOCATE(ranum)
 END IF
 !write(*,*) 'ran2'
 iy = ranum%iy
 idum2 = ranum%idum2
 iv(:) = ranum%iv(:)
 !write(*,*) 'ran2: iy = ', iy, ' idum2 = ', idum2

 if (idum.le.0) then
  idum=max(-idum,1)
  idum2=idum
  do j=NTAB+8,1,-1
   k=idum/IQ1
   idum=IA1*(idum-k*IQ1)-k*IR1
   if (idum.lt.0) idum=idum+IM1
   if (j.le.NTAB) iv(j)=idum
  end do
  iy=iv(1)
 endif
 k=idum/IQ1
 idum=IA1*(idum-k*IQ1)-k*IR1
 if (idum.lt.0) idum=idum+IM1
 k=idum2/IQ2
 idum2=IA2*(idum2-k*IQ2)-k*IR2
 if (idum2.lt.0) idum2=idum2+IM2
 j=1+iy/NDIV
 iy=iv(j)-idum2
 iv(j)=idum
 if(iy.lt.1)iy=iy+IMM1
 ran2 = min(AM*iy,RNMX)

 ! save new walues to the class variables
 ranum%iy = iy
 ranum%idum2 = idum2
 ranum%iv(:) = iv(:)

 RETURN 

END FUNCTION 

