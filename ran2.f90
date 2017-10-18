DOUBLE PRECISION   FUNCTION  ran2(idum)

USE ran2_class

 IMPLICIT NONE 
 CLASS(ranvar), POINTER      :: this

 INTEGER PARAMETER           :: IM1 = 2147483563, IM2 = 2147483399, IMM1 = IM1 - 1, &
                                IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791, &
                                NDIV=1+IMM1/NTAB
 INTEGER                     :: j,k
 INTEGER                     :: NTAB, iy, idum2
 INTEGER, DIMENSION(NTAB)    :: iv

 IF(.NOT. ASSOCIATED(this)) this = initialize_values()

 NTAB = this%NTAB
 iy = this%iy
 idum2 = this%idum2
 iv(:) = this%iv(:)

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
 ran2=min(AM*iy,RNMX)

 RETURN

END FUNCTION 

