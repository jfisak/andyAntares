SUBROUTINE  ran2(idum)

 USE ran2_class
 IMPLICIT NONE 


 INTEGER, PARAMETER          :: IM1 = 2147483563, IM2 = 2147483399, IMM1 = IM1 - 1, &
                                IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791, &
                                NDIV=1+IMM1/NTAB
 INTEGER                     :: j,k
 INTEGER                     :: NTAB, iy, idum2
 INTEGER, DIMENSION(NTAB)    :: iv
 CLASS(ranvar)      :: ranum
 DOUBLE PRECISION            :: out_ran

 IF(.NOT. ASSOCIATED(ranum)) ranum%init()

 NTAB = ranum%NTAB
 iy = ranum%iy
 idum2 = ranum%idum2
 iv(:) = ranum%iv(:)

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
 out_ran=min(AM*iy,RNMX)

 RETURN out_ran

END SUBROUTINE 

