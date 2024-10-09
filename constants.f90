MODULE constants

IMPLICIT NONE

!! Physical constants
  DOUBLE PRECISION, PARAMETER        :: const_pi=3.1415926535897932D+00,&
                                        const_me_g=9.1093837015D-28,&
                                        const_mp_g=1.67262192369D-24,&
                                        const_sigma_e=6.6524587321D-25,&
                                        const_h=6.62607015D-27,&
                                        const_c=2.99792458D+10,&
                                        const_e=4.8032068D-10,&
                                        const_rsun=6957.D+07,&
                                        const_kB=1.380649D-16,&
                                        ! m_sun=1.988409870698051D+33,&
                                        const_stefbolz =5.67037442D-05 !ergcm^(-2)s(-1)K(-4) !D. H. Cohen et al.2012
  DOUBLE PRECISION, PARAMETER        :: const_pc=30.85677814913674D17,&
                                        const_ev = 1.60217646D-12
                                        ! saha_const=2.0706839D-16,&
                                        ! saha_const=4.1414D-16,&


  DOUBLE PRECISION                      :: saha_const














END MODULE constants
