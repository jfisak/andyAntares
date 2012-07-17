SUBROUTINE do_event(pack_index, tau, tau_rand)

  USE types

  IMPLICIT NONE    

    INTEGER                           :: pack_index
    DOUBLE PRECISION                  :: dist, rand_numb, ran2, tau, tau_rand
    DOUBLE PRECISION, DIMENSION(3)    :: direction


!    print*, rand_numb
    rand_numb = 0.1D0
!    rand_numb = ran2(idum)     
    IF (rand_numb .LT. 0.5D0) THEN
        package(pack_index)%active = 0
!!        print*, 'DO ABSORPTION'
    ELSE
        CALL emit_rpackage(pack_index)
        tau = 0.D0
10      rand_numb = ran2(idum)     
        IF (rand_numb .EQ. 0.D0) GOTO 10    
        tau_rand = -LOG(rand_numb)
!!        print*, 'DO SCATTERING'
    END IF
    
END

! For absorption only  

!     package(pack_index)%active = 0

!END
