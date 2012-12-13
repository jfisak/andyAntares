SUBROUTINE init_photonpack(n_pack)

  USE types

  IMPLICIT NONE

    INTEGER                           :: I, n_pack, ind_cell_numb, ind_x, ind_y, ind_z
    DOUBLE PRECISION                  :: L_star
    DOUBLE PRECISION, PARAMETER       :: delta_t=1.D0
    DOUBLE PRECISION, DIMENSION(3)    :: direction

    L_star = 4.D0*pi*(R_star*r_sun)**2*sigma*T_eff
    R_star = 0.

!    ind_x = nx_cell/2 + 1
!    ind_y = ny_cell/2 + 1
!    ind_z = nz_cell/2 + 1
!    ind_cell_numb = (ind_x-1)*ny_cell*nz_cell + (ind_y-1)*nz_cell + ind_z
!    print*, ind_cell_numb

    DO I=1,n_pack
       package(I)%pos = 10.D-5
       CALL random_unitvector(direction)
       package(I)%dir = direction
       package(I)%e_rf = (L_star/n_pack) * delta_t  
!      Only for regular grid!!!!
       ind_x = FLOOR(package(I)%pos(1)/cell_width + DBLE(nx_cell)/2) + 1
       ind_y = FLOOR(package(I)%pos(2)/cell_width + DBLE(ny_cell)/2) + 1
       ind_z = FLOOR(package(I)%pos(3)/cell_width + DBLE(nz_cell)/2) + 1
       ind_cell_numb = (ind_x-1)*ny_cell*nz_cell + (ind_y-1)*nz_cell + ind_z
       IF ((ind_cell_numb .GT. nx_cell*ny_cell*nz_cell) .OR. (ind_cell_numb .LT. 1)) &
           STOP 'ERROR in cell_number'
       package(I)%cell_numb = ind_cell_numb
       package(I)%active = 1
       package(I)%typ = type_rpkt
       package(I)%last_cross = NONE
!       print*, package(I)%cell_numb,package(I)%dir !,  package(I)%pos, package(I)% e_rf
!      e_cmf, freq_cmf, freq_rf, cell_numb, pack_numb, active
    END DO

!    PRINT*, ind_x, ind_y, ind_z, ind_cell_numb

!    I=10
!    package(I)%pos = 95.
!    package(I)%dir = -1.
!    package(I)%e_rf = 0.

END 
