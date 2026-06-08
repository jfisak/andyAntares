! Estimate Voronoi cell volumes in 3D for a set of sites.
! The volumes are computed inside an enclosing bounding box around the points.
SUBROUTINE voronoi_volume(n_points, positions, volumes)

USE types, ONLY : const_dimofspace, ind_x, ind_y, ind_z
USE constants, ONLY : const_pi

IMPLICIT NONE

INTEGER, INTENT(IN)                                  :: n_points
DOUBLE PRECISION, INTENT(IN)                         :: positions(:,:)
DOUBLE PRECISION, INTENT(OUT)                        :: volumes(n_points)

INTEGER, PARAMETER                                   :: n_samples_min = 30000
DOUBLE PRECISION, PARAMETER                          :: padding_frac = 5.D-2

INTEGER                                              :: n_samples
INTEGER                                              :: sample_idx, point_idx, dim_idx, nearest_idx
INTEGER                                              :: n_dim
DOUBLE PRECISION                                     :: dist2, best_dist2, box_volume, box_area
DOUBLE PRECISION                                     :: r_sample, z_sample, weight
DOUBLE PRECISION, DIMENSION(const_dimofspace)        :: box_min, box_max, box_width, cur_pos
INTEGER, ALLOCATABLE                                 :: hit_count(:)
DOUBLE PRECISION, ALLOCATABLE                        :: hit_weight(:)

IF(n_points <= 0) RETURN
IF(SIZE(positions, 1) /= n_points) STOP 'voronoi_volume: positions first dimension must be n_points'

n_dim = SIZE(positions, 2)
IF(n_dim /= 2 .AND. n_dim /= 3) STOP 'voronoi_volume: positions must have 2 (r,z) or 3 (x,y,z) columns'

volumes(:) = 0.D0

n_samples = MAX(n_samples_min, 100*n_points)

IF(n_dim == 3) THEN

    ALLOCATE(hit_count(n_points))
    hit_count(:) = 0.D0

    ! Build a padded 3D box in Cartesian coordinates x, y, z.
    DO dim_idx = 1, 3
     box_min(dim_idx) = MINVAL(positions(:, dim_idx))
     box_max(dim_idx) = MAXVAL(positions(:, dim_idx))
     box_width(dim_idx) = box_max(dim_idx) - box_min(dim_idx)

     IF(box_width(dim_idx) <= 0.D0) THEN
      box_width(dim_idx) = 1.D0
     END IF

     box_min(dim_idx) = box_min(dim_idx) - padding_frac*box_width(dim_idx)
     box_max(dim_idx) = box_max(dim_idx) + padding_frac*box_width(dim_idx)
     box_width(dim_idx) = box_max(dim_idx) - box_min(dim_idx)
    END DO

    box_volume = box_width(ind_x)*box_width(ind_y)*box_width(ind_z)

    ! Sample the box with a deterministic low-discrepancy sequence.
    DO sample_idx = 1, n_samples
     cur_pos(ind_x) = box_min(ind_x) + box_width(ind_x)*halton_coord(sample_idx, 2)
     cur_pos(ind_y) = box_min(ind_y) + box_width(ind_y)*halton_coord(sample_idx, 3)
     cur_pos(ind_z) = box_min(ind_z) + box_width(ind_z)*halton_coord(sample_idx, 5)

     nearest_idx = 1
     best_dist2 = HUGE(1.D0)

     ! Find nearest site in Euclidean metric.
     DO point_idx = 1, n_points
      dist2 = (cur_pos(ind_x) - positions(point_idx, ind_x))**2 + &
          & (cur_pos(ind_y) - positions(point_idx, ind_y))**2 + &
          & (cur_pos(ind_z) - positions(point_idx, ind_z))**2

      IF(dist2 < best_dist2) THEN
       best_dist2 = dist2
       nearest_idx = point_idx
      END IF
     END DO

     hit_count(nearest_idx) = hit_count(nearest_idx) + 1
    END DO

    ! Convert occupancy fractions to physical volumes.
    DO point_idx = 1, n_points
     volumes(point_idx) = box_volume*DBLE(hit_count(point_idx))/DBLE(n_samples)
    END DO

    DEALLOCATE(hit_count)

ELSE IF(n_dim == 2) THEN

    ALLOCATE(hit_weight(n_points))
    hit_weight(:) = 0.D0

    ! Build a padded 2D box in axisymmetric coordinates r, z.
    DO dim_idx = 1, 2
     box_min(dim_idx) = MINVAL(positions(:, dim_idx))
     box_max(dim_idx) = MAXVAL(positions(:, dim_idx))
     box_width(dim_idx) = box_max(dim_idx) - box_min(dim_idx)

     IF(box_width(dim_idx) <= 0.D0) THEN
      box_width(dim_idx) = 1.D0
     END IF

     box_min(dim_idx) = box_min(dim_idx) - padding_frac*box_width(dim_idx)
     box_max(dim_idx) = box_max(dim_idx) + padding_frac*box_width(dim_idx)
     box_width(dim_idx) = box_max(dim_idx) - box_min(dim_idx)
    END DO

    ! For axisymmetric 2D input (r,z), avoid negative radii after padding.
    box_min(ind_x) = MAX(0.D0, box_min(ind_x))
    box_width(ind_x) = box_max(ind_x) - box_min(ind_x)
    IF(box_width(ind_x) <= 0.D0) box_width(ind_x) = 1.D0

    box_area = box_width(ind_x)*box_width(ind_y)

    ! Sample in (r,z); each sample represents volume element 2*pi*r*dr*dz.
    DO sample_idx = 1, n_samples
     r_sample = box_min(ind_x) + box_width(ind_x)*halton_coord(sample_idx, 2)
     z_sample = box_min(ind_y) + box_width(ind_y)*halton_coord(sample_idx, 3)

     nearest_idx = 1
     best_dist2 = HUGE(1.D0)

     ! Distance in axisymmetric geometry equals distance in the meridional plane.
     DO point_idx = 1, n_points
      dist2 = (r_sample - positions(point_idx, ind_x))**2 + &
          & (z_sample - positions(point_idx, ind_y))**2

      IF(dist2 < best_dist2) THEN
       best_dist2 = dist2
       nearest_idx = point_idx
      END IF
     END DO

     weight = 2.D0*const_pi*r_sample
     hit_weight(nearest_idx) = hit_weight(nearest_idx) + weight
    END DO

    ! Weighted average in (r,z) gives 3D rotated Voronoi volumes.
    DO point_idx = 1, n_points
     volumes(point_idx) = box_area*hit_weight(point_idx)/DBLE(n_samples)
    END DO

    DEALLOCATE(hit_weight)

END IF

CONTAINS

 DOUBLE PRECISION FUNCTION halton_coord(index_in, base)

 IMPLICIT NONE

 INTEGER, INTENT(IN)                           :: index_in, base
 INTEGER                                       :: idx
 DOUBLE PRECISION                              :: factor

 halton_coord = 0.D0
 idx = index_in
 factor = 1.D0/DBLE(base)

 ! Radical-inverse expansion in the requested base.
 DO WHILE(idx > 0)
  halton_coord = halton_coord + factor*DBLE(MOD(idx, base))
  idx = idx/base
  factor = factor/DBLE(base)
 END DO

 END FUNCTION halton_coord

END SUBROUTINE voronoi_volume