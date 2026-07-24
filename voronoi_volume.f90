! Estimate Voronoi cell volumes from the current model grid.
!
! High-level method:
! 1) Draw many deterministic quasi-random samples over a padded domain.
! 2) Assign each sample to the nearest model point (Voronoi ownership).
! 3) Convert the ownership fraction into physical volume.
!
! Geometry handling by model type:
! - model_type == 3: full 3D Cartesian Voronoi (x, y, z from vec_pos).
! - model_type == 2: axisymmetric case using (radius, theta), then revolve
!   around the symmetry axis via Jacobian weight 2*pi*r^2*sin(theta).
SUBROUTINE voronoi_volume(n_points, volumes)

USE types, ONLY : const_dimofspace, ind_x, ind_y, ind_z, model_grid, model_type
USE constants, ONLY : const_pi

IMPLICIT NONE

INTEGER, INTENT(IN)                                  :: n_points
DOUBLE PRECISION, INTENT(OUT)                        :: volumes(n_points)

INTEGER, PARAMETER                                   :: n_samples_min = 30000
DOUBLE PRECISION, PARAMETER                          :: padding_frac = 5.D-2

INTEGER                                              :: n_samples
INTEGER                                              :: sample_idx, point_idx, dim_idx, nearest_idx
integer                                              :: cur_mgi
DOUBLE PRECISION                                     :: dist2, best_dist2, box_volume, box_area
DOUBLE PRECISION                                     :: r_sample, theta_sample, weight
DOUBLE PRECISION                                     :: sample_mer_r, sample_mer_z
DOUBLE PRECISION                                     :: point_mer_r, point_mer_z
DOUBLE PRECISION, DIMENSION(const_dimofspace)        :: box_min, box_max, box_width, cur_pos
INTEGER, ALLOCATABLE                                 :: hit_count(:)
DOUBLE PRECISION, ALLOCATABLE                        :: hit_weight(:)

! Guard: nothing to do.
IF(n_points <= 0) RETURN

! Guard: model grid data must be available before sampling.
IF(.NOT. ALLOCATED(model_grid)) STOP 'voronoi_volume: model_grid is not allocated'

! Guard: requested number of points cannot exceed available model entries.
IF(SIZE(model_grid) < n_points) STOP 'voronoi_volume: n_points exceeds size(model_grid)'

! Guard: this routine currently supports only 2D axisymmetric and 3D models.
IF(model_type /= 2 .AND. model_type /= 3) STOP 'voronoi_volume: only model_type 2 or 3 is supported'

! Explicit initialization of output array.
volumes(:) = 0.D0

! Sampling budget:
! - keep a robust minimum for small models,
! - scale linearly with number of points for larger models.
n_samples = MAX(n_samples_min, 100*n_points)

IF(model_type == 3) THEN

    ! Integer ownership counter for each site in Cartesian 3D.
    ALLOCATE(hit_count(n_points))
    hit_count(:) = 0

    ! Build a padded 3D bounding box in Cartesian coordinates (x, y, z).
    ! The padding reduces edge truncation artifacts in finite sampling.
    DO dim_idx = 1, 3
     box_min(dim_idx) = model_grid(1)%vec_pos(dim_idx)
     box_max(dim_idx) = model_grid(1)%vec_pos(dim_idx)

     ! Find coordinate extrema over all model points for this dimension.
     DO point_idx = 2, n_points
      box_min(dim_idx) = MIN(box_min(dim_idx), model_grid(point_idx)%vec_pos(dim_idx))
      box_max(dim_idx) = MAX(box_max(dim_idx), model_grid(point_idx)%vec_pos(dim_idx))
     END DO

     box_width(dim_idx) = box_max(dim_idx) - box_min(dim_idx)

     ! Degenerate axis protection (all points identical in this coordinate).
     IF(box_width(dim_idx) <= 0.D0) THEN
      box_width(dim_idx) = 1.D0
     END IF

     ! Apply symmetric padding and recompute width.
     box_min(dim_idx) = box_min(dim_idx) - padding_frac*box_width(dim_idx)
     box_max(dim_idx) = box_max(dim_idx) + padding_frac*box_width(dim_idx)
     box_width(dim_idx) = box_max(dim_idx) - box_min(dim_idx)
    END DO

    ! Physical volume represented by the sampled Cartesian box.
    box_volume = box_width(ind_x)*box_width(ind_y)*box_width(ind_z)

    ! Sample with Halton low-discrepancy sequence (bases 2,3,5).
    ! This gives more uniform coverage than pseudo-random points at equal N.
    DO sample_idx = 1, n_samples
     cur_pos(ind_x) = box_min(ind_x) + box_width(ind_x)*halton_coord(sample_idx, 2)
     cur_pos(ind_y) = box_min(ind_y) + box_width(ind_y)*halton_coord(sample_idx, 3)
     cur_pos(ind_z) = box_min(ind_z) + box_width(ind_z)*halton_coord(sample_idx, 5)

     ! Track nearest Voronoi site for this sample.
     nearest_idx = 1
     best_dist2 = HUGE(1.D0)

     ! Nearest-neighbor query in standard Euclidean metric.
     DO point_idx = 1, n_points
      dist2 = (cur_pos(ind_x) - model_grid(point_idx)%vec_pos(ind_x))**2 + &
          & (cur_pos(ind_y) - model_grid(point_idx)%vec_pos(ind_y))**2 + &
          & (cur_pos(ind_z) - model_grid(point_idx)%vec_pos(ind_z))**2

      ! Keep the current best candidate.
      IF(dist2 < best_dist2) THEN
       best_dist2 = dist2
       nearest_idx = point_idx
      END IF
     END DO

     ! Increase ownership count for the winning site.
     hit_count(nearest_idx) = hit_count(nearest_idx) + 1
    END DO

    ! Convert occupancy fraction to Voronoi volume estimate.
    ! volume_i = box_volume * (hits_i / total_samples)
    DO cur_mgi = 1, n_points
     volumes(cur_mgi) = box_volume*DBLE(hit_count(cur_mgi))/DBLE(n_samples)
     model_grid(cur_mgi)%voronoi_volume = volumes(cur_mgi)
    END DO

    ! Cleanup temporary ownership array.
    DEALLOCATE(hit_count)

ELSE IF(model_type == 2) THEN

    ! Weighted ownership accumulator for axisymmetric geometry.
    ALLOCATE(hit_weight(n_points))
    hit_weight(:) = 0.D0

    ! Build a padded 2D box in (radius, theta).
    ! This is the sampling domain before applying axisymmetric weights.
    box_min(ind_x) = MINVAL(model_grid(:)%rwind)
    box_max(ind_x) = MAXVAL(model_grid(:)%rwind)
    box_min(ind_y) = MINVAL(model_grid(:)%angle)
    box_max(ind_y) = MAXVAL(model_grid(:)%angle)

    DO dim_idx = 1, 2
     ! Re-evaluate extrema explicitly over the requested point range.
     DO point_idx = 2, n_points
      IF(dim_idx == ind_x) THEN
       box_min(dim_idx) = MIN(box_min(dim_idx), model_grid(point_idx)%rwind)
       box_max(dim_idx) = MAX(box_max(dim_idx), model_grid(point_idx)%rwind)
      ELSE
       box_min(dim_idx) = MIN(box_min(dim_idx), model_grid(point_idx)%angle)
       box_max(dim_idx) = MAX(box_max(dim_idx), model_grid(point_idx)%angle)
      END IF
     END DO

     box_width(dim_idx) = box_max(dim_idx) - box_min(dim_idx)

     ! Degenerate range protection in either radius or theta dimension.
     IF(box_width(dim_idx) <= 0.D0) THEN
      box_width(dim_idx) = 1.D0
     END IF

     ! Apply symmetric padding around the data extent.
     box_min(dim_idx) = box_min(dim_idx) - padding_frac*box_width(dim_idx)
     box_max(dim_idx) = box_max(dim_idx) + padding_frac*box_width(dim_idx)
     box_width(dim_idx) = box_max(dim_idx) - box_min(dim_idx)
    END DO

    ! Clamp domain to physical ranges:
    ! radius >= 0, theta in [0, pi].
    box_min(ind_x) = MAX(0.D0, box_min(ind_x))
    box_min(ind_y) = MAX(0.D0, box_min(ind_y))
    box_max(ind_y) = MIN(const_pi, box_max(ind_y))
    box_width(ind_x) = box_max(ind_x) - box_min(ind_x)
    box_width(ind_y) = box_max(ind_y) - box_min(ind_y)
    ! Final safeguard against zero-width sampling interval.
    IF(box_width(ind_x) <= 0.D0) box_width(ind_x) = 1.D0
    IF(box_width(ind_y) <= 0.D0) box_width(ind_y) = 1.D0

    ! Area of the sampled (r, theta) rectangle.
    box_area = box_width(ind_x)*box_width(ind_y)

    ! Sample in (radius, theta) using Halton bases 2 and 3.
    ! Each sample carries Jacobian weight:
    ! dV = 2*pi*r^2*sin(theta) * dr * dtheta,
    ! where 2*pi comes from revolution around the symmetry axis.
    DO sample_idx = 1, n_samples
     r_sample = box_min(ind_x) + box_width(ind_x)*halton_coord(sample_idx, 2)
     theta_sample = box_min(ind_y) + box_width(ind_y)*halton_coord(sample_idx, 3)

     ! Map sample from spherical (r, theta) to meridional (R,z) plane.
     sample_mer_r = r_sample*SIN(theta_sample)
     sample_mer_z = r_sample*COS(theta_sample)

     ! Track nearest Voronoi site for this sample.
     nearest_idx = 1
     best_dist2 = HUGE(1.D0)

     ! Compare distances in meridional plane coordinates.
     ! This defines the Voronoi partition used in axisymmetric mode.
     DO point_idx = 1, n_points
      point_mer_r = model_grid(point_idx)%rwind*SIN(model_grid(point_idx)%angle)
      point_mer_z = model_grid(point_idx)%rwind*COS(model_grid(point_idx)%angle)
      dist2 = (sample_mer_r - point_mer_r)**2 + &
          & (sample_mer_z - point_mer_z)**2

      ! Keep the nearest candidate.
      IF(dist2 < best_dist2) THEN
       best_dist2 = dist2
       nearest_idx = point_idx
      END IF
     END DO

     ! Jacobian-based volume contribution for this sampled direction/radius.
     weight = 2.D0*const_pi*r_sample**2*SIN(theta_sample)

     ! Add weighted ownership to the winning site.
     hit_weight(nearest_idx) = hit_weight(nearest_idx) + weight
    END DO

    ! Weighted Monte Carlo average converted to physical 3D volume estimate.
    ! volume_i = box_area * (sum(weights for site i) / total_samples)
    DO cur_mgi = 1, n_points
     volumes(cur_mgi) = box_area*hit_weight(cur_mgi)/DBLE(n_samples)
     model_grid(cur_mgi)%voronoi_volume = volumes(cur_mgi)
    END DO

    ! Cleanup temporary weighted ownership array.
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

 ! Halton coordinate via radical-inverse expansion in the chosen base:
 ! index digits are reflected across the decimal point in that base.
 DO WHILE(idx > 0)
  halton_coord = halton_coord + factor*DBLE(MOD(idx, base))
  idx = idx/base
  factor = factor/DBLE(base)
 END DO

 END FUNCTION halton_coord

END SUBROUTINE voronoi_volume