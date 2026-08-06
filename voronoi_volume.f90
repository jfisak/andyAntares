! Estimate Voronoi cell volumes from the current model grid.
!
! High-level method:
! 1) Draw many deterministic quasi-random samples over a padded domain.
! 2) Assign each sample to the nearest model point (Voronoi ownership).
! 3) Convert the ownership fraction into physical volume.
!
! Geometry handling by model type:
! - model_type == 3: full 3D Cartesian Voronoi (x, y, z from vec_pos).
! - model_type == 2:
!   * inputmodel /= 2: axisymmetric case using (radius, theta), then revolve
!     around the symmetry axis via Jacobian weight 2*pi*r^2*sin(theta).
!   * inputmodel == 2: axisymmetric case using cylindrical meridional
!     coordinates (R,z), then revolve with Jacobian weight 2*pi*R.
SUBROUTINE voronoi_volume()

USE types, ONLY : const_dimofspace, ind_x, ind_y, ind_z, model_grid, model_type, n_modelgrid, inputmodel
USE constants, ONLY : const_pi

IMPLICIT NONE

DOUBLE PRECISION                                     :: volumes(n_modelgrid)

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
IF(n_modelgrid <= 0) RETURN

! Guard: model grid data must be available before sampling.
IF(.NOT. ALLOCATED(model_grid)) STOP 'voronoi_volume: model_grid is not allocated'

! Guard: requested number of points cannot exceed available model entries.
IF(SIZE(model_grid) < n_modelgrid) STOP 'voronoi_volume: n_modelgrid exceeds size(model_grid)'

! Guard: this routine currently supports only 2D axisymmetric and 3D models.
IF(model_type /= 2 .AND. model_type /= 3) STOP 'voronoi_volume: only model_type 2 or 3 is supported'

! Explicit initialization of output array.
volumes(:) = 0.D0

! Sampling budget:
! - keep a robust minimum for small models,
! - scale linearly with number of points for larger models.
n_samples = MAX(n_samples_min, 100*n_modelgrid)

IF(model_type == 3) THEN

    ! Integer ownership counter for each site in Cartesian 3D.
    ALLOCATE(hit_count(n_modelgrid))
    hit_count(:) = 0

    ! Build a padded 3D bounding box in Cartesian coordinates (x, y, z).
    ! The padding reduces edge truncation artifacts in finite sampling.
    DO dim_idx = 1, 3
     box_min(dim_idx) = model_grid(1)%vec_pos(dim_idx)
     box_max(dim_idx) = model_grid(1)%vec_pos(dim_idx)

     ! Find coordinate extrema over all model points for this dimension.
     DO point_idx = 2, n_modelgrid
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
     DO point_idx = 1, n_modelgrid
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
    DO cur_mgi = 1, n_modelgrid
     volumes(cur_mgi) = box_volume*DBLE(hit_count(cur_mgi))/DBLE(n_samples)
     model_grid(cur_mgi)%voronoi_volume = volumes(cur_mgi)
    END DO

    ! Cleanup temporary ownership array.
    DEALLOCATE(hit_count)

ELSE IF(model_type == 2) THEN

    ! Weighted ownership accumulator for axisymmetric geometry.
    ALLOCATE(hit_weight(n_modelgrid))
    hit_weight(:) = 0.D0

    IF(inputmodel == 2) THEN

     ! Peku-like geometry: sample directly in (R,z) and revolve around z-axis.
     box_min(ind_x) = MINVAL(model_grid(1:n_modelgrid)%rxywind)
     box_max(ind_x) = MAXVAL(model_grid(1:n_modelgrid)%rxywind)
     box_min(ind_y) = MINVAL(model_grid(1:n_modelgrid)%zwind)
     box_max(ind_y) = MAXVAL(model_grid(1:n_modelgrid)%zwind)

     DO dim_idx = 1, 2
      box_width(dim_idx) = box_max(dim_idx) - box_min(dim_idx)
      IF(box_width(dim_idx) <= 0.D0) box_width(dim_idx) = 1.D0
      box_min(dim_idx) = box_min(dim_idx) - padding_frac*box_width(dim_idx)
      box_max(dim_idx) = box_max(dim_idx) + padding_frac*box_width(dim_idx)
      box_width(dim_idx) = box_max(dim_idx) - box_min(dim_idx)
     END DO

     ! Cylindrical radius must be non-negative.
     box_min(ind_x) = MAX(0.D0, box_min(ind_x))
     box_width(ind_x) = box_max(ind_x) - box_min(ind_x)
     box_width(ind_y) = box_max(ind_y) - box_min(ind_y)
     IF(box_width(ind_x) <= 0.D0) box_width(ind_x) = 1.D0
     IF(box_width(ind_y) <= 0.D0) box_width(ind_y) = 1.D0

     box_area = box_width(ind_x)*box_width(ind_y)

     DO sample_idx = 1, n_samples
      sample_mer_r = box_min(ind_x) + box_width(ind_x)*halton_coord(sample_idx, 2)
      sample_mer_z = box_min(ind_y) + box_width(ind_y)*halton_coord(sample_idx, 3)

      nearest_idx = 1
      best_dist2 = HUGE(1.D0)

      DO point_idx = 1, n_modelgrid
       point_mer_r = model_grid(point_idx)%rxywind
       point_mer_z = model_grid(point_idx)%zwind
       dist2 = (sample_mer_r - point_mer_r)**2 + &
           & (sample_mer_z - point_mer_z)**2
       IF(dist2 < best_dist2) THEN
        best_dist2 = dist2
        nearest_idx = point_idx
       END IF
      END DO

      ! dV = 2*pi*R*dR*dz
      weight = 2.D0*const_pi*sample_mer_r
      hit_weight(nearest_idx) = hit_weight(nearest_idx) + weight
     END DO

    ELSE

     ! Original axisymmetric geometry using (radius, theta).
     box_min(ind_x) = MINVAL(model_grid(1:n_modelgrid)%rwind)
     box_max(ind_x) = MAXVAL(model_grid(1:n_modelgrid)%rwind)
     box_min(ind_y) = MINVAL(model_grid(1:n_modelgrid)%angle)
     box_max(ind_y) = MAXVAL(model_grid(1:n_modelgrid)%angle)

     DO dim_idx = 1, 2
      DO cur_mgi = 2, n_modelgrid
       IF(dim_idx == ind_x) THEN
        box_min(dim_idx) = MIN(box_min(dim_idx), model_grid(cur_mgi)%rwind)
        box_max(dim_idx) = MAX(box_max(dim_idx), model_grid(cur_mgi)%rwind)
       ELSE
        box_min(dim_idx) = MIN(box_min(dim_idx), model_grid(cur_mgi)%angle)
        box_max(dim_idx) = MAX(box_max(dim_idx), model_grid(cur_mgi)%angle)
       END IF
       write(*,*) 'voronoi_volume: cur_mgi = ', cur_mgi, ' box_min(', dim_idx, ') = ', box_min(dim_idx), &
           & ' box_max(', dim_idx, ') = ', box_max(dim_idx)
      END DO

      box_width(dim_idx) = box_max(dim_idx) - box_min(dim_idx)
      IF(box_width(dim_idx) <= 0.D0) box_width(dim_idx) = 1.D0
      box_min(dim_idx) = box_min(dim_idx) - padding_frac*box_width(dim_idx)
      box_max(dim_idx) = box_max(dim_idx) + padding_frac*box_width(dim_idx)
      box_width(dim_idx) = box_max(dim_idx) - box_min(dim_idx)
     END DO

     ! Clamp domain to physical ranges: radius >= 0, theta in [0, pi].
     box_min(ind_x) = MAX(0.D0, box_min(ind_x))
     box_min(ind_y) = MAX(0.D0, box_min(ind_y))
     box_max(ind_y) = MIN(const_pi, box_max(ind_y))
     box_width(ind_x) = box_max(ind_x) - box_min(ind_x)
     box_width(ind_y) = box_max(ind_y) - box_min(ind_y)
     IF(box_width(ind_x) <= 0.D0) box_width(ind_x) = 1.D0
     IF(box_width(ind_y) <= 0.D0) box_width(ind_y) = 1.D0

     box_area = box_width(ind_x)*box_width(ind_y)

     DO sample_idx = 1, n_samples
      r_sample = box_min(ind_x) + box_width(ind_x)*halton_coord(sample_idx, 2)
      theta_sample = box_min(ind_y) + box_width(ind_y)*halton_coord(sample_idx, 3)

      sample_mer_r = r_sample*SIN(theta_sample)
      sample_mer_z = r_sample*COS(theta_sample)

      nearest_idx = 1
      best_dist2 = HUGE(1.D0)

      DO cur_mgi = 1, n_modelgrid
       point_mer_r = model_grid(cur_mgi)%rwind*SIN(model_grid(cur_mgi)%angle)
       point_mer_z = model_grid(cur_mgi)%rwind*COS(model_grid(cur_mgi)%angle)
       write(*,*) 'voronoi_volume: cur_mgi = ', cur_mgi, ' point_mer_r = ', point_mer_r, &
           & ' point_mer_z = ', point_mer_z
       dist2 = (sample_mer_r - point_mer_r)**2 + &
           & (sample_mer_z - point_mer_z)**2
       IF(dist2 < best_dist2) THEN
        best_dist2 = dist2
        nearest_idx = point_idx
       END IF
      END DO

      weight = 2.D0*const_pi*r_sample**2*SIN(theta_sample)
      hit_weight(nearest_idx) = hit_weight(nearest_idx) + weight
     END DO

    END IF

    ! Weighted Monte Carlo average converted to physical 3D volume estimate.
    ! volume_i = box_area * (sum(weights for site i) / total_samples)
    DO cur_mgi = 1, n_modelgrid
     volumes(cur_mgi) = box_area*hit_weight(cur_mgi)/DBLE(n_samples)
     write(*,*) 'voronoi_volume: mgi = ', cur_mgi, ' volume = ', volumes(cur_mgi)
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