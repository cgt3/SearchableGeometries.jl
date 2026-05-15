## Types

```@docs
SearchableGeometry
BoundingVolume
Ball
Hyperplane
```

## Generic operations

```@docs
is_contained
intersects
get_intersection
get_closest_point
get_furthest_point
```

## Bounding-volume helpers

```@docs
face_index_to_spatial_index
get_face_bounding_volume
```

## Ball helpers

```@docs
get_reduced_dim_ball
tighten_bv_bounds!
```

## Hyperplane helpers

```@docs
get_furthest_pt
get_antifurthest_pt
tighten_bv_bounds
```