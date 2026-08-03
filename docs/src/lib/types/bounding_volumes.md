# [BoundingVolume](@id def_BoundingVolume)

```@meta
CurrentModule = SearchableGeometries.GeometricPrimitives
```

```@docs
BoundingVolume
```

## Constructors

```julia
BoundingVolume()
BoundingVolume(lb, ub; tol=DEFAULT_BV_POINT_TOL)
```

## Operations

```@docs
get_closest_point(::BoundingVolume, ::Vector{<:Real})
get_furthest_point(::BoundingVolume, ::Vector{<:Real})
is_contained(::BoundingVolume, ::Vector{<:Real})
is_contained(::BoundingVolume, ::BoundingVolume)
intersects(::BoundingVolume, ::BoundingVolume)
get_intersection(::BoundingVolume, ::BoundingVolume)
```

## Faces

```@docs
face_index_to_spatial_index
get_face_bounding_volume
```


