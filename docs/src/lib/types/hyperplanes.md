# [Hyperplane](@id def_Hyperplane)

```@meta
CurrentModule = SearchableGeometries.GeometricPrimitives
```

```@docs
Hyperplane
```

## Constructors

```julia
Hyperplane(point, n)
```

## Point operations

```@docs
is_contained(::Hyperplane, ::Vector{<:Real})
get_closest_point(::Vector{<:Real}, ::Hyperplane)
```

## Bounding-volume operations

```@docs
intersects(::BoundingVolume, ::Hyperplane)
get_intersection(::BoundingVolume, ::Hyperplane)
get_closest_point(::BoundingVolume, ::Hyperplane)
get_furthest_point(::BoundingVolume, ::Hyperplane)
tighten_bv_bounds
```

## Linear extrema

```@docs
get_furthest_pt
get_antifurthest_pt
```


