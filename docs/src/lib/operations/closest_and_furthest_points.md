# [Closest and Furthest Points](@id closest_furthest_points)

```@meta
CurrentModule = SearchableGeometries.GeometricPrimitives
```

Closest-point and furthest-point queries identify representative points in a geometric object relative to a query point or query geometry.

## Generic functions

```@docs
get_closest_point
get_furthest_point
```

## Bounding volume and point

```@docs; canonical=false
get_closest_point(::BoundingVolume, ::Vector{<:Real})
get_furthest_point(::BoundingVolume, ::Vector{<:Real})
```

## Point and hyperplane

```@docs; canonical=false
get_closest_point(::Vector{<:Real}, ::Hyperplane)
```

## Bounding volume and hyperplane

```@docs; canonical=false
get_closest_point(::BoundingVolume, ::Hyperplane)
get_furthest_point(::BoundingVolume, ::Hyperplane)
```

## Linear extrema helpers

```@docs; canonical=false
get_extreme_point
get_antiextreme_point
```