# Closest and Furthest Points

```@meta
CurrentModule = SearchableGeometries
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
get_furthest_pt
get_antifurthest_pt
```