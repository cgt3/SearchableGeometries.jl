# Balls

```@meta
CurrentModule = SearchableGeometries
```

A ball is a distance-based geometric object. It represents all points within a fixed radius of a center point.

In SearchableGeometries.jl, balls are defined using a ``p``-norm. A full-dimensional ball with center `c`, radius `r`, and norm parameter `p` represents

```math
\{x : \|x - c\|_p \le r\}.
```

Balls are useful in search because many search rules are distance-based. For example, an epsilon-ball query asks for all points within a fixed radius of a query point.

## Type

```@docs
Ball
```

## Bounding-volume enclosure

```@docs
BoundingVolume(::Ball)
```

## Containment

```@docs
is_contained(::Ball, ::Vector{<:Real})
is_contained(::BoundingVolume, ::Ball)
is_contained(::Ball, ::BoundingVolume)
```

## Intersection with bounding volumes

```@docs
intersects(::BoundingVolume, ::Ball)
get_intersection(::BoundingVolume, ::Ball)
```

## Dimension reduction and tightening

```@docs
get_reduced_dim_ball
tighten_bv_bounds!
```