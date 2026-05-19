# [Balls](@id balls_manual)

```@meta
CurrentModule = SearchableGeometries
```

A ball is a distance-based region used to describe all points within a fixed radius of a center point.

In two dimensions, a ball may look like a disk when using the Euclidean norm. In three dimensions, it may look like a solid sphere. In higher dimensions, it is a higher-dimensional distance region.

The exact shape of the ball depends on the chosen ``p``-norm.

## Mathematical definition

A ball with center `c`, radius `r`, and norm parameter `p` represents

```math
\{x \in \mathbb{R}^n : \|x - c\|_p \le r\}.
```

For example:

- when ``p = 1``, the ball is based on Manhattan distance;
- when ``p = 2``, the ball is based on Euclidean distance;
- when ``p = \infty``, the ball is based on maximum coordinate distance.

## Why balls are useful for search

Balls are useful when search is based on distance.

For example, an epsilon-ball query asks for all points within distance ``\epsilon`` of a query point. A ball gives a natural geometric representation of that search region.

Balls can also be compared with bounding volumes. This allows search algorithms to quickly decide whether a region is too far away, fully contained, or possibly intersecting the query region.

## Basic usage

```julia
using SearchableGeometries

ball = Ball([0.0, 0.0], 1.0; p=2)

is_contained(ball, [0.5, 0.5])
is_contained(ball, [2.0, 0.0])

bv = BoundingVolume(ball)
```

A ball can also be intersected with a bounding volume:

```julia
using SearchableGeometries

ball = Ball([0.0, 0.0], 1.0; p=2)
bv = BoundingVolume([0.5, -0.5], [2.0, 0.5])

intersects(bv, ball)
get_intersection(bv, ball)
```

## See also

- [`Ball`](@ref def_Ball)
- [`BoundingVolume`](@ref def_BoundingVolume)
- [`is_contained`](@ref)
- [`intersects`](@ref)
- [`get_intersection`](@ref)
- [`get_reduced_dim_ball`](@ref)
- [`tighten_bv_bounds!`](@ref)
