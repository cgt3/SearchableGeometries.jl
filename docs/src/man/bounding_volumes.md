# [Bounding Volumes](@id bounding_volumes_manual)

```@meta
CurrentModule = SearchableGeometries
```

A bounding volume is an axis-aligned region used to enclose a set of points or another geometric object.

In two dimensions, a bounding volume is a rectangle. In three dimensions, it is a rectangular box. In higher dimensions, it is a hyperrectangle.

## Mathematical definition

A bounding volume with lower bounds `lb` and upper bounds `ub` represents

```math
\{x \in \mathbb{R}^n : lb_i \le x_i \le ub_i,\quad i=1,\dots,n\}.
```

## Why bounding volumes are useful for search

Bounding volumes are cheap to query. A search algorithm can use them to quickly rule out regions that cannot contain the answer.

For example, if a query point is far from a bounding volume, then every point inside that bounding volume is also far from the query in a controlled way.

## Basic usage

```julia
using SearchableGeometries

bv = BoundingVolume([0.0, 0.0], [2.0, 1.0])

is_contained(bv, [1.0, 0.5])
get_closest_point(bv, [3.0, 0.5])
get_furthest_point(bv, [3.0, 0.5])
```

## See also
- [BoundingVolume](@ref def_BoundingVolume)
- [is_contained](@ref)
- [intersects](@ref)
- [get_intersection](@ref)
