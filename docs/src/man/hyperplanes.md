# [Hyperplanes](@id hyperplanes_manual)

```@meta
CurrentModule = SearchableGeometries
```

A hyperplane is a flat affine surface used to split or query space.

In two dimensions, a hyperplane is a line. In three dimensions, it is an ordinary plane. In higher dimensions, it is a codimension-one affine subspace.

A hyperplane is represented by a point on the hyperplane and a normal vector.

## Mathematical definition

A hyperplane with point `p₀` and normal vector `n` represents

```math
\{x \in \mathbb{R}^n : n^T(x - p_0) = 0\}.
```

The vector `n` is perpendicular to the hyperplane. The expression

```math
n^T(x - p_0)
```

measures the signed offset of a point `x` from the hyperplane. If the value is zero, then the point lies on the hyperplane.

## Why hyperplanes are useful for search

Hyperplanes are useful because they can split space into two sides.

This makes them important in geometric search, spatial partitioning, and tree-based methods. For example, a kd-tree repeatedly splits space using coordinate-aligned hyperplanes. More general search structures may use arbitrary hyperplanes to separate regions.

Hyperplanes are also useful for testing whether a bounding volume crosses a decision boundary, lies completely on one side, or touches the boundary.

## Basic usage

```julia
using SearchableGeometries

# The line x - y = 0 in R²
plane = Hyperplane([0.0, 0.0], [1.0, -1.0])

is_contained(plane, [2.0, 2.0])
is_contained(plane, [2.0, 0.0])

get_closest_point([2.0, 0.0], plane)
```

A hyperplane can also be queried against a bounding volume:

```julia
using SearchableGeometries

bv = BoundingVolume([0.0, 0.0], [2.0, 2.0])

# The line x + y = 1
plane = Hyperplane([1.0, 0.0], [1.0, 1.0])

intersects(bv, plane)
get_intersection(bv, plane)
get_closest_point(bv, plane)
get_furthest_point(bv, plane)
```

## See also

- [`Hyperplane`](@ref def_Hyperplane)
- [`BoundingVolume`](@ref def_BoundingVolume)
- [`is_contained`](@ref)
- [`intersects`](@ref)
- [`get_intersection`](@ref)
- [`get_closest_point`](@ref)
- [`get_furthest_point`](@ref)