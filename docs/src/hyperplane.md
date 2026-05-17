# Hyperplanes

```@meta
CurrentModule = SearchableGeometries
```

A hyperplane is a flat affine surface. In two dimensions, a hyperplane is a line. In three dimensions, it is an ordinary plane. In higher dimensions, it is a codimension-one affine subspace.

A hyperplane is represented by a point `p0` on the hyperplane and a normal vector `n`:

```math
\{x \in \mathbb{R}^n : n^T(x - p_0) = 0\}.
```

Hyperplanes are useful in search because they can split space, define decision boundaries, and intersect bounding volumes.

## Type

```@docs
Hyperplane
```

## Point containment and projection

```@docs
is_contained(::Hyperplane, ::Vector{<:Real})
get_closest_point(::Vector{<:Real}, ::Hyperplane)
```

## Linear extrema over bounding volumes

```@docs
get_furthest_pt
get_antifurthest_pt
```

## Bounding-volume and hyperplane operations

```@docs
intersects(::BoundingVolume, ::Hyperplane)
tighten_bv_bounds
get_intersection(::BoundingVolume, ::Hyperplane)
get_closest_point(::BoundingVolume, ::Hyperplane)
get_furthest_point(::BoundingVolume, ::Hyperplane)
```