# Bounding Volumes

```@meta
CurrentModule = SearchableGeometries
```

A bounding volume is an axis-aligned geometric region. In two dimensions, it is a rectangle. In three dimensions, it is a rectangular box. In higher dimensions, it is a hyperrectangle.

Bounding volumes are useful in search because they provide inexpensive geometric checks. Before a search algorithm performs a more expensive computation, it can first ask whether a bounding volume contains a point, intersects another bounding volume, or can be tightened around another geometric object.

## Mathematical definition

A bounding volume with lower bounds `lb` and upper bounds `ub` represents

```math
\{x \in \mathbb{R}^n : lb_i \le x_i \le ub_i,\quad \forall i \in \{1,\dots,n\}\}.
```

## Type

```@docs
BoundingVolume
```

## Point queries

```@docs
get_closest_point(::BoundingVolume, ::Vector{<:Real})
get_furthest_point(::BoundingVolume, ::Vector{<:Real})
```

## Containment

```@docs
is_contained(::BoundingVolume, ::Vector{<:Real})
is_contained(::BoundingVolume, ::BoundingVolume)
```

## Intersection with another bounding volume

```@docs
intersects(::BoundingVolume, ::BoundingVolume)
get_intersection(::BoundingVolume, ::BoundingVolume)
```

## Faces of a bounding volume

```@docs
face_index_to_spatial_index
get_face_bounding_volume
```