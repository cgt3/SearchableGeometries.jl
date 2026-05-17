# SearchableGeometries.jl

```@meta
CurrentModule = SearchableGeometries
```

```@docs
SearchableGeometries
```

SearchableGeometries.jl is designed around the idea of making geometric objects searchable.

The long-term goal is to support spatial and geometric search structures, such as graph-based search and tree-based search methods. For these search routines to work, the package first needs reliable geometric primitives: objects that can answer containment, intersection, closest-point, furthest-point, and bounding queries.

At the moment, the package focuses on the geometry layer.

## Geometry Types

SearchableGeometries.jl currently provides:

- [`BoundingVolume`](@ref)
- [`Ball`](@ref)
- [`Hyperplane`](@ref)

## Common Geometric Operations

The main geometric operations are:

- [`is_contained`](@ref): test whether a point or object lies inside another object.
- [`intersects`](@ref): test whether two objects intersect.
- [`get_intersection`](@ref): compute a bounding representation of an intersection.
- [`get_closest_point`](@ref): find a closest point in a bounding volume or on a hyperplane.
- [`get_furthest_point`](@ref): find a furthest point in a bounding volume.