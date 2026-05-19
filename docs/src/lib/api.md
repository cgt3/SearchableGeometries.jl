# API Reference

```@meta
CurrentModule = SearchableGeometries
```

This page lists the public API documented throughout the manual.

## Abstract type

```@docs
SearchableGeometry
```

## Geometry types

- [`BoundingVolume`](@ref)
- [`Ball`](@ref)
- [`Hyperplane`](@ref)
- [`Cones`]

## Generic operations

- [`is_contained`](@ref)
- [`intersects`](@ref)
- [`get_intersection`](@ref)
- [`get_closest_point`](@ref)
- [`get_furthest_point`](@ref)

## Helpers

- [`face_index_to_spatial_index`](@ref)
- [`get_face_bounding_volume`](@ref)
- [`get_reduced_dim_ball`](@ref)
- [`tighten_bv_bounds!`](@ref)
- [`get_furthest_pt`](@ref)
- [`get_antifurthest_pt`](@ref)
- [`tighten_bv_bounds`](@ref)