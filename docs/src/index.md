# SearchableGeometries.jl

```@meta
CurrentModule = SearchableGeometries
```

SearchableGeometries.jl provides geometric tools for building searchable geometric structures.

The long-term goal is to support spatial and graph-based search structures such as kd-trees, and bounding-volume hierarchies. The current package focuses on the geometric primitives needed to make those search routines possible.

## Introduction

Search algorithms often need to answer geometric questions quickly:

- Is this point inside this region?
- Do these two regions intersect?
- What is the closest point in a region to a query point?
- What is the furthest point in a region from a query object?
- Can a bounding volume be tightened around an intersection?

SearchableGeometries.jl provides types and methods for answering these questions.

## Example

```julia
using SearchableGeometries

bv = BoundingVolume([0.0, 0.0], [2.0, 1.0])
pt = [3.0, 0.5]

closest_pt = get_closest_point(bv, pt) # returns [2.0, 0.5]
```

## Features

The current package provides:

- axis-aligned bounding volumes;
- ``p``-norm balls;
- affine hyperplanes;
- containment tests;
- intersection tests;
- closest-point and furthest-point queries;
- bounding-volume tightening routines.

## Manual Outline

- [Getting Started](@ref getting_started)
- [Bounding Volumes](@ref bounding_volumes_manual)
- [Balls](@ref balls_manual)
- [Hyperplanes](@ref hyperplanes_manual)
- [Cones](@ref cones_manual)
- [Geometric Search](@ref geometric_search)

## Library Outline

- [API Reference](@ref)
- [BoundingVolume](@ref def_BoundingVolume)
- [Ball](@ref def_Ball)
- [Hyperplane](@ref def_Hyperplane)
- [Containment](@ref containment)
- [Intersection](@ref intersection)
- [Closest and Furthest Points](@ref closest_furthest_points)

