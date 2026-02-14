# SearchableGeometries.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://cgt3.github.io/SearchableGeometries.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://cgt3.github.io/SearchableGeometries.jl/dev/)
[![Build Status](https://github.com/cgt3/SearchableGeometries.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/cgt3/SearchableGeometries.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/cgt3/SearchableGeometries.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/cgt3/SearchableGeometries.jl)

A Julia package for defining and querying geometric primitives, specifically designed for efficient spatial indexing and searching operations.

## Features

- **Geometric Primitives**: Use `BoundingVolume` (Axis-Aligned Bounding Boxes) and `Ball` (defined by center, radius, and p-norm).
- **Dimension Agnostic**: Supports geometries in arbitrary dimensions.
- **Active Dimensions**: Ability to handle "active" and "inactive" dimensions, useful for subspace clustering or feature selection.
- **Efficient Queries**:
    - `isContained(geometry, query)`: Check containment.
    - `intersects(bv1, bv2)`: Check for intersection.
    - `getIntersection(bv1, bv2)`: Compute the intersection of two BoundingVolumes.
    - `getClosestPoint(bv, pt)` / `getFurthestPoint(bv, pt)`: Find points on the geometry relative to a query point.

## Installation

You can install `SearchableGeometries.jl` using the Julia package manager:

```julia
using Pkg
Pkg.add("SearchableGeometries")
```

## Usage

### 1. Working with Bounding Volumes

A `BoundingVolume` is defined by its lower bounds (`lb`) and upper bounds (`ub`).

```julia
using SearchableGeometries

# Define a 2D Bounding Volume from [0,0] to [1,1]
lb = [0.0, 0.0]
ub = [1.0, 1.0]
bv = BoundingVolume(lb, ub)

# Check points
pt_inside = [0.5, 0.5]
pt_outside = [1.5, 1.5]

println(isContained(bv, pt_inside))  # true
println(isContained(bv, pt_outside)) # false

# Get closest point on boundary
closest = getClosestPoint(bv, pt_outside)
println(closest) # [1.0, 1.0]
```

### 2. Working with Balls

A `Ball` is defined by a center, radius, and a p-norm (default Euclidean, `p=2`).

```julia
# Define a Ball at [0.5,0.5] with radius 0.25
center = [0.5, 0.5]
radius = 0.25
ball = Ball(center, radius)

# Check containment
println(isContained(ball, [6.0, 6.0])) # true (distance is sqrt(2) ≈ 1.414 < 2)
```

### 3. Interactions between Geometries

You can check if a `Ball` is contained within a `BoundingVolume`.

```julia
bv = BoundingVolume([0.0, 0.0], [1.0, 1.0])
ball_inside = Ball([0.5, 0.5], 0.25)
ball_outside = Ball([1.5, 1.5], 1.0)

println(isContained(bv, ball_inside))  # true
println(isContained(bv, ball_outside)) # false
```

### 4. Advanced: Active Dimensions

You can define geometries that only consider specific dimensions.

```julia
# A Ball in 3D space, but only active in dimensions 1 and 3 (index 2 is inactive)
# Effectively a cylinder along the y-axis
ball_subspace = Ball([0, 0, 0], 1; active_indices=true, indices=[1, 3])

# Point [0.5, 100, 0.5] is "contained" because axis 2 is ignored
pt = [0.5, 100, 0.5] 
println(isContained(ball_subspace, pt)) # true
```
