```@meta
CurrentModule = SearchableGeometries
```

# SearchableGeometries

[SearchableGeometries](https://github.com/cgt3/SearchableGeometries.jl) provides geometric primitives and query operations for working with axis-aligned bounding volumes and balls.

## Installation

```julia
using Pkg
Pkg.add("SearchableGeometries")
```

## Usage

```julia
using SearchableGeometries

# Create a bounding volume
bv = BoundingVolume([0.0, 0.0], [1.0, 1.0])

# Create a ball
ball = Ball([0.5, 0.5], 0.25)

# Check if a point is contained in the bounding volume
point = [0.5, 0.5]
println(isContained(bv, point))

# Check if a point is contained in the ball
println(isContained(ball, point))

# Check if two bounding volumes intersect
bv2 = BoundingVolume([0.5, 0.5], [1.5, 1.5])
println(intersects(bv, bv2))

# Get the intersection of two bounding volumes
intersection = getIntersection(bv, bv2)
println(intersection)
```

## API

```@autodocs
Modules = [SearchableGeometries]
Order = [:type, :function]
```