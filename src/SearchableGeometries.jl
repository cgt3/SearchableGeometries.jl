"""
    SearchableGeometries

Tools for representing and querying simple geometric objects used in search,
intersection, containment, and bounding-volume computations.

The package currently provides three main geometric objects:

- [`BoundingVolume`](@ref): an axis-aligned hyperrectangle described by lower
  and upper coordinate bounds.
- [`Ball`](@ref): a possibly lower-dimensional ``p``-norm ball embedded in a
  coordinate space.
- [`Hyperplane`](@ref): an affine hyperplane described by a point and a normal
  vector.

The public API is organized around a small set of generic operations:

- [`is_contained`](@ref): test whether a point or object lies inside another object.
- [`intersects`](@ref): test whether two objects intersect.
- [`get_intersection`](@ref): compute a bounding representation of an intersection.
- [`get_closest_point`](@ref): find a closest point in a bounding volume or on a hyperplane.
- [`get_furthest_point`](@ref): find a furthest point in a bounding volume.

# Boundary convention

Many predicates accept `include_boundary=true`. When `true`, points on the
boundary are treated as contained or intersecting. When `false`, strict interior
containment/intersection is required.

# Tolerance convention

Several methods accept `tol=DEFAULT_BV_POINT_TOL`. The tolerance is used to
handle floating-point roundoff, detect nearly inactive dimensions, and decide
whether values are close enough to geometric boundaries.
"""
module SearchableGeometries

using LinearAlgebra

# Submodels
export GeometricPrimitives
include("GeometricPrimitives.jl")
using .GeometricPrimitives

end # module SearchableGeometries