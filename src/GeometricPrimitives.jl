module GeometricPrimitives

using LinearAlgebra
import Base: ==, getindex

# Data types
export SearchableGeometry, BoundingVolume, Ball, Hyperplane, Line, Cone

# General functions
export is_contained, intersects, get_intersection, get_closest_point, get_furthest_point

const DEFAULT_BV_POINT_TOL = 1e-15

"""
    SearchableGeometry

Abstract supertype for geometric objects that can be queried by the package.

Concrete subtypes currently include:

- [`BoundingVolume`](@ref)
- [`Ball`](@ref)
- [`Hyperplane`](@ref)
"""
abstract type SearchableGeometry end

# ----------------------------------------------------------
# Generic geometric operations for documentation
# ----------------------------------------------------------

"""
    is_contained(container, query; kwargs...)

Test whether `query` is contained in `container`.

This is a generic containment predicate. Different methods are defined for
different geometric objects, such as bounding volumes, balls, and hyperplanes.

The meaning of containment depends on the types of the arguments. For example:

- `is_contained(bv, pt)` tests whether a point lies inside a bounding volume;
- `is_contained(ball, pt)` tests whether a point lies inside a ball;
- `is_contained(plane, pt)` tests whether a point lies on a hyperplane.

Many methods accept keyword arguments such as `include_boundary` and `tol`.

# See also

[`BoundingVolume`](@ref), [`Ball`](@ref), [`Hyperplane`](@ref)
"""
function is_contained end


"""
    intersects(a, b; kwargs...)

Test whether two geometric objects have a nonempty intersection.

This is a generic intersection predicate. Different methods are defined for
different pairs of objects, such as two bounding volumes, a bounding volume and
a ball, or a bounding volume and a hyperplane.

Many methods accept `include_boundary=true`, which controls whether boundary-only
contact counts as an intersection.

# See also

[`BoundingVolume`](@ref), [`Ball`](@ref), [`Hyperplane`](@ref)
"""
function intersects end


"""
    get_intersection(a, b; kwargs...)

Compute a representation of the intersection of two geometric objects.

For the current methods in this package, the return value is usually a
[`BoundingVolume`](@ref) that encloses the true intersection. If the objects do
not intersect, methods typically return `BoundingVolume()`, the empty bounding
volume.

# See also

[`intersects`](@ref), [`is_contained`](@ref), [`BoundingVolume`](@ref)
"""
function get_intersection end


"""
    get_closest_point(object, query; kwargs...)

Return a point associated with `object` that is closest to `query`.

The meaning of `query` depends on the method. For example:

- `get_closest_point(bv, pt)` returns the closest point in a bounding volume to a point;
- `get_closest_point(pt, plane)` projects a point onto a hyperplane;
- `get_closest_point(bv, plane)` returns a point in a bounding volume closest to a hyperplane.

# See also

[`get_furthest_point`](@ref), [`BoundingVolume`](@ref), [`Hyperplane`](@ref)
"""
function get_closest_point end


"""
    get_furthest_point(object, query; kwargs...)

Return a point associated with `object` that is furthest from `query`.

The meaning of `query` depends on the method. For example:

- `get_furthest_point(bv, pt)` returns a point in a bounding volume furthest from a point;
- `get_furthest_point(bv, plane)` returns a point in a bounding volume furthest from a hyperplane.

# See also

[`get_closest_point`](@ref), [`get_extreme_point`](@ref), [`get_antiextreme_point`](@ref)
"""
function get_furthest_point end

# ----------------------------------------------------------
# Lines
# ----------------------------------------------------------

struct Line
    source::Vector
    dir::Vector

    function Line(source::Vector, dir::Vector)
        if length(dir) != length(source)
            throw("SearchableGeometries.GeometricPrimitives.Line: direction vector and source point dimension must match")
        end

        return new(source, dir)# ./ norm(dir))
    end
end

function (line::Line)(s::Real)
    return line.source + line.dir * s
end

# include files for geometric primitives
include("boundingvolumes.jl")

include("balls.jl")

include("hyperplanes.jl")

include("cones.jl")
end # module GeometricPrimitives