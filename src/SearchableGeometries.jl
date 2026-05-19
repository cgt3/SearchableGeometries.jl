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

# Data Types:
export SearchableGeometry, Ball, BoundingVolume, Hyperplane

# BV only Functions:
export get_closest_point, get_furthest_point, face_index_to_spatial_index, get_face_bounding_volume

# General Functions:
export is_contained, intersects, get_intersection, tighten_bv_bounds!

# Ball only Functions:
export get_reduced_dim_ball

# Hyperplane only Functions:
export get_furthest_pt, get_antifurthest_pt, tighten_bv_bounds

import Base.getindex

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

# Generic geometric operations ----------------------------------------------------------

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

[`intersects`](@ref), [`get_intersection`](@ref), [`BoundingVolume`](@ref),
[`Ball`](@ref), [`Hyperplane`](@ref)
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

[`is_contained`](@ref), [`get_intersection`](@ref), [`BoundingVolume`](@ref),
[`Ball`](@ref), [`Hyperplane`](@ref)
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

[`get_closest_point`](@ref), [`get_furthest_pt`](@ref), [`get_antifurthest_pt`](@ref)
"""
function get_furthest_point end

# Bounding Volumes ----------------------------------------------------------------------
@doc raw"""
    BoundingVolume()
    BoundingVolume(lb, ub; tol=DEFAULT_BV_POINT_TOL)

Construct an axis-aligned bounding volume.

A `BoundingVolume` represents the hyperrectangle

```math
\{x \in \mathbb{R}^n : lb_i \le x_i \le ub_i,\quad \forall i \in \{1, \dots, n\}\}.
```

The lower and upper bounds are stored in the fields `lb` and `ub`. Dimensions
where `lb[i]` and `ub[i]` are equal up to `tol` are treated as inactive. This
allows the same type to represent full-dimensional boxes, lower-dimensional
faces, line segments, and points.

Calling `BoundingVolume()` constructs an empty bounding volume. Empty bounding
volumes are used as the return value for intersections that do not exist.

# Arguments

- `lb::Vector{<:Real}`: lower coordinate bounds.
- `ub::Vector{<:Real}`: upper coordinate bounds.
- `tol::Real`: tolerance used to identify inactive dimensions.

# Fields

- `lb`: lower bounds.
- `ub`: upper bounds.
- `is_empty`: whether the object represents the empty set.
- `dim`: active geometric dimension.
- `active_dim`: indices whose lower and upper bounds differ by at least `tol`.
- `inactive_dim`: indices whose lower and upper bounds are equal up to `tol`.
- `is_active`: Boolean mask indicating active coordinates.

# Throws

Throws an error if `lb` and `ub` have different lengths, or if any coordinate
satisfies `lb[i] > ub[i]`.

# Examples

```julia
using SearchableGeometries

bv = BoundingVolume([0.0, 0.0], [2.0, 1.0])
bv.dim          # 2
bv.active_dim   # [1, 2]
```

A lower-dimensional bounding volume can represent a line segment:

```julia
segment = BoundingVolume([0.0, 0.0], [0.0, 1.0])
segment.dim          # 1
segment.inactive_dim # [1]
```

# See also

[`Ball`](@ref), [`Hyperplane`](@ref)
"""
struct BoundingVolume <: SearchableGeometry
    lb::Vector                  # lower bounds
    ub::Vector                  # upper bounds
    is_empty::Bool              # is the bounding volume empty?
    dim::Integer                # dimension of the bounding volume
    active_dim::Vector          # active dimensions
    inactive_dim::Vector        # inactive dimensions
    is_active::Vector{Bool}     # is the dimension active?

    # For invalid/empty Bounding Volumes (note an empty BV differs from a non-empty BV with dimension 0 (a point))
    function BoundingVolume()
        return new([Inf], [-Inf], true, 0, [], [], Bool[])
    end

    function BoundingVolume(
        lb::Vector{<:Real}, ub::Vector{<:Real}; tol::Real=DEFAULT_BV_POINT_TOL
    )
        if length(lb) != length(ub)
            throw("SearchableGeometries: BoundingVolume: lb (length: $(length(lb))) and ub (length: $(length(ub))) have different dimensions.")
        elseif any(lb .> ub)
            throw("SearchableGeometries: BoundingVolume: Cannot construct bounding volume with lb (=$lb) > ub (=$ub).")
        end

        dim = length(lb)
        is_active = ones(Bool, length(lb))
        for d in eachindex(lb)
            if abs(lb[d] - ub[d]) < tol
                dim -= 1
                is_active[d] = false
            end
        end
        all_dim = [eachindex(lb)...]
        active_dim = all_dim[is_active]
        inactive_dim = all_dim[is_active.!=true]

        return new(lb, ub, false, dim, active_dim, inactive_dim, is_active)
    end
end

import Base.==
function Base.:(==)(bv1::BoundingVolume, bv2::BoundingVolume)
    return all(bv1.lb .== bv2.lb) &&
           all(bv1.ub .== bv2.ub) &&
           bv1.is_empty == bv2.is_empty &&
           bv1.dim == bv2.dim &&
           all(bv1.active_dim .== bv2.active_dim) &&
           all(bv1.inactive_dim .== bv2.inactive_dim) &&
           all(bv1.is_active .== bv2.is_active)
end

@doc raw"""
    get_closest_point(bv::BoundingVolume, query_pt)

Return the point in `bv` closest to `query_pt` in Euclidean coordinate distance.

This operation is the coordinate-wise projection of `query_pt` onto the interval
`[bv.lb[i], bv.ub[i]]` in each dimension. If the point is already inside the
bounding volume, the returned point equals `query_pt`.

# Arguments

- `bv::BoundingVolume`: bounding volume to project onto.
- `query_pt::Vector{<:Real}`: point in the same embedding dimension as `bv`.

# Returns

A vector representing the closest point in `bv`.

# Examples

```julia
using SearchableGeometries

bv = BoundingVolume([0.0, 0.0], [1.0, 2.0])
get_closest_point(bv, [3.0, 0.5]) # returns [1.0, 0.5]
get_closest_point(bv, [0.25, 1.0]) # returns [0.25, 1.0]
```

# See also 

[`get_furthest_point(::BoundingVolume, ::Vector{<:Real})`](@ref)
"""
function get_closest_point(bv::BoundingVolume, query_pt::Vector{<:Real})
    closest_pt = copy(query_pt)

    I_lb = query_pt .< bv.lb
    I_ub = query_pt .> bv.ub

    closest_pt[I_lb] .= bv.lb[I_lb]
    closest_pt[I_ub] .= bv.ub[I_ub]

    return closest_pt
end

"""
    get_furthest_point(bv::BoundingVolume, query_pt)

Return a point in `bv` farthest from `query_pt`.

For an axis-aligned bounding volume, a farthest point is obtained by choosing,
in each coordinate, the bound farthest from the corresponding coordinate of
`query_pt`. If both bounds are equally far in a coordinate, the implementation
uses a deterministic tie-breaking rule.

# Arguments

- `bv::BoundingVolume`: bounding volume to search.
- `query_pt::Vector{<:Real}`: query point.

# Returns

A vector representing one farthest point in `bv`.

# Examples

```julia
using SearchableGeometries

bv = BoundingVolume([0.0, 0.0], [2.0, 3.0])
get_furthest_point(bv, [0.25, 2.8]) # approximately [2.0, 0.0]
```

# See also

[`get_closest_point(::BoundingVolume, ::Vector{<:Real})`](@ref)
"""
function get_furthest_point(bv::BoundingVolume, query_pt::Vector{<:Real})
    furthest_pt = similar(query_pt)

    ub_is_closer = 0.5 * (bv.ub + bv.lb) .<= query_pt
    lb_is_closer = ub_is_closer .== false

    furthest_pt[ub_is_closer] .= bv.lb[ub_is_closer]
    furthest_pt[lb_is_closer] .= bv.ub[lb_is_closer]

    return furthest_pt
end

@doc raw"""
    is_contained(bv::BoundingVolume, query_pt; include_boundary=true)

Test whether `query_pt` lies inside the bounding volume `bv`.

When `include_boundary=true`, the test is

```math
bv.lb_i \le query\_pt_i \le bv.ub_i \quad \forall i \in \{1, \dots, n\}
```

When `include_boundary=false`, strict inequalities are used.

# Arguments

- `bv::BoundingVolume`: containing bounding volume.
- `query_pt::Vector{<:Real}`: point to test.
- `include_boundary::Bool`: whether boundary points count as contained.

# Returns

`true` if the point is contained in `bv`; otherwise `false`.

# Examples

```julia
using SearchableGeometries

bv = BoundingVolume([0.0, 0.0], [1.0, 1.0])
is_contained(bv, [1.0, 0.5])                         # true
is_contained(bv, [1.0, 0.5]; include_boundary=false) # false
```

# See also

[`is_contained(::BoundingVolume, ::BoundingVolume)`](@ref)
"""
function is_contained(bv::BoundingVolume, query_pt::Vector{<:Real}; include_boundary::Bool=true)
    if (include_boundary && all(bv.lb .<= query_pt .<= bv.ub)) ||
       (!include_boundary && all(bv.lb .< query_pt .< bv.ub))
        return true
    else
        return false
    end
end

"""
    is_contained(bv::BoundingVolume, query_bv::BoundingVolume; include_boundary=true)

Test whether one bounding volume is contained in another.

This returns `true` when every point of `query_bv` lies inside `bv`. With
`include_boundary=true`, touching boundaries count as contained. With
`include_boundary=false`, `query_bv` must lie strictly inside `bv`.

# Examples

```julia
using SearchableGeometries

outer = BoundingVolume([0.0, 0.0], [2.0, 2.0])
inner = BoundingVolume([0.5, 0.5], [1.0, 1.0])
is_contained(outer, inner) # true
```

# See also

[`is_contained(::BoundingVolume, ::Vector{<:Real})`](@ref)
"""
function is_contained(bv::BoundingVolume, query_bv::BoundingVolume; include_boundary::Bool=true)
    if (!include_boundary && (all(query_bv.ub .< bv.ub) && all(query_bv.lb .> bv.lb))) ||
       (include_boundary && (all(query_bv.ub .<= bv.ub) && all(query_bv.lb .>= bv.lb)))
        return true
    else
        return false
    end
end

"""
    intersects(bv1::BoundingVolume, bv2::BoundingVolume; include_boundary=true)

Test whether two axis-aligned bounding volumes intersect.

Two bounding volumes intersect when their coordinate intervals overlap in every
dimension. If `include_boundary=true`, touching at a boundary, edge, or corner
counts as an intersection. If `include_boundary=false`, the interiors must
overlap in every active coordinate.

# Arguments

- `bv1::BoundingVolume`: first bounding volume.
- `bv2::BoundingVolume`: second bounding volume.
- `include_boundary::Bool`: whether boundary-only contact counts as intersection.

# Returns

`true` if the bounding volumes intersect; otherwise `false`.

# Examples

```julia
using SearchableGeometries

bv1 = BoundingVolume([0.0, 0.0], [1.0, 1.0])
bv2 = BoundingVolume([1.0, 0.25], [2.0, 0.75])

intersects(bv1, bv2)                         # true
intersects(bv1, bv2; include_boundary=false) # false
```

# See also

[`get_intersection(::BoundingVolume, ::BoundingVolume)`](@ref)
"""
function intersects(bv1::BoundingVolume, bv2::BoundingVolume; include_boundary::Bool=true)
    if (include_boundary && (any(bv1.lb .> bv2.ub) || any(bv1.ub .< bv2.lb))) ||
       (!include_boundary && (any(bv1.lb .>= bv2.ub) || any(bv1.ub .<= bv2.lb)))
        return false
    else
        return true
    end
end

@doc raw"""
    get_intersection(bv1::BoundingVolume, bv2::BoundingVolume; tol=DEFAULT_BV_POINT_TOL)

Return the bounding volume representing the intersection of two bounding volumes.

The intersection is computed coordinate-wise using

```math
new\_lb_i = \max(bv1.lb_i, bv2.lb_i), \qquad
new\_ub_i = \min(bv1.ub_i, bv2.ub_i).
```

If the two bounding volumes do not intersect, this method returns
`BoundingVolume()`, the empty bounding volume.

# Arguments

- `bv1::BoundingVolume`: first bounding volume.
- `bv2::BoundingVolume`: second bounding volume.
- `tol::Real`: tolerance passed to the returned [`BoundingVolume`](@ref).

# Returns

A [`BoundingVolume`](@ref). The returned object may be empty.

# Examples

```julia
using SearchableGeometries

bv1 = BoundingVolume([0.0, 0.0], [2.0, 2.0])
bv2 = BoundingVolume([1.0, -1.0], [3.0, 1.0])

get_intersection(bv1, bv2)
# BoundingVolume with lb = [1.0, 0.0], ub = [2.0, 1.0]
```

# See also

[`intersects(::BoundingVolume, ::BoundingVolume)`](@ref)
"""
function get_intersection(bv1::BoundingVolume, bv2::BoundingVolume; tol::Real=DEFAULT_BV_POINT_TOL)
    if bv1.is_empty || bv2.is_empty
        return BoundingVolume()
    end

    new_lb = max.(bv1.lb, bv2.lb)
    new_ub = min.(bv1.ub, bv2.ub)
    if any(new_lb .> new_ub)
        return BoundingVolume()
    end

    return BoundingVolume(new_lb, new_ub; tol=tol)
end

"""
    face_index_to_spatial_index(face_index, num_dim)

Map a bounding-volume face index to its corresponding coordinate index.

For an `n`-dimensional bounding volume, faces are indexed from `1` to `2n`.
The first `n` indices correspond to lower-bound faces, and the next `n` indices
correspond to upper-bound faces.

# Arguments

- `face_index::Integer`: face index in `1:2*num_dim`.
- `num_dim::Integer`: embedding dimension of the bounding volume.

# Returns

The coordinate index associated with the face.

# Examples

```julia
using SearchableGeometries

face_index_to_spatial_index(2, 3) # 2, lower face in coordinate 2
face_index_to_spatial_index(5, 3) # 2, upper face in coordinate 2
```

# See also

[`get_face_bounding_volume`](@ref)
"""
function face_index_to_spatial_index(face_index::Integer, num_dim::Integer)
    return face_index <= num_dim ? face_index : face_index - num_dim
end

"""
    get_face_bounding_volume(face_index, bv::BoundingVolume; tol=DEFAULT_BV_POINT_TOL)

Return the lower- or upper-bound face of a bounding volume as a new
[`BoundingVolume`](@ref).

For a bounding volume embedded in `n` dimensions, faces are indexed as follows:

- `1:n`: lower-bound faces, where coordinate `i` is fixed to `bv.lb[i]`.
- `n+1:2n`: upper-bound faces, where coordinate `i-n` is fixed to `bv.ub[i-n]`.

The returned face has one coordinate fixed and therefore may have lower
geometric dimension than the original bounding volume.

# Arguments

- `face_index::Integer`: face index.
- `bv::BoundingVolume`: bounding volume whose face is requested.
- `tol::Real`: tolerance passed to the returned [`BoundingVolume`](@ref).

# Returns

A [`BoundingVolume`](@ref) representing the requested face.

# Examples

```julia
using SearchableGeometries

bv = BoundingVolume([0.0, 0.0], [2.0, 3.0])
left_face = get_face_bounding_volume(1, bv)
top_face = get_face_bounding_volume(4, bv)
```

# See also

[`face_index_to_spatial_index`](@ref), [`tighten_bv_bounds!`](@ref)
"""
function get_face_bounding_volume(face_index::Integer, bv::BoundingVolume; tol::Real=DEFAULT_BV_POINT_TOL)
    face_lb, face_ub = copy(bv.lb), copy(bv.ub)

    if face_index <= length(bv.lb) # Lower bound face
        face_ub[face_index] = bv.lb[face_index]
    else # Upper bound face
        d = face_index - length(bv.lb)
        face_lb[d] = bv.ub[d]
    end

    return BoundingVolume(face_lb, face_ub; tol=tol)
end


# Balls ---------------------------------------------------------------------------------

@doc raw"""
    Ball(center, radius; p=2, active_indices=true, indices=eachindex(center))

Construct a possibly lower-dimensional ``p``-norm ball.

A full-dimensional ball represents the set

```math
\{x : \|x - c\|_p \le r\},
```

where `c` is `center`, `r` is `radius`, and `p` is the norm parameter.
The ball may also be restricted to a subset of active coordinates. Inactive
coordinates are fixed to the corresponding entries of `center`.

# Arguments

- `center::Vector{<:Real}`: center of the ball in the embedding space.
- `radius::Real`: nonnegative radius.
- `p::Real`: norm parameter. Common values are `1`, `2`, and `Inf`.
- `active_indices::Bool`: controls how `indices` is interpreted.
- `indices::Vector{Int}`: active or inactive coordinate indices.

# Active and inactive dimensions

If `active_indices=true`, then `indices` gives the active dimensions of the
ball. If `active_indices=false`, then `indices` gives the inactive dimensions.
Inactive coordinates are fixed to the ball center.

# Fields

- `center`: center of the ball.
- `radius`: radius of the ball.
- `p`: norm parameter.
- `dim`: active geometric dimension.
- `active_dim`: active coordinate indices.
- `inactive_dim`: inactive coordinate indices.
- `is_active`: Boolean mask indicating active coordinates.
- `embedding_dim`: dimension of the ambient coordinate space.

# Throws

Throws an error if `radius < 0`. A zero-radius ball is allowed and represents a
point.

# Examples

```julia
using SearchableGeometries

ball = Ball([0.0, 0.0], 1.0; p=2)

# A one-dimensional ball embedded in 2D: the line segment
# {(x, 0): |x| <= 1}
segment_ball = Ball([0.0, 0.0], 1.0; p=2, active_indices=true, indices=[1])
```

# See also

[`BoundingVolume`](@ref), [`get_reduced_dim_ball`](@ref), [`is_contained`](@ref)
"""
struct Ball <: SearchableGeometry
    center::Vector                      # center of the ball
    radius::Real                        # radius of the ball
    p::Real                             # p-norm
    dim::Integer                        # dimension of the ball
    active_dim::Vector                  # active dimensions
    inactive_dim::Vector                # inactive dimensions
    is_active::Vector{Bool}             # is the dimension active?
    embedding_dim::Integer              # embedding dimension

    function Ball(
        center::Vector{<:Real}, radius::Real; p::Real=2::Real,
        active_indices::Bool=true::Bool, indices=(active_indices ? [eachindex(center)...] : Vector{Int}[])::Vector{Int}
    )
        if radius < 0
            throw("SearchableGeometries.Ball: Cannot construct ball with negative radius.")
        elseif radius == 0
            return new(center, radius, p, 0, [], [eachindex(center)...], zeros(Bool, length(center)), length(center))
        end

        unique_indices = unique(indices)
        if active_indices
            is_active = zeros(Bool, length(center))
            is_active[unique_indices] .= true

            all_dim = [eachindex(center)...]
            inactive_dim = all_dim[is_active.==false]
            dim = sum(is_active)
            return new(center, radius, p, dim, unique_indices, inactive_dim, is_active, length(center))

        else
            is_active = ones(Bool, length(center))
            is_active[unique_indices] .= false

            all_dim = [eachindex(center)...]
            active_dim = all_dim[is_active]
            dim = sum(is_active)
            return new(center, radius, p, dim, active_dim, unique_indices, is_active, length(center))
        end
    end
end

import Base.==
function Base.:(==)(ball1::Ball, ball2::Ball)
    return all(ball1.center .== ball2.center) &&
           ball1.radius == ball2.radius &&
           ball1.p == ball2.p &&
           ball1.dim == ball2.dim &&
           all(ball1.active_dim .== ball2.active_dim) &&
           all(ball1.inactive_dim .== ball2.inactive_dim) &&
           all(ball1.is_active .== ball2.is_active) &&
           ball1.embedding_dim == ball2.embedding_dim
end

"""
    BoundingVolume(ball::Ball; tol=DEFAULT_BV_POINT_TOL)

Construct the axis-aligned bounding volume enclosing `ball`.

For a ball with center `c` and radius `r`, the returned bounding volume has
bounds `c .- r` and `c .+ r`. This is an enclosing box for the ball, not
necessarily an exact representation of the ball unless `p=Inf` and the ball is
full-dimensional.

# Arguments

- `ball::Ball`: ball to enclose.
- `tol::Real`: tolerance passed to the returned [`BoundingVolume`](@ref).

# Returns

A [`BoundingVolume`](@ref) enclosing `ball`.

# Examples

```julia
using SearchableGeometries

ball = Ball([1.0, 2.0], 0.5)
bv = BoundingVolume(ball) # lb = [0.5, 1.5], ub = [1.5, 2.5]
```

# See also

[`Ball`](@ref), [`get_intersection`](@ref)
"""
function BoundingVolume(ball::Ball; tol::Real=DEFAULT_BV_POINT_TOL)
    lb = ball.center .- ball.radius
    ub = ball.center .+ ball.radius
    return BoundingVolume(lb, ub; tol=tol)
end

@doc raw"""
    is_contained(ball::Ball, query_pt; include_boundary=true, tol=DEFAULT_BV_POINT_TOL)

Test whether `query_pt` lies inside `ball`.

For a full-dimensional ball, this checks whether

```math
\|query\_pt - ball.center\|_p \le ball.radius.
```

For a lower-dimensional ball, inactive coordinates must match the ball center
up to `tol`, and the norm is computed only on active coordinates.

# Arguments

- `ball::Ball`: containing ball.
- `query_pt::Vector{<:Real}`: point to test.
- `include_boundary::Bool`: whether points exactly on the boundary count as contained.
- `tol::Real`: tolerance for inactive coordinates and boundary checks.

# Returns

`true` if the point is contained in the ball; otherwise `false`.

# Throws

Throws an error if the point dimension does not match `ball.embedding_dim`.

# Examples

```julia
using SearchableGeometries

ball = Ball([0.0, 0.0], 1.0; p=2)
is_contained(ball, [0.5, 0.5]) # true
is_contained(ball, [1.0, 0.0]; include_boundary=false) # false
```

# See also

[`Ball`](@ref), [`BoundingVolume`](@ref), [`intersects`](@ref)
"""
function is_contained(ball::Ball, query_pt::Vector{<:Real}; include_boundary::Bool=true, tol::Real=DEFAULT_BV_POINT_TOL)
    if length(query_pt) != ball.embedding_dim
        throw("Point dimension($(length(query_pt))) does not match ball embedding dimension($(ball.embedding_dim))")
    end

    if ball.dim < length(ball.center) # The ball does not have full dimension
        for d_fixed in ball.inactive_dim
            if abs(query_pt[d_fixed] - ball.center[d_fixed]) > tol
                return false
            end
        end
        R_query = norm(query_pt[ball.active_dim] - ball.center[ball.active_dim], ball.p)
    else # The ball has full dimension
        R_query = norm(query_pt - ball.center, ball.p)
    end

    return include_boundary ? R_query <= ball.radius : R_query < ball.radius
end

"""
    is_contained(bv::BoundingVolume, query_ball::Ball; include_boundary=true, tol=DEFAULT_BV_POINT_TOL)

Test whether `query_ball` is contained in the bounding volume `bv`.

This method first converts the ball to its enclosing [`BoundingVolume`](@ref)
and then tests bounding-volume containment. Therefore, it is exact for deciding
whether the entire ball is inside `bv`, because a ball is contained in an
axis-aligned box if and only if its axis-aligned bounding box is contained in
that box.

# Examples

```julia
using SearchableGeometries

bv = BoundingVolume([-2.0, -2.0], [2.0, 2.0])
ball = Ball([0.0, 0.0], 1.0)
is_contained(bv, ball) # true
```

# See also

[`Ball`](@ref), [`BoundingVolume`](@ref), [`intersects`](@ref)
"""
function is_contained(bv::BoundingVolume, query_ball::Ball; include_boundary::Bool=true, tol::Real=DEFAULT_BV_POINT_TOL)
    return is_contained(bv, BoundingVolume(query_ball; tol=tol); include_boundary=include_boundary)
end

"""
    is_contained(ball::Ball, query_bv::BoundingVolume; include_boundary=true)

Test whether a bounding volume is contained in a ball.

For a convex ball and an axis-aligned bounding volume, the maximum distance from
the ball center to the bounding volume occurs at a corner. This method finds a
farthest point of `query_bv` from `ball.center` and checks whether that point is
contained in `ball`.

# Examples

```julia
using SearchableGeometries

ball = Ball([0.0, 0.0], 2.0)
bv = BoundingVolume([-0.5, -0.5], [0.5, 0.5])
is_contained(ball, bv) # true
```

# See also

[`get_furthest_point`](@ref), [`intersects`](@ref)
"""
function is_contained(ball::Ball, query_bv::BoundingVolume; include_boundary::Bool=true)
    furthest_pt = get_furthest_point(query_bv, ball.center)
    return is_contained(ball, furthest_pt; include_boundary=include_boundary)
end

"""
    intersects(bv::BoundingVolume, ball::Ball; include_boundary=true, tol=DEFAULT_BV_POINT_TOL)

Test whether a bounding volume and a ball intersect.

The method first performs inexpensive checks using the ball's enclosing bounding
volume. If those checks are inconclusive, it projects the ball center onto the
bounding volume and tests whether that closest point lies in the ball.

# Arguments

- `bv::BoundingVolume`: bounding volume.
- `ball::Ball`: ball.
- `include_boundary::Bool`: whether boundary-only contact counts as intersection.
- `tol::Real`: tolerance used in geometric comparisons.

# Returns

`true` if `bv` and `ball` intersect; otherwise `false`.

# Examples

```julia
using SearchableGeometries

bv = BoundingVolume([1.0, -0.5], [2.0, 0.5])
ball = Ball([0.0, 0.0], 1.0)
intersects(bv, ball) # true, touching at x = 1
```

# See also

[`get_intersection`](@ref), [`get_closest_point`](@ref)
"""
function intersects(bv::BoundingVolume, ball::Ball; include_boundary::Bool=true, tol::Real=DEFAULT_BV_POINT_TOL)
    # First, do the easy checks against the ball's BV:
    if !intersects(bv, BoundingVolume(ball; tol=tol); include_boundary=include_boundary)
        # The two are completely disjoint
        return false
    elseif is_contained(ball, bv; include_boundary=include_boundary)
        # The ball is completely contained in the BV
        return true
    end
    # Note: the below check could catch all of the above cases but its use of the
    #       lp-norm makes it more expensive. The above checks are to avoid having
    #       to compute norms.

    # Check if the closest point on the BV to the ball's center is inside the ball
    closest_pt = get_closest_point(bv, ball.center)
    return is_contained(ball, closest_pt; include_boundary=include_boundary)
end

@doc raw"""
    get_reduced_dim_ball(removal_dim, x_d, ball::Ball)

Slice `ball` by the coordinate plane `x[removal_dim] = x_d`.

The returned ball is embedded in the same coordinate space as `ball`, but the
coordinate `removal_dim` is made inactive and fixed to `x_d`. For finite `p`,
the new radius is

```math
\left(r^p - |x_d - c_d|^p\right)^{1/p}.
```

For `p == Inf`, the radius is unchanged whenever the slice intersects the ball.

# Arguments

- `removal_dim::Integer`: coordinate to fix.
- `x_d::Real`: fixed coordinate value.
- `ball::Ball`: ball to slice.

# Returns

A lower-dimensional [`Ball`](@ref) representing the coordinate slice.

# Throws

Throws an error if the coordinate plane does not intersect the ball.

# Examples

```julia
using SearchableGeometries

ball = Ball([0.0, 0.0], 1.0; p=2)
slice = get_reduced_dim_ball(2, 0.0, ball)
slice.dim          # 1
slice.inactive_dim # [2]
```

# See also

[`Ball`](@ref), [`tighten_bv_bounds!`](@ref)
"""
function get_reduced_dim_ball(removal_dim::Integer, x_d::Real, ball::Ball)
    if x_d < ball.center[removal_dim] - ball.radius || ball.center[removal_dim] + ball.radius < x_d
        throw("SearchableGeometries.Ball: coordinate plane defined by x_$removal_dim = $x_d does not intersect the ball (center=$(ball.center), radius=$(ball.radius))")
    end

    new_center = copy(ball.center)
    new_center[removal_dim] = x_d

    new_radius = ball.p === Inf ? ball.radius : (ball.radius^ball.p - abs(x_d - ball.center[removal_dim])^ball.p)^(1 / ball.p)
    inactive_dim = [ball.inactive_dim..., removal_dim]

    return Ball(new_center, new_radius; p=ball.p, active_indices=false, indices=inactive_dim)
end

"""
    tighten_bv_bounds!(bv::BoundingVolume, ball::Ball; tol=DEFAULT_BV_POINT_TOL)

Mutate `bv` by tightening its bounds around the intersection with `ball`.

This method shrinks the axis-aligned bounding volume so that it more tightly
encloses `bv ∩ ball`. It does not return the exact curved intersection; instead,
it returns the indices of bounds that were changed and mutates `bv` in place.

For one-dimensional balls, the update is direct. For higher-dimensional balls,
the method recursively examines faces of the bounding volume and tightens bounds
when faces do not intersect the ball.

# Arguments

- `bv::BoundingVolume`: bounding volume to mutate.
- `ball::Ball`: ball used to tighten `bv`.
- `tol::Real`: tolerance used in intersection checks.

# Returns

A pair `(altered_lb_indices, altered_ub_indices)` where each entry is a vector
of coordinate indices whose lower or upper bounds were modified.

# Mutating behavior

This function modifies `bv.lb` and/or `bv.ub` in place.

# Examples

```julia
using SearchableGeometries

bv = BoundingVolume([-3.0, -3.0], [3.0, 3.0])
ball = Ball([1.0, 1.0], 1.0)
changed_lb, changed_ub = tighten_bv_bounds!(bv, ball)
```

# See also

[`get_intersection`](@ref), [`get_reduced_dim_ball`](@ref)
"""
function tighten_bv_bounds!(bv::BoundingVolume, ball::Ball; tol::Real=DEFAULT_BV_POINT_TOL)
    if ball.dim == 1
        d = ball.active_dim[1]
        lb_ball = ball.center[d] - ball.radius
        ub_ball = ball.center[d] + ball.radius

        if bv.lb[d] < lb_ball < bv.ub[d]
            bv.lb[d] = lb_ball
            altered_lb_indices = [d]
        else
            altered_lb_indices = []
        end

        if bv.lb[d] < ub_ball < bv.ub[d]
            bv.ub[d] = ub_ball
            altered_ub_indices = [d]
        else
            altered_ub_indices = []
        end
        return altered_lb_indices, altered_ub_indices
    end

    # For non-simple intersections
    ub_pt_projected = copy(ball.center)
    lb_pt_projected = copy(ball.center)

    # For every face with no intersection with the ball, recurse
    altered_lb_indices = []
    altered_ub_indices = []
    num_dim = length(ball.center)
    for f_target in 1:2*num_dim # for each face
        face_bv = get_face_bounding_volume(f_target, bv, tol=tol)

        non_simple = false
        if !intersects(face_bv, ball, include_boundary=true, tol=tol)
            adjacent_faces = [1:f_target-1..., f_target+1:2*num_dim...]
            d_target = face_index_to_spatial_index(f_target, num_dim)

            # Check if this face needs to be updated using a non-simple intersection
            if f_target <= num_dim # f_target is a lb face
                lb_pt_projected[d_target] = face_bv.lb[d_target]
                if is_contained(face_bv, lb_pt_projected, include_boundary=true)
                    push!(altered_lb_indices, d_target)
                    bv.lb[d_target] = ball.center[d_target] - ball.radius
                    non_simple = true
                end
                lb_pt_projected[d_target] = ball.center[d_target]
            else # f_target is an ub face
                ub_pt_projected[d_target] = face_bv.ub[d_target]
                if is_contained(face_bv, ub_pt_projected, include_boundary=true)
                    push!(altered_ub_indices, d_target)
                    bv.ub[d_target] = ball.center[d_target] + ball.radius
                    non_simple = true
                end
                ub_pt_projected[d_target] = ball.center[d_target]
            end

            # For simple intersections
            if !non_simple
                for f_adjacent in adjacent_faces
                    adjacent_face_bv = get_face_bounding_volume(f_adjacent, bv, tol=tol)

                    if intersects(adjacent_face_bv, ball, include_boundary=true)
                        d_fixed = face_index_to_spatial_index(f_adjacent, num_dim)
                        reduced_ball = get_reduced_dim_ball(d_fixed, adjacent_face_bv.lb[d_fixed], ball)

                        # This will modify face_adjacent's bounds 
                        altered_lb_indices_new, altered_ub_indices_new = tighten_bv_bounds!(adjacent_face_bv, reduced_ball, tol=tol)

                        # Update the higher-dim BV with the new bounds on face_adjacent
                        bv.lb[altered_lb_indices_new] .= adjacent_face_bv.lb[altered_lb_indices_new]
                        bv.ub[altered_ub_indices_new] .= adjacent_face_bv.ub[altered_ub_indices_new]

                        altered_lb_indices = vcat(altered_lb_indices, altered_lb_indices_new)
                        altered_ub_indices = vcat(altered_ub_indices, altered_ub_indices_new)
                    end
                end # for
            end # if simple
        end
    end

    return altered_lb_indices, altered_ub_indices
end

"""
    get_intersection(bv::BoundingVolume, ball::Ball; tol=DEFAULT_BV_POINT_TOL)

Return a bounding volume enclosing the intersection `bv ∩ ball`.

If the bounding volume and ball do not intersect, this method returns the empty
bounding volume `BoundingVolume()`. Otherwise, it first intersects `bv` with the
ball's axis-aligned bounding volume and then tightens the result when possible.

The returned object is a [`BoundingVolume`](@ref), not a curved geometric
region. It is intended to be a tight axis-aligned enclosure of the true
intersection.

# Arguments

- `bv::BoundingVolume`: bounding volume.
- `ball::Ball`: ball.
- `tol::Real`: tolerance used in geometric comparisons.

# Returns

A [`BoundingVolume`](@ref), possibly empty.

# Examples

```julia
using SearchableGeometries

bv = BoundingVolume([0.0, 0.0], [3.0, 3.0])
ball = Ball([1.0, 1.0], 1.0)
intersection_bv = get_intersection(bv, ball)
```

# See also

[`intersects`](@ref), [`tighten_bv_bounds!`](@ref)
"""
function get_intersection(bv::BoundingVolume, ball::Ball; tol::Real=DEFAULT_BV_POINT_TOL)
    if !intersects(bv, ball; include_boundary=true, tol=tol)
        return BoundingVolume()
    end

    bv_ball = BoundingVolume(ball; tol=tol)
    cropped_bv = get_intersection(bv, bv_ball, tol=tol)

    # Check if the ball's center is in the BV or the BV is completely contained in the ball
    if is_contained(ball, cropped_bv)
        return cropped_bv
    end

    # The ball's center is not contained in the BV, so it
    # may be possible to crop the BV further
    tighten_bv_bounds!(cropped_bv, ball, tol=tol)
    return cropped_bv
end


# Hyperplanes ----------------------------------------------------------------------

@doc raw"""
    Hyperplane(point, n)

Construct an affine hyperplane from a point and a normal vector.

The hyperplane is the set

```math
\{x \in \mathbb{R}^n : n^T(x - point) = 0\}.
```

The constructor normalizes the normal vector, so the stored field `n` has unit
Euclidean norm. This makes signed offsets `dot(n, x - point)` equal to signed
Euclidean distances from the hyperplane.

# Arguments

- `point::Vector`: any point on the hyperplane.
- `n::Vector`: nonzero normal vector.

# Fields

- `point`: a point on the hyperplane.
- `n`: unit normal vector.
- `dim`: number of coordinates with nonzero normal entries.
- `embedding_dim`: dimension of the ambient coordinate space.
- `active_dim`: coordinates that affect the hyperplane equation.
- `inactive_dim`: coordinates with zero normal coefficient.
- `is_active`: Boolean mask indicating nonzero normal coefficients.

# Throws

Throws an error if `point` and `n` have different lengths or if `n` is the zero
vector.

# Examples

```julia
using SearchableGeometries

# The line x - y = 0 in R².
plane = Hyperplane([0.0, 0.0], [1.0, -1.0])

# The plane z = 2 in R³.
zplane = Hyperplane([0.0, 0.0, 2.0], [0.0, 0.0, 1.0])
```

# See also

[`is_contained`](@ref), [`intersects`](@ref), [`get_closest_point`](@ref)
"""
struct Hyperplane <: SearchableGeometry
    point::Vector
    n::Vector
    dim::Integer
    embedding_dim::Integer
    active_dim::Vector
    inactive_dim::Vector
    is_active::Vector{Bool}

    function Hyperplane(point::Vector, n::Vector)
        if length(point) != length(n)
            throw("SearchableGeometries.Hyperplane: point and normal vector must have the same dimension")
        end

        if iszero(norm(n))
            throw("SearchableGeometries.Hyperplane: normal vector must be nonzero")
        end

        all_indices = [eachindex(point)...]
        is_active = n .!= 0
        return new(point, n ./ norm(n), sum(is_active), length(point), all_indices[is_active], all_indices[.!is_active], is_active)
    end
end

import Base.==
function Base.:(==)(plane1::Hyperplane, plane2::Hyperplane; tol::Real=DEFAULT_BV_POINT_TOL)
    return isapprox(plane1.n' * plane1.point, plane2.n' * plane2.point, atol=tol) &&
           all(plane1.n .== plane2.n) &&
           plane1.dim == plane2.dim &&
           plane1.embedding_dim == plane2.embedding_dim &&
           all(plane1.active_dim .== plane2.active_dim) &&
           all(plane1.inactive_dim .== plane2.inactive_dim) &&
           all(plane1.is_active .== plane2.is_active)
end

@doc raw"""
    is_contained(plane::Hyperplane, query_pt; tol=DEFAULT_BV_POINT_TOL)

Test whether `query_pt` lies on `plane`.

A point is contained in a hyperplane if its signed offset from the plane is zero
up to tolerance:

```math
|n^T(query\_pt - point)| \le tol.
```

Because the constructor normalizes `n`, this signed offset is a signed Euclidean
distance.

# Arguments

- `plane::Hyperplane`: hyperplane.
- `query_pt::Vector{<:Real}`: point to test.
- `tol::Real`: tolerance for deciding whether the point lies on the plane.

# Returns

`true` if the point lies on the hyperplane up to `tol`; otherwise `false`.

# Throws

Throws an error if the point dimension does not match `plane.embedding_dim`.

# Examples

```julia
using SearchableGeometries

plane = Hyperplane([0.0, 0.0], [1.0, -1.0])
is_contained(plane, [2.0, 2.0]) # true
is_contained(plane, [2.0, 0.0]) # false
```

# See also

[`Hyperplane`](@ref), [`get_closest_point`](@ref)
"""
function is_contained(plane::Hyperplane, query_pt::Vector{<:Real}; tol::Real=DEFAULT_BV_POINT_TOL)
    if length(query_pt) != plane.embedding_dim
        throw("SearchableGeometries.Hyperplane: point dimension($(length(query_pt))) does not match hyperplane embedding dimension($(plane.embedding_dim))")
    end

    return abs(dot(plane.n, query_pt - plane.point)) <= tol
end

@doc raw"""
    get_closest_point(pt, query_plane::Hyperplane)

Project a point onto a hyperplane.

For a hyperplane with unit normal `n` and point `p0`, the closest point is

```math
pt - n^T(pt - p0)n.
```

This is the orthogonal projection of `pt` onto `query_plane`.

# Arguments

- `pt::Vector{<:Real}`: point to project.
- `query_plane::Hyperplane`: target hyperplane.

# Returns

The closest point on `query_plane` to `pt`.

# Throws

Throws an error if the point dimension does not match the hyperplane embedding
dimension.

# Examples

```julia
using SearchableGeometries

plane = Hyperplane([0.0, 0.0], [1.0, -1.0])
get_closest_point([2.0, 0.0], plane) # approximately [1.0, 1.0]
```

# See also

[`Hyperplane`](@ref), [`is_contained`](@ref)
"""
function get_closest_point(pt::Vector{<:Real}, query_plane::Hyperplane)
    if length(pt) != query_plane.embedding_dim
        throw("SearchableGeometries.Hyperplane: point dimension($(length(pt))) does not match hyperplane embedding dimension($(query_plane.embedding_dim))")
    end

    return pt - dot(query_plane.n, pt - query_plane.point) * query_plane.n
end

"""
    get_furthest_pt(bv::BoundingVolume, n; tol=DEFAULT_BV_POINT_TOL)

Return a point in `bv` maximizing `dot(n, x)`.

For each coordinate, this method chooses the upper bound when `n[i]` is positive
and the lower bound when `n[i]` is negative. Coordinates with approximately zero
coefficient do not affect the dot product; the implementation uses a
deterministic tie-breaking rule.

This helper is useful for computing extrema of linear functions over an
axis-aligned bounding volume.

# Arguments

- `bv::BoundingVolume`: bounding volume to search.
- `n::AbstractVector{<:Real}`: nonzero direction vector.
- `tol::Real`: tolerance used to treat small coefficients as zero.

# Returns

A point in `bv` maximizing `dot(n, x)`.

# Throws

Throws an error if the dimension of `n` does not match `bv`, or if `n` is zero.

# Examples

```julia
using SearchableGeometries

bv = BoundingVolume([0.0, 0.0], [2.0, 3.0])
get_furthest_pt(bv, [1.0, -1.0]) # [2.0, 0.0]
```

# See also

[`get_antifurthest_pt`](@ref), [`intersects`](@ref), [`Hyperplane`](@ref)
"""
function get_furthest_pt(bv::BoundingVolume, n::Vector{<:Real})
    if length(n) != length(bv.lb)
        throw("SearchableGeometries.get_furthest_pt: normal vector dimension($(length(n))) does not match bounding volume dimension($(length(bv.lb)))")
    end

    T = promote_type(eltype(bv.lb), eltype(bv.ub), eltype(n))
    pt = similar(bv.lb, T)

    for i in eachindex(n)
        pt[i] = n[i] >= 0 ? bv.ub[i] : bv.lb[i]
    end

    return pt
end

"""
    get_antifurthest_pt(bv::BoundingVolume, n; tol=DEFAULT_BV_POINT_TOL)

Return a point in `bv` minimizing `dot(n, x)`.

This is the opposite linear extreme from [`get_furthest_pt`](@ref). For each
coordinate, it chooses the lower bound when `n[i]` is positive and the upper
bound when `n[i]` is negative.

# Arguments

- `bv::BoundingVolume`: bounding volume to search.
- `n::AbstractVector{<:Real}`: direction vector.
- `tol::Real`: tolerance used to treat small coefficients as zero.

# Returns

A point in `bv` minimizing `dot(n, x)`.

# Examples

```julia
using SearchableGeometries

bv = BoundingVolume([0.0, 0.0], [2.0, 3.0])
get_antifurthest_pt(bv, [1.0, -1.0]) # [0.0, 3.0]
```

# See also

[`get_furthest_pt`](@ref), [`intersects`](@ref), [`Hyperplane`](@ref)
"""
function get_antifurthest_pt(bv::BoundingVolume, n::Vector{<:Real})
    return get_furthest_pt(bv, -n)
end

"""
    intersects(bv::BoundingVolume, query_plane::Hyperplane; include_boundary=true, tol=DEFAULT_BV_POINT_TOL)

Test whether a bounding volume intersects a hyperplane.

The method evaluates the minimum and maximum signed offsets of the bounding
volume from the hyperplane. If zero lies between these two extrema, then the
hyperplane intersects the bounding volume.

With `include_boundary=true`, boundary contact counts as intersection. With
`include_boundary=false`, the hyperplane must pass through the strict interior
of the bounding volume.

# Arguments

- `bv::BoundingVolume`: bounding volume.
- `query_plane::Hyperplane`: hyperplane.
- `include_boundary::Bool`: whether boundary-only contact counts as intersection.
- `tol::Real`: tolerance used in signed-offset comparisons.

# Returns

`true` if the bounding volume intersects the hyperplane; otherwise `false`.

# Throws

Throws an error if the bounding volume dimension does not match the hyperplane
embedding dimension.

# Examples

```julia
using SearchableGeometries

bv = BoundingVolume([0.0, 0.0], [1.0, 1.0])
plane = Hyperplane([0.5, 0.0], [1.0, 0.0])
intersects(bv, plane) # true
```

# See also

[`get_intersection`](@ref), [`get_furthest_pt`](@ref), [`get_antifurthest_pt`](@ref)
"""
function intersects(bv::BoundingVolume, query_plane::Hyperplane; include_boundary::Bool=true, tol::Real=DEFAULT_BV_POINT_TOL)
    # If the bounding volume is empty, it cannot intersect with anything
    if bv.is_empty
        return false
    end

    # The dimension of the bounding volume and the hyperplane must match
    if length(bv.lb) != query_plane.embedding_dim
        throw("SearchableGeometries.Hyperplane: bounding volume dimension($(length(bv.lb))) does not match hyperplane embedding dimension($(query_plane.embedding_dim))")
    end

    pt = query_plane.point
    n = query_plane.n

    # Compute the furthest and antifurthest points
    furthest_pt = get_furthest_pt(bv, n)
    antifurthest_pt = get_antifurthest_pt(bv, n)

    # Compute the minimum and maximum signed offsets over the BV
    smin = dot(n, antifurthest_pt - pt)
    smax = dot(n, furthest_pt - pt)

    if include_boundary
        # The hyperplane intersects the bounding volume if 0 is in the interval [smin, smax]
        return smin <= tol && smax >= -tol
    else
        # The hyperplane intersects the bounding volume if 0 is in the open interval (smin, smax)
        return smin < -tol && smax > tol
    end
end

"""
    tighten_bv_bounds(bv::BoundingVolume, query_plane::Hyperplane; tol=DEFAULT_BV_POINT_TOL)

Return a new bounding volume tightened around `bv ∩ query_plane`.

This is a non-mutating operation. The returned bounding volume is the smallest
axis-aligned box obtained by tightening each coordinate interval using the plane
equation and the remaining coordinate bounds.

If the hyperplane does not intersect `bv`, this method throws an error because
there is no nonempty tightened bounding volume representing the intersection.
Use [`get_intersection`](@ref) when you want a non-throwing method that returns
`BoundingVolume()` for empty intersections.

# Arguments

- `bv::BoundingVolume`: original bounding volume.
- `query_plane::Hyperplane`: hyperplane used for tightening.
- `tol::Real`: tolerance used in intersection and bound checks.

# Returns

A new [`BoundingVolume`](@ref) enclosing `bv ∩ query_plane`.

# Throws

Throws an error if dimensions do not match or if `bv` does not intersect the
hyperplane.

# Examples

```julia
using SearchableGeometries

bv = BoundingVolume([0.0, 0.0], [2.0, 2.0])
plane = Hyperplane([1.0, 0.0], [1.0, 1.0]) # x + y = 1
tight = tighten_bv_bounds(bv, plane)
```

# See also

[`get_intersection`](@ref), [`intersects`](@ref), [`Hyperplane`](@ref)
"""
function tighten_bv_bounds(bv::BoundingVolume, query_plane::Hyperplane; tol=DEFAULT_BV_POINT_TOL::Real)
    # Empty BV: return an empty/tightened copy
    if bv.is_empty
        return BoundingVolume(copy(bv.lb), copy(bv.ub))
    end

    # The BV and hyperplane must live in the same embedding space
    if length(bv.lb) != query_plane.embedding_dim
        throw("SearchableGeometries.Hyperplane: bounding volume dimension($(length(bv.lb))) does not match hyperplane embedding dimension($(query_plane.embedding_dim))")
    end

    # If the BV does not intersect the hyperplane, there is no tightened BV
    # representing BV ∩ Hyperplane.
    if !intersects(bv, query_plane; include_boundary=true, tol=tol)
        throw("SearchableGeometries.Hyperplane: cannot tighten BoundingVolume because it does not intersect the hyperplane")
    end

    # Plane equation: dot(n, x) = c0
    T = promote_type(eltype(bv.lb), eltype(query_plane.point), eltype(query_plane.n), typeof(tol))

    n = query_plane.n
    is_active = query_plane.is_active
    c0 = dot(n, query_plane.point)

    new_lb = T.(copy(bv.lb))
    new_ub = T.(copy(bv.ub))

    # Compute the furthest and antifurthest points
    f_pt = get_furthest_pt(bv, n)
    af_pt = get_antifurthest_pt(bv, n)

    # Compute the minimum and maximum signed offsets over the BV
    total_min = dot(n, af_pt)
    total_max = dot(n, f_pt)

    for i in eachindex(bv.lb)
        # If the plane does not depend on coordinate i, skip
        if !is_active[i]
            continue
        end

        ni = n[i]

        # Compute min/max contribution from all coordinates except i
        rest_min = total_min - ni * af_pt[i]
        rest_max = total_max - ni * f_pt[i]

        # We want values of x_i = t such that:
        #   ni * t + rest = c0
        # where rest is in [rest_min, rest_max]
        t1 = (c0 - rest_max) / ni
        t2 = (c0 - rest_min) / ni

        t_lb = ni > 0 ? t1 : t2
        t_ub = ni > 0 ? t2 : t1

        new_lb[i] = max(bv.lb[i], t_lb)
        new_ub[i] = min(bv.ub[i], t_ub)
    end

    return BoundingVolume(new_lb, new_ub; tol=tol)
end

"""
    get_intersection(bv::BoundingVolume, query_plane::Hyperplane; tol=DEFAULT_BV_POINT_TOL)

Return a bounding volume enclosing `bv ∩ query_plane`.

If the bounding volume does not intersect the hyperplane, this method returns
`BoundingVolume()`. Otherwise, it returns the result of
[`tighten_bv_bounds`](@ref), which is a tightened axis-aligned enclosure of the
flat intersection.

# Arguments

- `bv::BoundingVolume`: bounding volume.
- `query_plane::Hyperplane`: hyperplane.
- `tol::Real`: tolerance used in intersection and tightening.

# Returns

A [`BoundingVolume`](@ref), possibly empty.

# Examples

```julia
using SearchableGeometries

bv = BoundingVolume([0.0, 0.0], [2.0, 2.0])
plane = Hyperplane([1.0, 0.0], [1.0, 1.0])
intersection_bv = get_intersection(bv, plane)
```

# See also

[`intersects`](@ref), [`tighten_bv_bounds`](@ref)
"""
function get_intersection(bv::BoundingVolume, query_plane::Hyperplane; tol::Real=DEFAULT_BV_POINT_TOL)
    # No intersection at all
    if !intersects(bv, query_plane; include_boundary=true, tol=tol)
        return BoundingVolume()
    end

    # Return the tightened BV enclosing BV ∩ query_plane.
    # tighten_bv_bounds is non-mutating, so it returns a new BoundingVolume.
    return tighten_bv_bounds(bv, query_plane; tol=tol)
end

"""
    get_closest_point(bv::BoundingVolume, query_plane::Hyperplane; tol=DEFAULT_BV_POINT_TOL)

Return a point in `bv` closest to `query_plane`.

If `bv` intersects the hyperplane, the closest distance is zero, and this method
returns a deterministic feasible point in `bv ∩ query_plane`. If the bounding
volume lies entirely on one side of the hyperplane, the method returns the
corner of `bv` closest to the hyperplane.

The returned point is chosen deterministically. In tie cases, free coordinates
are selected using a lexicographic rule.

# Arguments

- `bv::BoundingVolume`: bounding volume to search.
- `query_plane::Hyperplane`: hyperplane.
- `tol::Real`: tolerance used in signed-distance and tie checks.

# Returns

A point in `bv` closest to `query_plane`.

# Throws

Throws an error if `bv` is empty or if dimensions do not match.

# Examples

```julia
using SearchableGeometries

bv = BoundingVolume([2.0, -1.0], [3.0, 1.0])
plane = Hyperplane([0.0, 0.0], [1.0, 0.0]) # x = 0
get_closest_point(bv, plane) # [2.0, -1.0] with lexicographic tie-breaking
```

# See also

[`get_furthest_point`](@ref), [`intersects`](@ref), [`Hyperplane`](@ref)
"""
function get_closest_point(bv::BoundingVolume, query_plane::Hyperplane; tol=DEFAULT_BV_POINT_TOL::Real)
    # If the bounding volume is empty, you cannot find a closest point
    if bv.is_empty
        throw("SearchableGeometries.Hyperplane: cannot compute closest point of an empty BoundingVolume")
    end

    # The dimension of the bounding volume and the hyperplane must match
    if length(bv.lb) != query_plane.embedding_dim
        throw("SearchableGeometries.Hyperplane: bounding volume dimension($(length(bv.lb))) does not match hyperplane embedding dimension($(query_plane.embedding_dim))")
    end

    n = query_plane.n
    pt = query_plane.point
    is_active = query_plane.is_active
    c = dot(n, pt)

    # Compute the furthest and antifurthest points
    f_pt = get_furthest_pt(bv, n)
    af_pt = get_antifurthest_pt(bv, n)

    # Compute the minimum and maximum signed offsets over the BV
    smin = dot(n, af_pt - pt)
    smax = dot(n, f_pt - pt)

    # Case 1: the BV intersects the plane.
    # Then the closest distance is 0, so we return the lexicographically
    # smallest point in BV ∩ plane.
    if smin <= tol && smax >= -tol # intersects(bv, query_plane; include_boundary=true, tol=tol)
        T = promote_type(eltype(bv.lb), eltype(query_plane.point), eltype(query_plane.n), typeof(tol))
        closest_pt = Vector{T}(undef, length(bv.lb))

        # Precompute min/max contribution from each coordinate.
        #
        # For active coordinates:
        #   min contribution is n[k] * af_pt[k]
        #   max contribution is n[k] * f_pt[k]
        #
        # For inactive coordinates:
        #   contribution is zero because the plane does not depend on them.
        d = length(bv.lb)

        min_contrib = zeros(T, d)
        max_contrib = zeros(T, d)

        for k in eachindex(bv.lb)
            if !is_active[k]
                continue
            end

            min_contrib[k] = n[k] * af_pt[k]
            max_contrib[k] = n[k] * f_pt[k]
        end

        # suffix_min[j] is the minimum contribution from coordinates j, ..., d.
        # suffix_max[j] is the maximum contribution from coordinates j, ..., d.
        suffix_min = zeros(T, d + 1)
        suffix_max = zeros(T, d + 1)

        for k in d:-1:1
            suffix_min[k] = suffix_min[k+1] + min_contrib[k]
            suffix_max[k] = suffix_max[k+1] + max_contrib[k]
        end

        # Contribution from coordinates that have already been fixed
        prefix_sum = zero(T)

        for j in eachindex(bv.lb)

            # Compute the smallest and largest possible contribution
            # from the remaining coordinates j+1, ..., end
            rem_min = suffix_min[j+1]
            rem_max = suffix_max[j+1]

            # If the plane does not depend on coordinate j, choose the smallest
            # allowed value to get the lexicographically smallest feasible point.
            if !is_active[j]
                closest_pt[j] = bv.lb[j]
            else
                nj = n[j]

                # Need:
                #   prefix_sum + nj*x[j] + remaining_contribution = c
                #
                # Since remaining_contribution can vary in [rem_min, rem_max],
                # this gives a feasible interval for x[j].
                rhs_low = c - prefix_sum - rem_max
                rhs_high = c - prefix_sum - rem_min

                t1 = rhs_low / nj
                t2 = rhs_high / nj

                t_lb = nj > 0 ? t1 : t2

                # Lexicographically smallest feasible choice
                closest_pt[j] = max(bv.lb[j], t_lb)
            end

            prefix_sum += n[j] * closest_pt[j]
        end

        return closest_pt
    end

    # Case 2: BV is entirely on the positive side of the plane.
    # Then the closest point is the point minimizing the signed offset.
    if smin > tol
        closest_pt = copy(af_pt)

        # Coordinates with n[i] ≈ 0 do not affect distance to the plane.
        # Choose lb to get the lexicographically smallest closest point.
        for i in eachindex(n)
            if !is_active[i]
                closest_pt[i] = bv.lb[i]
            end
        end

        return closest_pt
    end

    # Case 3: BV is entirely on the negative side of the plane.
    # Then the closest point is the point maximizing the signed offset.
    closest_pt = copy(f_pt)

    # Coordinates with n[i] ≈ 0 do not affect distance to the plane.
    # Choose lb to get the lexicographically smallest closest point.
    for i in eachindex(n)
        if !is_active[i]
            closest_pt[i] = bv.lb[i]
        end
    end

    return closest_pt
end

"""
    get_furthest_point(bv::BoundingVolume, query_plane::Hyperplane; tol=DEFAULT_BV_POINT_TOL)

Return a point in `bv` farthest from `query_plane`.

The method compares the signed-distance extrema of `bv` relative to the
hyperplane and returns the point with the largest absolute signed distance. If
both sides are equally far, a deterministic lexicographic tie-breaking rule is
used.

Because [`Hyperplane`](@ref) normalizes its normal vector, absolute signed
offsets are Euclidean distances to the hyperplane.

# Arguments

- `bv::BoundingVolume`: bounding volume to search.
- `query_plane::Hyperplane`: hyperplane.
- `tol::Real`: tolerance used in signed-distance and tie checks.

# Returns

A point in `bv` farthest from `query_plane`.

# Throws

Throws an error if `bv` is empty or if dimensions do not match.

# Examples

```julia
using SearchableGeometries

bv = BoundingVolume([0.0, 0.0], [2.0, 1.0])
plane = Hyperplane([0.0, 0.0], [1.0, 0.0]) # x = 0
get_furthest_point(bv, plane) # [2.0, 0.0] with lexicographic tie-breaking
```

# See also

[`get_closest_point`](@ref), [`get_furthest_pt`](@ref), [`get_antifurthest_pt`](@ref)
"""
function get_furthest_point(bv::BoundingVolume, query_plane::Hyperplane; tol=DEFAULT_BV_POINT_TOL::Real)
    # If the bounding volume is empty, you cannot find a furthest point
    if bv.is_empty
        throw("SearchableGeometries.Hyperplane: cannot compute furthest point of an empty BoundingVolume")
    end

    # The dimension of the bounding volume and the hyperplane must match
    if length(bv.lb) != query_plane.embedding_dim
        throw("SearchableGeometries.Hyperplane: bounding volume dimension($(length(bv.lb))) does not match hyperplane embedding dimension($(query_plane.embedding_dim))")
    end

    n = query_plane.n
    pt = query_plane.point
    is_active = query_plane.is_active

    # Compute the furthest and antifurthest points
    f_pt = get_furthest_pt(bv, n)
    af_pt = get_antifurthest_pt(bv, n)

    # Compute the minimum and maximum signed offsets over the BV
    smin = dot(n, af_pt - pt)
    smax = dot(n, f_pt - pt)

    # Case 1: If the antifurthest point is strictly farther from the plane,
    # return it as the furthest point from the plane.
    if abs(smin) > abs(smax) + tol
        furthest_pt = copy(af_pt)

        # Coordinates with n[i] == 0 do not affect distance to the plane.
        # Choose lb to get the lexicographically smallest furthest point.
        for i in eachindex(n)
            if !is_active[i]
                furthest_pt[i] = bv.lb[i]
            end
        end

        return furthest_pt
    end

    # Case 2: If the furthest point in the normal direction is strictly farther
    # from the plane, return it.
    if abs(smax) > abs(smin) + tol
        furthest_pt = copy(f_pt)

        # Coordinates with n[i] == 0 do not affect distance to the plane.
        # Choose lb to get the lexicographically smallest furthest point.
        for i in eachindex(n)
            if !is_active[i]
                furthest_pt[i] = bv.lb[i]
            end
        end

        return furthest_pt
    end

    # Case 3: Tie case - both sides are equally far from the plane.
    # First make both candidates lexicographically smallest in free coordinates.
    af_candidate = copy(af_pt)
    f_candidate = copy(f_pt)

    for i in eachindex(n)
        if !is_active[i]
            af_candidate[i] = bv.lb[i]
            f_candidate[i] = bv.lb[i]
        end
    end

    # Break ties lexicographically.
    for i in eachindex(af_candidate)
        if af_candidate[i] < f_candidate[i] - tol
            return af_candidate
        elseif af_candidate[i] > f_candidate[i] + tol
            return f_candidate
        end
    end

    # Equal up to tolerance: either is fine
    return af_candidate
end

end # module SearchableGeometries