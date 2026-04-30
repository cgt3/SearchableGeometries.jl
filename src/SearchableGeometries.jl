module SearchableGeometries

using LinearAlgebra

# Data Types:
export SearchableGeometry, Ball, BoundingVolume, Hyperplane

# BV only Functions:
export getClosestPoint, getFurthestPoint, faceIndex2SpatialIndex, getFaceBoundingVolume

# General Functions:
export isContained, intersects, getIntersection

# Ball only Functions:
export getReducedDimBall, tightenBVBounds!

# Hyperplane only Functions:

import Base.getindex

const DEFAULT_BV_POINT_TOL = 1e-15

"""
SearchableGeometry

Abstract supertype for searchable geometric objects supported by 'SearchableGeometries.jl'.
Concreate subtypes currently include `BoundingVolume` and `Ball`.
"""
abstract type SearchableGeometry end


# Bounding Volumes ----------------------------------------------------------------------
"""
BoundingVolume(lb, ub; tol=DEFAULT_BV_POINT_TOL)

Construct a bounding volume with lower bounds `lb` and upper bounds `ub`.

# Arguments
- `lb::Vector{<:Real}`: Lower bounds of the bounding volume.
- `ub::Vector{<:Real}`: Upper bounds of the bounding volume.
- `tol::Real`: Tolerance for determining active dimensions.

# Returns
- `BoundingVolume`: The constructed bounding volume.

# Throws
- `ArgumentError`: If `lb` and `ub` have different dimensions.
- `ArgumentError`: If `lb` > `ub`.

# Examples
```julia
using SearchableGeometries

bv = BoundingVolume([0, 0], [1, 1])
```
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

function getClosestPoint(bv::BoundingVolume, query_pt::Vector{<:Real})
    closest_pt = copy(query_pt)

    I_lb = query_pt .< bv.lb
    I_ub = query_pt .> bv.ub

    closest_pt[I_lb] .= bv.lb[I_lb]
    closest_pt[I_ub] .= bv.ub[I_ub]

    return closest_pt
end

function getFurthestPoint(bv::BoundingVolume, query_pt::Vector{<:Real})
    furthest_pt = similar(query_pt)

    ub_is_closer = 0.5 * (bv.ub + bv.lb) .<= query_pt
    lb_is_closer = ub_is_closer .== false

    furthest_pt[ub_is_closer] .= bv.lb[ub_is_closer]
    furthest_pt[lb_is_closer] .= bv.ub[lb_is_closer]

    return furthest_pt
end

function isContained(bv::BoundingVolume, query_pt::Vector{<:Real}; include_boundary::Bool=true)
    if (include_boundary && all(bv.lb .<= query_pt .<= bv.ub)) ||
       (!include_boundary && all(bv.lb .< query_pt .< bv.ub))
        return true
    else
        return false
    end
end

function isContained(bv::BoundingVolume, query_bv::BoundingVolume; include_boundary::Bool=true)
    if (!include_boundary && (all(query_bv.ub .< bv.ub) && all(query_bv.lb .> bv.lb))) ||
       (include_boundary && (all(query_bv.ub .<= bv.ub) && all(query_bv.lb .>= bv.lb)))
        return true
    else
        return false
    end
end

function intersects(bv1::BoundingVolume, bv2::BoundingVolume; include_boundary::Bool=true)
    if (include_boundary && (any(bv1.lb .> bv2.ub) || any(bv1.ub .< bv2.lb))) ||
       (!include_boundary && (any(bv1.lb .>= bv2.ub) || any(bv1.ub .<= bv2.lb)))
        return false
    else
        return true
    end
end

function getIntersection(bv1::BoundingVolume, bv2::BoundingVolume; tol::Real=DEFAULT_BV_POINT_TOL)
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

function faceIndex2SpatialIndex(face_index::Integer, num_dim::Integer)
    return face_index <= num_dim ? face_index : face_index - num_dim
end

function getFaceBoundingVolume(face_index::Integer, bv::BoundingVolume; tol::Real=DEFAULT_BV_POINT_TOL)
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

"""
Ball(center, radius; p=2, active_indices=true, indices)

Construct a ball with center `center` and radius `radius`.

# Arguments
- `center::Vector{<:Real}`: Center of the ball.
- `radius::Real`: Radius of the ball.
- `p::Real`: p-norm.
- `active_indices::Bool`: Whether to use active indices.
- `indices::Vector{Int}`: Indices of the active dimensions.

# Returns
- `Ball`: The constructed ball.

# Throws
- `ArgumentError`: If `radius` is negative.
- `ArgumentError`: If `p` is not a positive real number.

# Examples
```julia
using SearchableGeometries

ball = Ball([0, 0], 1)
```
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

function BoundingVolume(ball::Ball; tol::Real=DEFAULT_BV_POINT_TOL)
    lb = ball.center .- ball.radius
    ub = ball.center .+ ball.radius
    return BoundingVolume(lb, ub; tol=tol)
end

function isContained(ball::Ball, query_pt::Vector{<:Real}; include_boundary::Bool=true, tol::Real=DEFAULT_BV_POINT_TOL)
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

function isContained(bv::BoundingVolume, query_ball::Ball; include_boundary::Bool=true, tol::Real=DEFAULT_BV_POINT_TOL)
    return isContained(bv, BoundingVolume(query_ball; tol=tol); include_boundary=include_boundary)
end

function isContained(ball::Ball, query_bv::BoundingVolume; include_boundary::Bool=true)
    furthest_pt = getFurthestPoint(query_bv, ball.center)
    return isContained(ball, furthest_pt; include_boundary=include_boundary)
end

function intersects(bv::BoundingVolume, ball::Ball; include_boundary::Bool=true, tol::Real=DEFAULT_BV_POINT_TOL)
    # First, do the easy checks against the ball's BV:
    if !intersects(bv, BoundingVolume(ball; tol=tol); include_boundary=include_boundary)
        # The two are completely disjoint
        return false
    elseif isContained(ball, bv; include_boundary=include_boundary)
        # The ball is completely contained in the BV
        return true
    end
    # Note: the below check could catch all of the above cases but its use of the
    #       lp-norm makes it more expensive. The above checks are to avoid having
    #       to compute norms.

    # Check if the closest point on the BV to the ball's center is inside the ball
    closest_pt = getClosestPoint(bv, ball.center)
    return isContained(ball, closest_pt; include_boundary=include_boundary)
end

function getReducedDimBall(removal_dim::Integer, x_d::Real, ball::Ball)
    if x_d < ball.center[removal_dim] - ball.radius || ball.center[removal_dim] + ball.radius < x_d
        throw("SearchableGeometries.Ball: coordinate plane defined by x_$removal_dim = $x_d does not intersect the ball (center=$(ball.center), radius=$(ball.radius))")
    end

    new_center = copy(ball.center)
    new_center[removal_dim] = x_d

    new_radius = ball.p === Inf ? ball.radius : (ball.radius^ball.p - abs(x_d - ball.center[removal_dim])^ball.p)^(1 / ball.p)
    inactive_dim = [ball.inactive_dim..., removal_dim]

    return Ball(new_center, new_radius; p=ball.p, active_indices=false, indices=inactive_dim)
end

function tightenBVBounds!(bv::BoundingVolume, ball::Ball; tol::Real=DEFAULT_BV_POINT_TOL)
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
        face_bv = getFaceBoundingVolume(f_target, bv, tol=tol)

        non_simple = false
        if !intersects(face_bv, ball, include_boundary=true, tol=tol)
            adjacent_faces = [1:f_target-1..., f_target+1:2*num_dim...]
            d_target = faceIndex2SpatialIndex(f_target, num_dim)

            # Check if this face needs to be updated using a non-simple intersection
            if f_target <= num_dim # f_target is a lb face
                lb_pt_projected[d_target] = face_bv.lb[d_target]
                if isContained(face_bv, lb_pt_projected, include_boundary=true)
                    push!(altered_lb_indices, d_target)
                    bv.lb[d_target] = ball.center[d_target] - ball.radius
                    non_simple = true
                end
                lb_pt_projected[d_target] = ball.center[d_target]
            else # f_target is an ub face
                ub_pt_projected[d_target] = face_bv.ub[d_target]
                if isContained(face_bv, ub_pt_projected, include_boundary=true)
                    push!(altered_ub_indices, d_target)
                    bv.ub[d_target] = ball.center[d_target] + ball.radius
                    non_simple = true
                end
                ub_pt_projected[d_target] = ball.center[d_target]
            end

            # For simple intersections
            if !non_simple
                for f_adjacent in adjacent_faces
                    adjacent_face_bv = getFaceBoundingVolume(f_adjacent, bv, tol=tol)

                    if intersects(adjacent_face_bv, ball, include_boundary=true)
                        d_fixed = faceIndex2SpatialIndex(f_adjacent, num_dim)
                        reduced_ball = getReducedDimBall(d_fixed, adjacent_face_bv.lb[d_fixed], ball)

                        # This will modify face_adjacent's bounds 
                        altered_lb_indices_new, altered_ub_indices_new = tightenBVBounds!(adjacent_face_bv, reduced_ball, tol=tol)

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

function getIntersection(bv::BoundingVolume, ball::Ball; tol::Real=DEFAULT_BV_POINT_TOL)
    if !intersects(bv, ball; include_boundary=true, tol=tol)
        return BoundingVolume()
    end

    bv_ball = BoundingVolume(ball; tol=tol)
    cropped_bv = getIntersection(bv, bv_ball, tol=tol)

    # Check if the ball's center is in the BV or the BV is completely contained in the ball
    if isContained(ball, cropped_bv)
        return cropped_bv
    end

    # The ball's center is not contained in the BV, so it
    # may be possible to crop the BV further
    tightenBVBounds!(cropped_bv, ball, tol=tol)
    return cropped_bv
end


# Hyperplanes ----------------------------------------------------------------------
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
    return all(plane1.point .== plane2.point) &&
           all(plane1.n .== plane2.n) &&
           plane1.dim == plane2.dim &&
           plane1.embedding_dim == plane2.embedding_dim &&
           all(plane1.active_dim .== plane2.active_dim) &&
           all(plane1.inactive_dim .== plane2.inactive_dim) &&
           all(plane1.is_active .== plane2.is_active)
end

function isContained(plane::Hyperplane, query_pt::Vector{<:Real}; tol::Real=DEFAULT_BV_POINT_TOL)
    if length(query_pt) != plane.embedding_dim
        throw("SearchableGeometries.Hyperplane: point dimension($(length(query_pt))) does not match hyperplane embedding dimension($(plane.embedding_dim))")
    end

    return abs(dot(plane.n, query_pt - plane.point)) <= tol
end

function getClosestPoint(pt::Vector{<:Real}, query_plane::Hyperplane)
    if length(pt) != query_plane.embedding_dim
        throw("SearchableGeometries.Hyperplane: point dimension($(length(pt))) does not match hyperplane embedding dimension($(query_plane.embedding_dim))")
    end

    return pt - dot(query_plane.n, pt - query_plane.point) * query_plane.n
end

function signedExtrema(bv::BoundingVolume, query_plane::Hyperplane)
    # We study the signed offset function
    #     h(x) = dot(n, x - point)
    # over the whole BV. Since h is linear, its minimum and maximum
    # occur at corners of the box.
    T = promote_type(eltype(bv.lb), eltype(bv.ub), eltype(query_plane.point), eltype(query_plane.n))
    smin = zero(T)
    smax = zero(T)


    for d in eachindex(bv.lb)
        nd = query_plane.n[d]

        if nd > 0
            # Positive normal component:
            # - Lower bound gives the minimum signed value,
            # - Upper bound gives the maximum signed value.
            smin += nd * (bv.lb[d] - query_plane.point[d])
            smax += nd * (bv.ub[d] - query_plane.point[d])
        elseif nd < 0
            # Negative normal component:
            # - Lower bound gives the maximum signed value,
            # - Upper bound gives the minimum signed value.
            smin += nd * (bv.ub[d] - query_plane.point[d])
            smax += nd * (bv.lb[d] - query_plane.point[d])
        end
        # If nd == 0, this coordinate does not affect the signed offset
    end

    return smin, smax
end

function intersects(bv::BoundingVolume, query_plane::Hyperplane; include_boundary::Bool=true, tol::Real=DEFAULT_BV_POINT_TOL)
    # If the bounding volume is empty, it cannot intersect with anything
    if bv.is_empty
        return false
    end

    # The dimension of the bounding volume and the hyperplane must match
    if length(bv.lb) != query_plane.embedding_dim
        throw("SearchableGeometries.Hyperplane: bounding volume dimension($(length(bv.lb))) does not match hyperplane embedding dimension($(query_plane.embedding_dim))")
    end

    # Compute the minimum and maximum signed offsets over the BV
    smin, smax = signedExtrema(bv, query_plane)

    if include_boundary
        # The hyperplane intersects the bounding volume if 0 is in the interval [smin, smax]
        return smin <= tol && smax >= -tol
    else
        # The hyperplane intersects the bounding volume if 0 is in the open interval (smin, smax)
        return smin < -tol && smax > tol
    end
end

end # module SearchableGeometries