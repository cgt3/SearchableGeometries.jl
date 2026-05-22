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

[`BoundingVolume`](@ref), [`Hyperplane`](@ref)
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

[`is_contained(::BoundingVolume, ::Ball)`](@ref), [`is_contained(::Ball, ::BoundingVolume)`](@ref)
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

[`is_contained(::Ball, ::Vector{<:Real})`](@ref), [`is_contained(::Ball, ::BoundingVolume)`](@ref)
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

[`is_contained(::Ball, ::Vector{<:Real})`](@ref), [`is_contained(::BoundingVolume, ::Ball)`](@ref)
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

[`get_intersection(::BoundingVolume, ::Ball)`](@ref)
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

[`tighten_bv_bounds!`](@ref)
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

[`get_reduced_dim_ball`](@ref)
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

[`intersects(::BoundingVolume, ::Ball)`](@ref)
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