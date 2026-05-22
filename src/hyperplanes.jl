@doc raw"""
    Hyperplane(point, n)

Construct an affine hyperplane from a point `point` and a normal vector `n`.

The hyperplane is the set

```math
\{x \in \mathbb{R}^n : n^T(x - \operatorname{point}) = 0\}.
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

# The line x - y = 0 in 2D space.
plane = Hyperplane([0.0, 0.0], [1.0, -1.0])

# The plane z = 2 in 3D space.
zplane = Hyperplane([0.0, 0.0, 2.0], [0.0, 0.0, 1.0])
```

# See also

[`BoundingVolume`](@ref), [`Ball`](@ref)
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
|n^T(\operatorname{query\_pt} - \operatorname{point})| \le \operatorname{tol}.
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

[`get_closest_point(::Vector{<:Real}, ::Hyperplane)`](@ref)
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

For a hyperplane with unit normal `n` and point `p`, the closest point is

```math
pt - n^T(pt - p)n.
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

[`get_furthest_pt`](@ref), [`get_antifurthest_pt`](@ref)
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
- `n::Vector{<:Real}`: nonzero direction vector.
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

[`get_antifurthest_pt`](@ref)
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
- `n::Vector{<:Real}`: direction vector.
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

[`get_furthest_pt`](@ref)
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

[`get_intersection(::BoundingVolume, ::Hyperplane)`](@ref)
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

@doc raw"""
    tighten_bv_bounds(bv::BoundingVolume, query_plane::Hyperplane; tol=DEFAULT_BV_POINT_TOL)

Return a new bounding volume tightened around ``bv \cap query_plane``.

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

A new [`BoundingVolume`](@ref) enclosing `bv` intersecting with `query_plane`.

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

[`get_intersection(::BoundingVolume, ::Hyperplane)`](@ref)
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

@doc raw"""
    get_intersection(bv::BoundingVolume, query_plane::Hyperplane; tol=DEFAULT_BV_POINT_TOL)

Return a bounding volume enclosing ``\operatorname{bv} \cap \operatorname{query\_plane}``.

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

[`intersects(::BoundingVolume, ::Hyperplane)`](@ref), [`tighten_bv_bounds`](@ref)
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

@doc raw"""
    get_closest_point(bv::BoundingVolume, query_plane::Hyperplane; tol=DEFAULT_BV_POINT_TOL)

Return a point in `bv` closest to `query_plane`.

If `bv` intersects the hyperplane, the closest distance is zero, and this method
returns a deterministic feasible point in ``\operatorname{bv} \cap \operatorname{query\_plane}``. If the bounding
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

[`get_furthest_point(::BoundingVolume, ::Hyperplane)`](@ref)
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

[`get_closest_point(::BoundingVolume, ::Hyperplane)`](@ref)
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