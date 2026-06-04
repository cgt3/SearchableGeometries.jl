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
            throw("SearchableGeometries.GeometricPrimitives.BoundingVolume: lb (length: $(length(lb))) and ub (length: $(length(ub))) have different dimensions.")
        elseif any(lb .> ub)
            throw("SearchableGeometries.GeometricPrimitives.BoundingVolume: Cannot construct bounding volume with lb (=$lb) > ub (=$ub).")
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

function intersects(bv::BoundingVolume, line::Line; include_boundary=true, tol::Real=DEFAULT_BV_POINT_TOL)
    if length(bv.lb) != length(line.source)
        throw("SearchableGeometries.GeometricPrimitives.Line: bounding volume dimension($(length(bv.lb))) does not match line dimension($(length(line.source)))")
    end
    
    # If the line is a point
    if all(abs.(line.dir) .<= tol)
        return is_contained(bv, line.source; include_boundary=include_boundary)
    end

    s_min = include_boundary ? 0.0 : tol
    s_max = Inf

    for i in eachindex(bv.lb)
        if abs(line.dir[i]) <= tol
            # Line is parallel to this coordinate slab.
            if (include_boundary && (line.source[i] < bv.lb[i] - tol || line.source[i] > bv.ub[i] + tol)) ||
                (!include_boundary && (line.source[i] <= bv.lb[i] + tol || line.source[i] >= bv.ub[i] - tol))
                return false
            end
        else
            s1 = (bv.lb[i] - line.source[i]) / line.dir[i]
            s2 = (bv.ub[i] - line.source[i]) / line.dir[i]

            s_min = max(s_min, min(s1, s2))
            s_max = min(s_max, max(s1, s2))

            if (include_boundary && s_min > s_max + tol) ||
               (!include_boundary && s_min >= s_max - tol)
                return false
            end
        end
    end
    return true
end