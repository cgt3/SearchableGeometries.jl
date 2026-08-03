@doc raw""" 
    Cone(vertex, axis, slope)

Construct a one-sided Euclidean cone. 

The constructor normalizes `axis`. 

For a point `x`, define 
```math 
R(x) = a^\mathsf{T}(x-v)
``` 

and 

```math 
r(x) = \left\|x-v-R(x)a\right\|_2,
```

where `v` is the cone vertex and `a` is the normalized cone axis. 

The cone contains points satisfying 
```math 
R(x) \geq 0 
``` 

and 

```math 
r(x) \leq slope\,R(x).
``` 

# Arguments 

- `vertex::Vector{<:Real}`: Vertex of the cone. 
- `axis::Vector{<:Real}`: Nonzero direction in which the cone opens. 
- `slope::Real`: Nonnegative rate at which the cone widens. 

# Fields 

- `vertex`: Cone vertex. 
- `axis`: Normalized cone axis. 
- `slope`: Cone slope. 

# Throws 

Throws an error when: 
- `vertex` and `axis` have different dimensions; 
- `axis` is the zero vector; 
- `slope` is negative. 

# Examples 

```julia
using SearchableGeometries.GeometricPrimitives

cone = Cone(
    [0.0, 0.0],
    [1.0, 0.0],
    1.0
)
```
# See also 

[`BoundingVolume`](@ref), [`is_contained`](@ref), [`intersects`](@ref), [`get_intersection`](@ref) 
"""
struct Cone <: SearchableGeometry
    vertex::Vector
    axis::Vector
    slope::Real

    function Cone(vertex::Vector{<:Real}, axis::Vector{<:Real}, slope::Real)
        if length(vertex) != length(axis)
            throw("SearchableGeometries.GeometricPrimitives.Cone: Vertex and axis vector do not have the same dimensions (dim(vertex)=$(length(vertex)), dim(axis)=$(length(axis))")
        elseif slope <= zero(slope)
            throw("SearchableGeometries.GeometricPrimitives.Cone: Cannot construct cone with non-positive  slope.")
        end

        axis_norm = norm(axis)
        if iszero(axis_norm)
            throw("SearchableGeometries.GeometricPrimitives.Cone: Cannot construct cone with zero axis vector.")
        end

        return new(vertex, axis ./ axis_norm, slope)
    end
end

function Base.:(==)(cone1::Cone, cone2::Cone; tol=DEFAULT_BV_POINT_TOL::Real)
    return all(abs.(cone1.vertex .- cone2.vertex) .< tol) &&
           all(abs.(cone1.axis .- cone2.axis) .< tol) &&
           abs(cone1.slope - cone2.slope) < tol
end

function get_R(cone::Cone, p::Vector{<:Real})
    return (p .- cone.vertex)' * cone.axis
end

function get_radii(cone::Cone, p::Vector{<:Real})
    dist2Vertex = p .- cone.vertex
    R = dist2Vertex' * cone.axis
    r = norm(dist2Vertex .- R * cone.axis)
    return R, r
end

@doc raw""" 
    is_contained(cone::Cone, query_pt::Vector{<:Real}; include_boundary=true)

Test whether `query_pt` lies inside `cone`. 

A point is contained when: 

1. its axial distance from the cone vertex is nonnegative; and 
2. its perpendicular distance from the cone axis does not exceed the cone radius at that axial distance. 

Set `include_boundary=false` to exclude points on the cone surface and the cone vertex. 

The dimensions of `query_pt` and `cone` must match. 

# Arguments 

- `cone::Cone`: Cone used for the containment test. 
- `query_pt::Vector{<:Real}`: Point being tested. 
- `include_boundary::Bool=true`: Whether points on the cone boundary are considered contained. 

# Returns 

Returns `true` when `query_pt` is contained in `cone`. Otherwise, returns `false`. 

# Examples 

```julia
using SearchableGeometries.GeometricPrimitives

cone = Cone(
    [0.0, 0.0],
    [1.0, 0.0],
    1.0
)
is_contained(cone, [2.0, 1.0]) 
# true

is_contained(cone, [2.0, 3.0]) 
# false
```
"""
function is_contained(cone::Cone, query_pt::Vector{<:Real}; include_boundary=true)
    R, r = get_radii(cone, query_pt)
    if include_boundary
        return R >= 0 ? r <= R * cone.slope : false
    else
        return R > 0 ? r < R * cone.slope : false
    end
end

function get_bound_lines(cone::Cone)
    if cone.slope == 0
        return Line(cone.vertex, cone.axis, normalize=false), Line(cone.vertex, cone.axis, normalize=false)
    end

    lb_scalings = similar(cone.vertex)
    ub_scalings = similar(cone.vertex)
    e_i = zeros(length(cone.vertex))

    for i in eachindex(cone.vertex)
        e_i[i] = 1.0
        if isapprox(abs(e_i' * cone.axis), 1.0)
            lb_scalings[i] = cone.axis[i]
            ub_scalings[i] = cone.axis[i]
        else
            u0 = cone.vertex + e_i
            R0, r = get_radii(cone, u0)

            R_bound = r / cone.slope
            u_ub_i = u0[i] + (R_bound - R0)*cone.axis[i]
            u_lb_i = (cone.vertex - e_i)[i] + (R_bound + R0)*cone.axis[i]

            lb_scalings[i] = ( u_lb_i - cone.vertex[i] ) / R_bound
            ub_scalings[i] = ( u_ub_i - cone.vertex[i] ) / R_bound
        end
        e_i[i] = 0.0
    end

    return Line(cone.vertex, lb_scalings; normalize=false), Line(cone.vertex, ub_scalings; normalize=false)
end

function get_bounding_radii(cone::Cone, query_bv::BoundingVolume)
    closest_R_pt  = get_antiextreme_point(query_bv, cone.axis)
    furthest_R_pt = get_extreme_point(query_bv, cone.axis)
    return get_R(cone, closest_R_pt), get_R(cone, furthest_R_pt)
end

function BoundingVolume(lb_func::Line, ub_func::Line, R1::Real, R2::Real; tol::Real=DEFAULT_BV_POINT_TOL)
    lb1 = lb_func(R1)
    ub1 = ub_func(R1)

    lb2 = lb_func(R2)
    ub2 = ub_func(R2)

    return BoundingVolume(min.(lb1, lb2), max.(ub1, ub2); tol=tol)
end

@doc raw""" 
    BoundingVolume(cone::Cone, R_min::Real, R_max::Real, tol=DEFAULT_BV_POINT_TOL)
    
Construct an axis-aligned bounding volume enclosing the section of `cone` between axial radii `R_min` and `R_max`. 

The axial radii must satisfy 

```math 
0 \leq R_{\min} \leq R_{\max}.
``` 

# Arguments 

- `cone::Cone`: Cone being bounded. 
- `R_min::Real`: Minimum axial radius of the truncated cone section. 
- `R_max::Real`: Maximum axial radius of the truncated cone section. 
- `tol::Real=DEFAULT_BV_POINT_TOL`: Numerical tolerance used by the construction. 

# Returns 

Returns a [`BoundingVolume`](@ref) enclosing the truncated cone section. 

# Throws 

Throws an error when: 
- `R_min` is negative; 
- `R_max` is smaller than `R_min`. 

# Examples 

```julia
using SearchableGeometries.GeometricPrimitives

cone = Cone(
    [0.0, 0.0],
    [1.0, 0.0],
    1.0
)
bv = BoundingVolume(
    cone,
    0.0,
    2.0
    )
``` 
"""
function BoundingVolume(cone::Cone, R_min::Real, R_max::Real, tol::Real=DEFAULT_BV_POINT_TOL)
    if R_max < R_min
        throw("SearchableGeometries.GeometricPrimitives.Cone: Cannot construct bounding volume with R_max < R_min (R_max=$(R_max), R_min=$(R_min))")
    end

    if R_min < 0
        throw("SearchableGeometries.GeometricPrimitives.Cone: R_min and R_max do not satisfy 0 <= R_min <= R_max (R_min=$(R_min), R_max=$(R_max))")
    end
    
    cone_lb, cone_ub = get_bound_lines(cone)
    return BoundingVolume(cone_lb, cone_ub, R_min, R_max; tol=tol)
end

function get_alpha(bv::BoundingVolume, vertex::Vector{<:Real}, axis::Vector{<:Real}; is_max::Bool=true)
    s = is_max ? 1.0 : -1.0
    e_i = zeros(length(vertex))
    p = 0.5 * (bv.lb + bv.ub)

    for i in eachindex(vertex)
        d_vec = p .- vertex
        R_vec = (d_vec' * axis) * axis
        r_vec = d_vec .- R_vec
        e_i[i] = 1.0

        if s *  (e_i' * r_vec) >= 0
            p[i] = bv.ub[i]

        else
            p[i] = bv.lb[i]
        end
        e_i[i] = 0.0
    end
    
    # Compute alpha from the selected point
    R = (p .- vertex)' * axis
    r = norm((p .- vertex) .- R * axis)

    return r / R
end

@doc raw""" 
    is_contained(cone::Cone, query_bv::BoundingVolume; include_boundary=true, tol=DEFAULT_BV_POINT_TOL,)
    
Test whether the complete bounding volume `query_bv` lies inside `cone`. 

The result is `true` only when every point in `query_bv` is contained in the cone. 

Set `include_boundary=false` to require strict containment. 

The dimensions of `cone` and `query_bv` must match. 

# Arguments 

- `cone::Cone`: Cone used for the containment test. 
- `query_bv::BoundingVolume`: Bounding volume being tested. 
- `include_boundary::Bool=true`: Whether the cone boundary is included. 
- `tol::Real=DEFAULT_BV_POINT_TOL`: Numerical tolerance. 

# Returns 

Returns `true` when `query_bv` is completely contained in `cone`. Otherwise, returns `false`. 

# See also 

[`is_contained(::Cone, ::Vector{<:Real})`](@ref), [`intersects(::Cone, ::BoundingVolume)`](@ref) 
"""
function is_contained(cone::Cone, query_bv::BoundingVolume; include_boundary::Bool=true, tol::Real=DEFAULT_BV_POINT_TOL)
    if length(query_bv.lb) != length(cone.vertex)
        throw("SearchableGeometries.GeometricPrimitives.Cone: bounding volume dimension($(length(query_bv.lb))) does not match cone dimension($(length(cone.vertex)))")
    end

    R_min, _ = get_bounding_radii(cone, query_bv)

    alpha_max = get_alpha(query_bv, cone.vertex, cone.axis; is_max=true)

    if include_boundary
        return R_min >= -tol && alpha_max <= cone.slope + tol
    else
        return R_min > tol && alpha_max < cone.slope - tol
    end
end

@doc raw""" 
    intersects(cone::Cone, bv::BoundingVolume; include_boundary=true, tol=DEFAULT_BV_POINT_TOL) 
    
Test whether `cone` and `bv` have a nonempty intersection. 

When `include_boundary=true`, boundary-only contact counts as an intersection. 

When `include_boundary=false`, the cone and bounding volume must have a strict interior intersection. 

The dimensions of `cone` and `bv` must match. 

# Arguments 

- `cone::Cone`: Cone used in the intersection test. 
- `bv::BoundingVolume`: Bounding volume used in the intersection test. 
- `include_boundary::Bool=true`: Whether boundary-only contact counts as an intersection. 
- `tol::Real=DEFAULT_BV_POINT_TOL`: Numerical tolerance. 

# Returns 

Returns `true` when `cone` and `bv` intersect. Otherwise, returns `false`. 

# See also 

[`get_intersection(::Cone, ::BoundingVolume)`](@ref), [`is_contained(::Cone, ::BoundingVolume)`](@ref) 
"""
function intersects(cone::Cone, bv::BoundingVolume; include_boundary::Bool=true, tol::Real=DEFAULT_BV_POINT_TOL)
    if length(bv.lb) != length(cone.vertex)
        throw("SearchableGeometries.GeometricPrimitives.Cone: bounding volume dimension($(length(bv.lb))) does not match cone dimension($(length(cone.vertex)))")
    end

    if intersects(bv, Line(cone.vertex, cone.axis; normalize=true); include_boundary=include_boundary, tol=tol)
        return true
    else
        _, R_max = get_bounding_radii(cone, bv)

        # If the whole BV is behind the cone vertex, it cannot intersect the cone.
        if include_boundary
            R_max < -tol && return false
        else
            R_max <= tol && return false
        end

        alpha_min = get_alpha(bv, cone.vertex, cone.axis; is_max=false)
        if include_boundary
            return alpha_min <= cone.slope + tol
        else
            return alpha_min < cone.slope - tol
        end
    end
end

@doc raw""" 
    get_intersection(cone::Cone, bv::BoundingVolume; tol=DEFAULT_BV_POINT_TOL, ) 
    
Return an axis-aligned bounding volume enclosing the intersection of `cone` and `bv`. 

The exact intersection may have curved or sloped boundaries and therefore may not itself be an axis-aligned bounding volume. 

This method returns a conservative bounding-volume enclosure of the exact intersection. 

If `cone` and `bv` do not intersect, the method returns an empty `BoundingVolume`. 

The dimensions of `cone` and `bv` must match. 

# Arguments 

- `cone::Cone`: Cone used in the intersection. 
- `bv::BoundingVolume`: Bounding volume used in the intersection. 
- `tol::Real=DEFAULT_BV_POINT_TOL`: Numerical tolerance. 

# Returns 

Returns a [`BoundingVolume`](@ref) enclosing the exact intersection. 

# See also 

[`intersects(::Cone, ::BoundingVolume)`](@ref), [`is_contained(::Cone, ::BoundingVolume)`](@ref)
"""
function get_intersection(cone::Cone, bv::BoundingVolume; tol::Real=DEFAULT_BV_POINT_TOL)
    if !intersects(cone, bv; include_boundary=true, tol=tol)
        return BoundingVolume()
    end
    
    intersection = bv
    cone_lbs, cone_ubs = get_bound_lines(cone)

    while !intersection.is_empty
        R_min, R_max = get_bounding_radii(cone, intersection)

        bv_cone = BoundingVolume(cone_lbs, cone_ubs, max(0.0, R_min), R_max; tol=tol)

        old_intersection = intersection
        intersection = get_intersection(old_intersection, bv_cone; tol=tol)

        if isapprox(old_intersection.lb, intersection.lb; atol=tol) &&
        isapprox(old_intersection.ub, intersection.ub; atol=tol)
            return intersection
        end
    end

    # If we reach here the intersection is empty, so we return the empty bounding volume.
    return BoundingVolume()
end