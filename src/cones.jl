struct Cone <: SearchableGeometry
    vertex::Vector
    axis::Vector
    slope::Real

    function Cone(vertex::Vector{<:Real}, axis::Vector{<:Real}, slope::Real)
        if length(vertex) != length(axis)
            throw("SearchableGeometries.GeometricPrimitives.Cone: Vertex and axis vector do not have the same dimensions (dim(vertex)=$(length(vertex)), dim(axis)=$(length(axis))")
        elseif slope < zero(slope)
            throw("SearchableGeometries.GeometricPrimitives.Cone: Cannot construct cone with negative slope.")
        end

        axis_norm = norm(axis)
        if iszero(axis_norm)
            throw("SearchableGeometries.GeometricPrimitives.Cone: Cannot construct cone with zero axis vector.")
        end

        return new(vertex, axis ./ axis_norm, slope)
    end
end

function Base.:(==)(cone1::Cone, cone2::Cone; tol=DEFAULT_BV_POINT_TOL::Real)
    return all(cone1.vertex .== cone2.vertex) &&
           all(cone1.axis .== cone2.axis) &&
           cone1.slope == cone2.slope
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

