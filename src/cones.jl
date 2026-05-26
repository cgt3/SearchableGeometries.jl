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