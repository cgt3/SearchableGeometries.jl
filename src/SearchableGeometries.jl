module SearchableGeometries

    using LinearAlgebra

    # Data types
    export SearchableGeometry, BoundingVolume

    import Base.getindex

    const DEFAULT_BV_POINT_TOL = 1e-15

    abstract type SearchableGeometry end

    struct BoundingVolume <: SearchableGeometry
        lb::Vector
        ub::Vector
        is_empty::Bool
        dim::Integer
        active_dim::Vector
        inactive_dim::Vector
        is_active::Vector{Bool}
        
        # For invalid/empty Bounding Volumes (note an empty BV differs from a non-empty BV with dimension 0 (a point))
        function BoundingVolume()
            return new([Inf], [-Inf], true, 0, [], [], Bool[])
        end

        function BoundingVolume(
            lb::Vector{<:Real}, ub::Vector{<:Real}; tol=DEFAULT_BV_POINT_TOL::Real 
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
                    dim -=1
                    is_active[d] = false
                end
            end
            all_dim = [eachindex(lb)...]
            active_dim = all_dim[is_active]
            inactive_dim = all_dim[is_active .!= true]

            return new(lb, ub, false, dim, active_dim, inactive_dim, is_active)
        end
    end

end