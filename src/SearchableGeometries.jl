module SearchableGeometries

    using LinearAlgebra

    # Data types
    export SearchableGeometry, Ball, BoundingVolume
    
    # BV only functions
    export getClosestPoint, getFurthestPoint, faceIndex2SpatialIndex, getFaceBoundingVolume
    
    # General functions
    export isContained, intersects, getIntersection

    import Base.getindex

    const DEFAULT_BV_POINT_TOL = 1e-15

    abstract type SearchableGeometry end

    # Balls ---------------------------------------------------------------------------------

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
            center::Vector{<:Real}, radius::Real; p=2::Real, 
            active_dim=true::Bool, indices=(active_dim ? [eachindex(center)...] : Vector{Int}[])::Vector{Int}
        )
            if radius < 0
                throw("SearchableGeometries.Ball: Cannot construct ball with negative radius.")
            elseif radius == 0
                return new(center, radius, p, 0, [], [eachindex(center)...], zeros(Bool, length(center)), length(center))
            end

            unique_indices = unique(indices)
            if active_dim
                is_active = zeros(Bool, length(center))
                is_active[unique_indices] .= true
                
                all_dim = [eachindex(center)...]
                inactive_dim = all_dim[is_active .== false]
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

    # Bounding Volumes ----------------------------------------------------------------------

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

    function getClosestPoint(bv::BoundingVolume, query_pt::Array)
        closest_pt = copy(query_pt)

        I_lb = query_pt .< bv.lb
        I_ub = query_pt .> bv.ub

        closest_pt[I_lb] .= bv.lb[I_lb]
        closest_pt[I_ub] .= bv.ub[I_ub]

        return closest_pt
    end

    function getFurthestPoint(bv::BoundingVolume, query_pt::Array)
        furthest_pt = similar(query_pt)

        ub_is_closer = 0.5*(bv.ub + bv.lb) .<= query_pt
        lb_is_closer = ub_is_closer .== false

        furthest_pt[ub_is_closer] .= bv.lb[ub_is_closer]
        furthest_pt[lb_is_closer] .= bv.ub[lb_is_closer]

        return furthest_pt
    end

    function isContained(bv::BoundingVolume, query_pt::Array; include_boundary=true::Bool)
        if (include_boundary && all(bv.lb .<= query_pt .<= bv.ub)) ||
            (!include_boundary && all(bv.lb .< query_pt .< bv.ub))
            return true
        else
            return false
        end
    end

    function isContained(bv::BoundingVolume, query_bv::BoundingVolume; include_boundary=true::Bool)
        if (!include_boundary && (all(query_bv.ub .< bv.ub) && all(query_bv.lb .> bv.lb))) ||
            (include_boundary && (all(query_bv.ub .<= bv.ub) && all(query_bv.lb .>= bv.lb)))
            return true
        else
            return false
        end
    end

    function intersects(bv1::BoundingVolume, bv2::BoundingVolume; include_boundary=true::Bool)
        if (include_boundary && (any(bv1.lb .> bv2.ub) || any(bv1.ub .< bv2.lb))) ||
            (!include_boundary && (any(bv1.lb .>= bv2.ub) || any(bv1.ub .<= bv2.lb)))
                return false
        else
            return true
        end
    end

    function getIntersection(bv1::BoundingVolume, bv2::BoundingVolume; tol=DEFAULT_BV_POINT_TOL::Real)
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

    function getFaceBoundingVolume(face_index::Integer, bv::BoundingVolume; tol=DEFAULT_BV_POINT_TOL::Real)
        face_lb, face_ub = copy(bv.lb), copy(bv.ub)
        
        if face_index <= length(bv.lb) # Lower bound face
            face_ub[face_index] = bv.lb[face_index]
        else # Upper bound face
            d = face_index - length(bv.lb)
            face_lb[d] = bv.ub[d]
        end

        return BoundingVolume(face_lb, face_ub; tol=tol)
    end
end