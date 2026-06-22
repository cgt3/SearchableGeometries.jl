import .GeometricPrimitives: BoundingVolume

const DEFAULT_NUM_LEAF_PTS = 40
const DEFAULT_PT_TOL = 1e-12

# Constants
export DEFAULT_NUM_LEAF_PTS

# Data types
export PointRange, Node, Watertight, Arbitrary

# Functions
export get_split, partition_data!, partition_bv

struct PointRange{T_data, dim}
    data::Matrix{T_data}
    first::Int
    last::Int
    n::Int

    function PointRange(data::Matrix{T}, first::Integer; n::Integer) where {T<:Real}
        n_int = Int(n)
        if n < 0
            throw("SearchableGeometries.PointRange: n must be nonnegative")
        end

        first_int = Int(first)
        n_int = Int(n)
        D = size(data, 1)

        return new{T, D}(data, first_int, first_int + n_int - 1, n_int)
    end

    function PointRange(data::Matrix{T}, first::Integer, last::Integer) where {T<:Real}
        first_int = Int(first)
        last_int = Int(last)
        n_int = last_int - first_int + 1

        if n_int < 0
            throw("SearchableGeometries.PointRange: last must be at least first - 1")
        end
        D = size(data, 1)

        return new{T, D}(data, first_int, last_int, n_int)
    end
end

Base.length(r::PointRange) = r.n
Base.isvalid(r::PointRange) = r.n > 0

abstract type BVH end

struct Watertight <: BVH end
struct Arbitrary <: BVH end

mutable struct Node{T_v, T_BVH, T_x}
	parent::Union{Node{T_v, T_BVH, T_x}, Nothing}
	left::Union{Node{T_v, T_BVH, T_x}, Nothing}
    right::Union{Node{T_v, T_BVH, T_x}, Nothing}
    split_dim::Int
    split_val::T_x
    is_leaf::Bool
    val::T_v
    bv::BoundingVolume
    # depth::Integer

    function Node{T_v, T_BVH, T_x}(
        val::T_v,
        bv::BoundingVolume;
        is_leaf::Function,
        partition_data!::Function,
        partition_bv::Function,
        parent::Union{Node{T_v, T_BVH, T_x}, Nothing}=nothing
    ) where {T_v, T_BVH<:BVH, T_x<:Real}
        
        # Create current internal node
        node = new{T_v, T_BVH, T_x}(
            parent,
            nothing,
            nothing,
            0,
            zero(T_x),
            true,
            val,
            bv
        )

        # Base case
        if is_leaf(val)
            return node
        end
        
        # Partition data
        split_dim, split_val, val_L, val_R = partition_data!(node)
        node.split_dim, node.split_val = split_dim, split_val
        
        # Partition bounding volumes
        bv_L, bv_R = partition_bv(node, val_L, val_R)
        
        # If both children are invalid, keep the node as a leaf.
        if !isvalid(val_L) && !isvalid(val_R)
            return node
        end

        # Otherwise ths is an internal node
        node.is_leaf = false

        # Build left child only if valid
        if isvalid(val_L)
            node.left = Node{T_v, T_BVH, T_x}(
                val_L,
                bv_L;
                is_leaf = is_leaf,
                partition_data! = partition_data!,
                partition_bv = partition_bv,
                parent = node
            )
        end
        
        # Build right child only if valid
        if isvalid(val_R)
            node.right = Node{T_v, T_BVH, T_x}(
                val_R,
                bv_R;
                is_leaf = is_leaf,
                partition_data! = partition_data!,
                partition_bv = partition_bv,
                parent = node
            )
        end

        return node
    end
end

function Node(
    val::PointRange,
    bv::BoundingVolume,
    ::Type{Arbitrary};
    is_leaf = is_leaf,
    partition_data! = partition_data!,
    partition_bv = partition_bv,
    parent::Union{Node{PointRange, Arbitrary, T_x}, Nothing}=nothing
) where {T_x<:Real}
    
    T = float(promote_type(eltype(bv.lb), eltype(bv.ub)))

    return Node{PointRange, Arbitrary, T}(
        val,
        bv;
        is_leaf,
        partition_data!,
        partition_bv,
        parent
    )
end

function Node(
    val::PointRange,
    bv::BoundingVolume,
    ::Type{Watertight};
    is_leaf = is_leaf,
    partition_data! = partition_data!,
    partition_bv = partition_bv,
    parent::Union{Node{PointRange, Watertight, T_x}, Nothing}=nothing
) where {T_x<:Real}

    T = float(promote_type(eltype(bv.lb), eltype(bv.ub)))

    return Node{PointRange, Watertight, T}(
        val,
        bv;
        is_leaf = is_leaf,
        partition_data!,
        partition_bv,
        parent
    )
end


function is_leaf(val::PointRange; leafsize::Int=DEFAULT_NUM_LEAF_PTS)
    # Stop if this node is already small enough
    if val.n <= leafsize
        return true
    end

    return false
end

function partition_data!(node::Node{PointRange, T_BVH, T_x}) where {T_BVH<:BVH, T_x<:Real}
    val = node.val
    data = val.data

    if val.n <= 1
        throw("SearchableGeometries.partition_data!: cannot partition a PointRange with $(val.n) point(s).")
    end

    lo = val.first
    hi = val.last
    n = val.n

    # Compute centroid of points in PointRange
    T, D = eltype(data), size(data, 1)
    centroid = zeros(float(T), D)

    @inbounds for i in lo:hi
        centroid .+= @view data[:, i]
    end

    centroid ./= n

    # Try dimensions from widest BV side to shortest BV side
    bv_lengths = node.bv.ub - node.bv.lb
    split_dims = sortperm(bv_lengths; rev=true)

    for split_dim in split_dims
        split_val = centroid[split_dim]

        # Check if the split actually seperates the data
        n_left = 0

        @inbounds for i in lo:hi
            if data[split_dim, i] <= split_val
                n_left += 1
            end
        end

        # If all points go left or all points go right,
        # this split is not useful. Try another dimension.
        if n_left == 0 || n_left == n
            continue
        end

        # Now partition the data in-place using your two-pointer logic.
        i_start = lo
        i_end = hi

        while i_start < i_end
            is_valid_L = data[split_dim, i_start] <= split_val
            is_valid_R = data[split_dim, i_end] > split_val

            if !is_valid_L && !is_valid_R
                tmp = copy(@view data[:, i_start])
                data[:, i_start] .= @view data[:, i_end]
                data[:, i_end] .= tmp
            end

            if is_valid_L && !is_valid_R
                i_start += 1
            elseif !is_valid_L && is_valid_R
                i_end -= 1
            else
                i_start += 1
                i_end -= 1
            end
        end

        last_left = i_start

        if data[split_dim, last_left] > split_val
            last_left -= 1
        end

        val_L = PointRange(data, lo, last_left)
        val_R = PointRange(data, last_left + 1, hi)

        # Safety check: if the two-pointer partition somehow produced
        # an invalid child, try another dimension.
        if !isvalid(val_L) || !isvalid(val_R)
            continue
        end

        return split_dim, split_val, val_L, val_R
    end

    # If no dimension gives a useful split, return invalid children.
    # The Node constructor should keep the current node as a leaf.
    invalid_L = PointRange(data, lo, lo - 1)
    invalid_R = PointRange(data, lo, lo - 1)

    return 0, zero(float(T)), invalid_L, invalid_R
end

function BoundingVolume(point_range::PointRange{T, D}) where {T<:Real, D}
    if !isvalid(point_range)
        throw("SearchableGeometries.BoundingVolume: cannot build a BoundingVolume from an empty PointRange.")
    end

    data = point_range.data

    lb = Vector{T}(undef, D)
    ub = Vector{T}(undef, D)

    cols = point_range.first:point_range.last

    for d in 1:D
        d_lb, d_ub = extrema(@view data[d, cols])

        lb[d] = d_lb
        ub[d] = d_ub
    end

    return BoundingVolume(lb, ub)
end

function partition_bv(
    node::Node{PointRange, Arbitrary, T_x},
    val_L::PointRange,
    val_R::PointRange
) where {T_x<:Real}

    return BoundingVolume(val_L), BoundingVolume(val_R)
end

function partition_bv(
    node::Node{PointRange, Watertight, T_x},
    val_L::PointRange,
    val_R::PointRange
) where {T_x<:Real}
    lb_L = copy(node.bv.lb)
    ub_L = copy(node.bv.ub)

    lb_R = copy(node.bv.lb)
    ub_R = copy(node.bv.ub)

    ub_L[node.split_dim] = node.split_val
    lb_R[node.split_dim] = node.split_val

    return BoundingVolume(lb_L, ub_L), BoundingVolume(lb_R, ub_R)
end