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

mutable struct Node{T_v, T_BVH, T_d}
	parent::Union{Node{T_v, T_BVH, T_d}, Nothing}
	left::Union{Node{T_v, T_BVH, T_d}, Nothing}
    right::Union{Node{T_v, T_BVH, T_d}, Nothing}
    split_dim::Int
    split_val::T_d
    is_leaf::Bool
    val::T_v
    bv::BoundingVolume
    # depth::Integer

    function Node{T_v, T_BVH, T_d}(
        val::T_v,
        bv::BoundingVolume;
        is_leaf::Function,
        partition_data!::Function,
        partition_bv::Function,
        parent::Union{Node{T_v, T_BVH, T_d}, Nothing}=nothing
    ) where {T_v, T_BVH<:BVH, T_d<:Real}
        
        # Create current internal node
        node = new{T_v, T_BVH, T_d}(
            parent,
            nothing,
            nothing,
            0,
            zero(T_d),
            true,
            val,
            bv
        )

        # Base case
        if is_leaf(val)
            return node
        end
        
        # Partition data
        node.split_dim, node.split_val, val_L, val_R = partition_data!(node)
        
        # Partition bounding volumes
        bv_L, bv_R = partition_bv(node, val_L, val_R)
        
        is_valid_L, is_valid_R = isvalid(val_L), isvalid(val_R)

        # If both children are invalid, keep the node as a leaf
        if !is_valid_L && !is_valid_R
            return node
        end

        # Otherwise ths is an internal node
        node.is_leaf = false

        # Build left child only if valid
        if is_valid_L
            node.left = Node{T_v, T_BVH, T_d}(
                val_L,
                bv_L;
                is_leaf = is_leaf,
                partition_data! = partition_data!,
                partition_bv = partition_bv,
                parent = node
            )
        end
        
        # Build right child only if valid
        if is_valid_R
            node.right = Node{T_v, T_BVH, T_d}(
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
    parent::Union{Node{PointRange, Arbitrary, T_d}, Nothing}=nothing
) where {T_d<:Real}

    T = eltype(bv.lb)

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
    parent::Union{Node{PointRange, Watertight, T_d}, Nothing}=nothing
) where {T_d<:Real}
    
    T = eltype(bv.lb)

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

function partition_data!(node::Node{PointRange, T_BVH, T_d}; tie_tol::Real=DEFAULT_PT_TOL) where {T_BVH<:BVH, T_d<:Real}
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
    centroid = zeros(T, D)

    @inbounds for i in lo:hi
        centroid .+= @view data[:, i]
    end

    centroid ./= n

    # Compute BV side lengths
    bv_lengths = node.bv.ub - node.bv.lb
    normalized_lengths = bv_lengths ./ maximum(bv_lengths)

    # Sort dimensions from the longest to shortest
    split_dims = sortperm(normalized_lengths; rev=true)

    # Keep only dimenions tied with the longest side
    tied_dims = [
        d for d in split_dims if normalized_lengths[d] >= one(T) - tie_tol
    ]

    # Compute n_left and n_right for each tied dimension
    n_lefts = Vector{Int}(undef, length(tied_dims))
    n_rights = Vector{Int}(undef, length(tied_dims))
    split_vals = Vector{T}(undef, length(tied_dims))

    for d in eachindex(tied_dims)
        split_dim = tied_dims[d]
        split_val = centroid[split_dim]

        n_left = 0

        @inbounds for i in lo:hi
            if data[split_dim, i] <= split_val
                n_left += 1
            end
        end
        
        n_right = n - n_left
        
        n_rights[d] = n_right
        n_lefts[d] = n_left
        split_vals[d] = split_val
    end

    best_d = argmin(abs.(n_lefts .- n_rights))
    split_dim = tied_dims[best_d]
    split_val = split_vals[best_d]
    n_left = n_lefts[best_d]

    # Now partition the data in-place using your the selected split_dim.
    i_start = lo
    i_end = hi

    while i_start < i_end
        is_valid_L = data[split_dim, i_start] <= split_val
        is_valid_R = data[split_dim, i_end] > split_val

        if !is_valid_L && !is_valid_R
            tmp = copy(@view data[:, i_start])
            data[:, i_start] .= @view data[:, i_end]
            data[:, i_end] .= tmp
        elseif is_valid_L && !is_valid_R
            i_start += 1
        elseif !is_valid_L && is_valid_R
            i_end -= 1
        else
            i_start += 1
            i_end -= 1
        end
    end

    # Since n_left was counted before partitioning, the left child
    # must contain exactly n_left points after partitioning.
    last_left = lo + n_left - 1

    val_L = PointRange(data, lo, last_left)
    val_R = PointRange(data, last_left + 1, hi)

    return split_dim, split_val, val_L, val_R
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
        lb[d], ub[d] = extrema(@view data[d, cols])
    end

    return BoundingVolume(lb, ub)
end

function partition_bv(
    node::Node{PointRange, Arbitrary, T_d},
    val_L::PointRange,
    val_R::PointRange
) where {T_d<:Real}

    return BoundingVolume(val_L), BoundingVolume(val_R)
end

function partition_bv(
    node::Node{PointRange, Watertight, T_d},
    val_L::PointRange,
    val_R::PointRange
) where {T_d<:Real}
    lb_L = copy(node.bv.lb)
    ub_L = copy(node.bv.ub)

    lb_R = copy(node.bv.lb)
    ub_R = copy(node.bv.ub)

    ub_L[node.split_dim] = node.split_val
    lb_R[node.split_dim] = node.split_val

    return BoundingVolume(lb_L, ub_L), BoundingVolume(lb_R, ub_R)
end

struct Tree{T_v, T_BVH, T_d}
    root::Node{T_v, T_BVH, T_d}

    function Tree{T_v, T_BVH, T_d}(val::T_v) where {T_v, T_BVH<:BVH, T_d<:Real}
        
        bv = BoundingVolume(val)
        
        root = Node{T_v, T_BVH, T_d}(val, bv)

        return new{T_v, T_BVH, T_d}(root)
    end
end

function tree_map_preorder(node::Node{T_v, T_BVH, T_d}, func::Function; left_first::Bool=true) where {T_v, T_BVH<:BVH, T_d<:Real}
    
    first_child, second_child = left_first ? (node.left, node.right) : (node.right, node.left)

    # If this is a leaf, there are no children to visit.
    if node.is_leaf
        return nothing
    end

    if left_first
        tree_map_preorder(node.left, func; left_first=left_first)
        tree_map_preorder(node.right, func; left_first=left_first)
    else
        tree_map_preorder(node.right, func; left_first=left_first)
        tree_map_preorder(node.left, func; left_first=left_first)
    end

    return nothing
end

function leaf_search(node::Node{T_v, T_BVH, T_d}, func::Function) where {T_v, T_BVH<:BVH, T_d<:Real}
    # If func is false, prune the whole subtree
    if !func(node)
        return Node{T_v, T_BVH, T_d}[]
    end

    # If the node passes the test and is a leaf, return it
    if node.is_leaf
        return Node{T_v, T_BVH, T_d}[node]
    end

    # Otherwise recurse on children
    nodes_L = leaf_search(node.left, func)
    nodes_R = leaf_search(node.right, func)

    return vcat(nodes_L, nodes_R)
end