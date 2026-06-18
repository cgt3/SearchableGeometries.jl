import .GeometricPrimitives: BoundingVolume

const DEFAULT_NUM_LEAF_PTS = 40
const DEFAULT_PT_TOL = 1e-12

const DEFAULT_PLOT_WIDTH_PER_NODE = 100
const DEFAULT_PLOT_HEIGHT_PER_NODE = 200
const DEFAULT_LEAF_FONTSIZE = 8
const DEFAULT_MAX_TREE_PLOT_FONTSIZE = 24

const MIN_POINT_PLOT_SIZE = 2
const MAX_POINT_PLOT_SIZE = 16
const MIN_SPLIT_LINE_WIDTH = 3
const MAX_SPLIT_LINE_WIDTH = 10
const TREE_MARKER_SIZE = 20


const DEFAULT_PARTITION_PALETTE = Plots.palette(
    [:darkred, :orangered3, :darkorange2, :goldenrod1,
    :chartreuse3, :forestgreen, :darkgreen, #:olivedrab2
    :navy, :blue, :deepskyblue, #:blue
    :mediumpurple1, :purple3, :purple4]
)

# Constants
export DEFAULT_NUM_LEAF_PTS

# Data types
export PointRange, Node, Watertight, Arbitrary

# Functions
export get_split

struct PointRange{T_data, dim}
    data::Array{T_data, dim}
    first::Int
    last::Int
    n::Int

    function PointRange(data::Array{T_data, dim}, first::Integer; n::Integer) where {T_data, dim}
        n < 0 && throw("SearchableGeometries.PointRange: n must be nonnegative")

        first_int = Int(first)
        n_int = Int(n)

        return new{T_data, dim}(data, first_int, first_int + n_int - 1, n_int)
    end

    function PointRange(data::Array{T_data, dim}, first::Integer, last::Integer) where {T_data, dim}
        first_int = Int(first)
        last_int = Int(last)
        n_int = last_int - first_int + 1

        n_int < 0 && throw("SearchableGeometries.PointRange: last must be at least first - 1")

        return new{T_data, dim}(data, first_int, last_int, n_int)
    end
end

Base.length(r::PointRange) = r.n
Base.isempty(r::PointRange) = r.n == 0

abstract type BVH end

struct Watertight <: BVH end
struct Arbitrary <: BVH end

mutable struct Node{T_v, T_BVH<:BVH, T_x}
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
        is_base_case::Function,
        get_split::Function,
        partition_data!::Function,
        partition_bv::Function,
        parent::Union{Node{T_v, T_BVH, T_x}, Nothing}=nothing
    ) where {T_v, T_BVH<:BVH, T_x}
        # Base case: return nothing
        if is_base_case(val, bv, parent)
            return new{T_v, T_BVH, T_x}(
                parent,
                nothing,
                nothing,
                0,
                zero(T_x),
                true,
                val,
                bv,
            )
        end
        
        # Create current internal node
        node = new{T_v, T_BVH, T_x}(
            parent,
            nothing,
            nothing,
            0,
            zero(T_x),
            false,
            val,
            bv
        )

        # Compute split
        split_dim, split_val = get_split(node)
        node.split_dim = split_dim
        node.split_val = convert(T_x, split_val)
        
        # Partition data
        val_L, val_R = partition_data!(node)

        # Important safety check:
        # if either child is empty, this split is not useful.
        # Stop here and make the current node a leaf.
        if isempty(val_L) || isempty(val_R)
            node.is_leaf = true
            node.left = nothing
            node.right = nothing
            return node
        end

        # Partition bounding volumes
        bv_L, bv_R = partition_bv(node, val_L, val_R)

        # Recursively build children
        node.left = Node{T_v, T_BVH, T_x}(
            val_L,
            bv_L;
            is_base_case=is_base_case,
            get_split = get_split,
            partition_data! = partition_data!,
            partition_bv = partition_bv,
            parent = node
        )
        node.right = Node{T_v, T_BVH, T_x}(
            val_R,
            bv_R;
            is_base_case = is_base_case,
            get_split = get_split,
            partition_data! = partition_data!,
            partition_bv = partition_bv,
            parent = node
        )

        return node
    end
end

function Node(
    val::T_v,
    bv::BoundingVolume,
    ::Type{T_BVH};
    kwargs...,
) where {T_v, T_BVH<:BVH}

    T_x = float(promote_type(eltype(bv.lb), eltype(bv.ub)))

    return Node{T_v, T_BVH, T_x}(
        val,
        bv;
        kwargs...,
    )
end

function is_base_case(val::PointRange, bv::BoundingVolume, parent; leafsize::Int=DEFAULT_NUM_LEAF_PTS)
    # Stop if this node is already small enough
    if val.n <= leafsize
        return true
    end

    # Root has no parent
    if isnothing(parent)
        return false
    end

    # Stop if the child has the same range as the parent.
    # This prevents infinite recursion.
    if val.first == parent.val.first && val.last == parent.val.last
        return true
    end

    return false
end

function get_split(node)
    val = node.val
    data = val.data

    val.n == 0 && throw("SearchableGeometries.get_split: cannot split an empty PointRange.")

    # Compute centroid of points in this node.
    x0 = data[val.first].x
    centroid = zeros(float(eltype(x0)), size(x0))

    @inbounds for i in val.first:val.last
        centroid .+= data[i].x
    end

    centroid ./= val.n

    # Choose the dimension that splits the BV into the most square pieces.
    bv_length = node.bv.ub .- node.bv.lb
    split_dim = argmax(bv_length)
    split_val = centroid[split_dim]

    return split_dim, split_val
end

function get_split(node::Node{<:PointRange{T, 2}}) where {T}
    val = node.val
    data = val.data

    val.n == 0 && throw("SearchableGeometries.get_split: cannot split an empty PointRange.")

    n_dims = size(data, 1)
    centroid = zeros(float(T), n_dims)

    @inbounds for i in val.first:val.last
        centroid .+= @view data[:, i]
    end

    centroid ./= val.n

    bv_length = node.bv.ub .- node.bv.lb
    split_dim = argmax(bv_length)
    split_val = centroid[split_dim]

    return split_dim, split_val
end

function partition_data!(node)
    val = node.val
    data = val.data

    val.n <= 1 && throw("SearchableGeometries.partition_data!: cannot partition a PointRange with $(val.n) point(s).")

    split_dim = node.split_dim
    split_val = node.split_val

    i_start = val.first
    i_end = val.last

    while i_start < i_end
        is_valid_L = data[i_start].x[split_dim] <= split_val
        is_valid_R = data[i_end].x[split_dim] > split_val

        if !is_valid_L && !is_valid_R
            data[i_start], data[i_end] = data[i_end], data[i_start]
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

    if data[last_left].x[split_dim] > split_val
        last_left -= 1
    end

    val_L = PointRange(data, val.first, last_left)
    val_R = PointRange(data, last_left + 1, val.last)

    return val_L, val_R
end

function partition_data!(node::Node{<:PointRange{T, 2}}) where {T}
    val = node.val
    data = val.data

    val.n <= 1 && throw("SearchableGeometries.partition_data!: cannot partition a PointRange with $(val.n) point(s).")

    split_dim = node.split_dim
    split_val = node.split_val

    i_start = val.first
    i_end = val.last

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

    val_L = PointRange(data, val.first, last_left)
    val_R = PointRange(data, last_left + 1, val.last)

    return val_L, val_R
end

function BoundingVolume(point_range::PointRange)
    isempty(point_range) && throw("SearchableGeometries.BoundingVolume: cannot build a BoundingVolume from an empty PointRange")
    
    data = point_range.data
	ub = copy(data[point_range.first].x)
	lb = copy(data[point_range.first].x)

	for i = point_range.first:point_range.last
		d_ub = ub .< data[i].x
		ub[d_ub] = data[i].x[d_ub]

		d_lb = lb .> data[i].x
		lb[d_lb] = data[i].x[d_lb]
	end

	return BoundingVolume(lb, ub)
end

function BoundingVolume(point_range::PointRange{T, 2}) where {T}
    isempty(point_range) && throw("SearchableGeometries.BoundingVolume: cannot build a BoundingVolume from an empty PointRange.")

    data = point_range.data

    lb = copy(@view data[:, point_range.first])
    ub = copy(@view data[:, point_range.first])

    @inbounds for i in point_range.first:point_range.last
        x = @view data[:, i]

        d_ub = ub .< x
        ub[d_ub] .= x[d_ub]

        d_lb = lb .> x
        lb[d_lb] .= x[d_lb]
    end

    return BoundingVolume(lb, ub)
end

function partition_bv(
    node::Node{T_v, Arbitrary, T_x},
    val_L::PointRange,
    val_R::PointRange
) where {T_v, T_x}

    return BoundingVolume(val_L), BoundingVolume(val_R)
end

function partition_bv(
    node::Node{T_v, Watertight, T_x},
    val_L::PointRange,
    val_R::PointRange
) where {T_v, T_x}
    lb_L = copy(node.bv.lb)
    ub_L = copy(node.bv.ub)

    lb_R = copy(node.bv.lb)
    ub_R = copy(node.bv.ub)

    ub_L[node.split_dim] = node.split_val
    lb_R[node.split_dim] = node.split_val

    return BoundingVolume(lb_L, ub_L), BoundingVolume(lb_R, ub_R)
end