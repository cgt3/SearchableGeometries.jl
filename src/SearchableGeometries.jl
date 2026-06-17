"""
    SearchableGeometries

Tools for representing and querying simple geometric objects used in search,
intersection, containment, and bounding-volume computations.

The package currently provides three main geometric objects:

- [`BoundingVolume`](@ref): an axis-aligned hyperrectangle described by lower
  and upper coordinate bounds.
- [`Ball`](@ref): a possibly lower-dimensional ``p``-norm ball embedded in a
  coordinate space.
- [`Hyperplane`](@ref): an affine hyperplane described by a point and a normal
  vector.

The public API is organized around a small set of generic operations:

- [`is_contained`](@ref): test whether a point or object lies inside another object.
- [`intersects`](@ref): test whether two objects intersect.
- [`get_intersection`](@ref): compute a bounding representation of an intersection.
- [`get_closest_point`](@ref): find a closest point in a bounding volume or on a hyperplane.
- [`get_furthest_point`](@ref): find a furthest point in a bounding volume.

# Boundary convention

Many predicates accept `include_boundary=true`. When `true`, points on the
boundary are treated as contained or intersecting. When `false`, strict interior
containment/intersection is required.

# Tolerance convention

Several methods accept `tol=DEFAULT_BV_POINT_TOL`. The tolerance is used to
handle floating-point roundoff, detect nearly inactive dimensions, and decide
whether values are close enough to geometric boundaries.
"""
module SearchableGeometries

using LinearAlgebra
using Plots
# using Printf

import Base.getindex

# Submodels
export GeometricPrimitives
include("GeometricPrimitives.jl")
using .GeometricPrimitives: BoundingVolume, Ball, Hyperplane, Cone, intersects, get_intersection, is_contained

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


const DEFAULT_PARTITION_PALETTE = Plots.palette([:darkred, :orangered3,  :darkorange2, :goldenrod1,
                                            :chartreuse3,   :forestgreen, :darkgreen,      #:olivedrab2
                                            :navy,          :blue,        :deepskyblue,    #:blue
                                            :mediumpurple1, :purple3,     :purple4  ])     #:purple1

# Data types
export IndexRange, Node, Tree

# Functions
export get_split, spatial_2_data, contains, get_containing_node, get_leaves, partition_plot, tree_plot

struct IndexRange
    first::Integer
    last::Integer
    n::Integer

    IndexRange(first; n) = new(first, first + n - 1, n)
    IndexRange(; first, last) = new(first, last, last - first + 1)
end

struct DataPoint{T}
    x::AbstractVector
    y::T
end

function spatial_2_data(data_num::Vector{Vector{T}}) where T <:Real
    return DataPoint.(data_num, nothing)
end

function Base.getindex(p::DataPoint{T}, i::Integer) where T
    return p.x[i]
end

function BoundingVolume(data::VDP, index_range::IndexRange) where {T, VDP<:Vector{DataPoint{T}}}
	ub = copy(data[index_range.first].x)
	lb = copy(data[index_range.first].x)

	for i = index_range.first:index_range.last
		d_ub = ub .< data[i].x
		ub[d_ub] = data[i].x[d_ub]

		d_lb = lv .> data[i].x
		lb[d_lb] = data[i].x[d_lb]
	end

	return BoundingVolume(lb, ub)
end

mutable struct Node{T, VDP<:Vector{DataPoint{T}}}
	data::VDP
	original_indices

end

function get_centroid(bv::BoundingVolume, data::VDP, index_range::IndexRange) where {T, VDP<:Vector{DataPoint{T}}}
	x0 = data[1].x
	centroid = zeros(eltype(x0), size(x0))

	for i in index_range.first:index_range.last
		centroid += data[i].x
	end

	centroid *= 1/index_range.n

	return centroid
end

function get_watertight_bounding_volume(parent::Union{Node, Nothing}, data::VDP, index_range::IndexRange) where {T, VDP<:Vector{DataPoint{T}}}
	if isnothing(parent)
		return BoundingVolume(data, index_range)
	else
		parent_bv = parent.bv_watertight
		if parent.index_range.first == index_range.first
			bv_ub = copy(parent_bv.ub)
			bv_ub[parent.split_dim] = parent.split_val
			return BoundingVolume(copy(parent_bv.lb), bv_ub)
		else
			bv_lb = copy(parent_bv.lb)
			bv_lb[parent.split_dim] = parent.split_value
			return BoundingVolume(bv_lb, copy(parent_bv.ub))
		end
	end
end

# function get_split(bv::BoundingVolume, bv_wt::BoundingVolume, data::VDP, index_range::IndexRange) where {T, VDP<:Vector{DataPoint{T}}}
# end

# function partition_data!(data::VDP, original_indices, index_range::IndexRange, split_dim, split_val) where {T, VDP<:Vector{DataPoint{T}}}
# end

# mutable struct Tree
# end

# function get_tree_stats(node::Union{Node, Nothing})
# end

# # User exposed function
# function get_containing_node(tree::Tree, query_pt; pt_tol=DEFAULT_NUM_LEAF_PTS)
# end

# function get_containing_node(node::Union{Node, Nothing}, query_pt; pt_tol=DEFAULT_PT_TOL)
# end

# function contains(tree::Tree, query_pr; pt_tol=DEFAULT_BV_POINT_TOL)
# end

# function contains(node::Union{Node, Nothing}, query_pt; pt_tol=DEFAULT_PT_TOL)
# end

# function tree_plot(tree::Tree;
#       x_spacing=DEFAULT_PLOT_WIDTH_PER_NODE,
#       y_spacing=DEFAULT_PLOT_HEIGHT_PER_NODE,
#       plot_text=true,
#       leaf_fontsize=DEFAULT_LEAF_FONTSIZE,
#       max_fontsize=DEFAULT_MAX_TREE_PLOT_FONTSIZE)
# end

# function tree_plot!(p, i_level::Integer, max_depth::Integer, node::Union{Node, Nothing}; 
#       plot_text=true,
#       leaf_fontsize=DEFAULT_TREE_PLOT_FONTSIZE,
#       max_fontsize=DEFAULT_MAX_TREE_PLOT_FONTSIZE )
# end

# function partition_plot(tree::Tree; 
#       size=(1500,1000), 
#       fontsize=20,
#       index=(1,2), 
#       watertight=true,
#       color_palette=DEFAULT_PARTITION_PALETTE,
#       linewidth=4)
# end

# function partition_plot!(p, node::Union{Node, Nothing}; 
#       index=(1,2),
#       watertight=true,
#       color_palette=DEFAULT_PARTITION_PALETTE,
#       linewidth=4)
# end

# function initialize_node_list(index_search::Bool)
# end

# function add_node!(nodes::Vector{Node{T, VDP}}, node::kdNode) where {T, VDP<:Vector{DataPoint{T}}}
# end

# function add_node!(nodes::Vector{IndexRange}, node::Node)
# end

# function get_leaves(tree::Tree; index_search=false::Bool)
# end

# function get_leaves!(leaves, node::Union{Node, Nothing}; index_search=false::Bool)
# end

end # module SearchableGeometries