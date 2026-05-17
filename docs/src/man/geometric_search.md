# [Geometric Search](@id geometric_search)

```@meta
CurrentModule = SearchableGeometries
```

Geometric search is the process of finding points, regions, or objects that satisfy a spatial relationship.

Examples include:

- finding points close to a query point;
- checking whether a point lies inside a region;
- deciding whether two regions intersect;
- pruning regions that cannot contain a search result;
- splitting space into smaller searchable pieces.

SearchableGeometries.jl is designed to provide the geometric tools needed for these kinds of search problems.

## Why geometry matters in search

Many search algorithms become efficient because they avoid checking every object directly.

Instead of comparing a query point with every point in a dataset, a search algorithm can group points into geometric regions and ask simpler questions first:

- Does this region contain the query point?
- Is this region close enough to matter?
- Does this region intersect the query region?
- Can this region be safely ignored?

If the answer shows that a region cannot contain the desired result, the algorithm can skip that region entirely. This is called pruning.

## Example: bounding-volume pruning

Suppose a set of points is enclosed by a bounding volume.

```julia
using SearchableGeometries

bv = BoundingVolume([0.0, 0.0], [2.0, 1.0])
query_pt = [5.0, 0.5]

closest_pt = get_closest_point(bv, query_pt)
```

The closest point in the bounding volume to the query point gives a lower bound on how close any point inside the bounding volume can be to the query point.

If this lower bound is already too large, then a search algorithm may be able to ignore the entire bounding volume.

This is one reason bounding volumes are useful in spatial search.

## Example: ball queries

A ball describes all points within a fixed distance of a center point.

```julia
ball = Ball([0.0, 0.0], 1.0; p=2)
```

This can represent a query region such as:

```math
\{x : \|x - c\|_2 \le 1\}.
```

A search algorithm can then ask whether another region intersects this ball:

```julia
bv = BoundingVolume([0.5, -0.5], [2.0, 0.5])

intersects(bv, ball)
```

If the bounding volume does not intersect the ball, then none of the points inside that bounding volume can be part of the ball query.

## Example: hyperplane splitting

A hyperplane can split space into two sides.

```julia
plane = Hyperplane([0.0, 0.0], [1.0, -1.0])
```

In two dimensions, this represents the line

```math
x - y = 0.
```

Hyperplanes are important in search because many spatial data structures split space recursively. For example, a kd-tree uses coordinate-aligned splitting hyperplanes to divide points into smaller regions.

SearchableGeometries.jl provides tools for checking how hyperplanes interact with bounding volumes:

```julia
bv = BoundingVolume([0.0, 0.0], [2.0, 2.0])

intersects(bv, plane)
get_intersection(bv, plane)
```

These operations are useful when a search structure needs to know whether a region crosses a splitting boundary.

## Current scope of the package

At the moment, SearchableGeometries.jl focuses on the geometric primitives needed for search:

- [Bounding Volumes](@ref "bounding_volumes_manual")
- [Balls](@ref "balls_manual")
- [Hyperplanes](@ref "hyperplanes_manual")

These objects support operations such as containment, intersection, closest-point queries, furthest-point queries, and bounding-volume tightening.

## Future direction

The long-term goal is to build on these geometric primitives to support searchable structures such as:

- kd-trees;
- ball trees;
- bounding-volume hierarchies;
- graph-based search structures;
- nearest-neighbor and range-search algorithms.

The geometry layer comes first because these search structures need reliable geometric operations before they can be implemented efficiently.

## See also

- [Getting Started](@ref "getting_started")
- [Bounding Volumes](@ref "bounding_volumes_manual")
- [Balls](@ref "balls_manual")
- [Hyperplanes](@ref "hyperplanes_manual")