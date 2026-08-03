# [Cones](@id cones_manual)

```@meta
CurrentModule = SearchableGeometries.GeometricPrimitives
```

A cone is a one-sided geometric region that expands outward from a vertex in the direction of an axis. 

A [`Cone`](@ref) is represented by:
- a `vertex`, which is the starting point of the cone;
- an `axis`, which determines the direction in which the cone opens;
- a nonnegative `slope`, which determines how quickly the cone widens.

The `Cone` constructor normalizes the axis vector automatically.

## Mathematical definition

Let ``v`` be the cone vertex, ``a`` its unit axis vector, and ``s`` its slope. For a point ``x``, define the axial distance 
```math 
R(x) = a^\mathsf{T}(x-v) 
``` 

and the perpendicular distance from the cone axis 

```math 
r(x) = \left\|x-v-R(x)a\right\|_2. 
``` 

A point is contained in the cone when 
```math 
R(x) \geq 0 
``` 

and 

```math 
r(x) \leq sR(x). 
``` 

Therefore, the cone is the set 
```math 
\left\{ x \in \mathbb{R}^n : R(x) \geq 0 \text{ and } r(x) \leq sR(x) \right\}. 
``` 

The condition ``R(x) \geq 0`` makes the cone one-sided. Points behind the vertex are not contained in the cone. When the cone has half-opening angle ``\theta``, its slope satisfies 

```math 
s = \tan(\theta). 
``` 

A cone with zero slope degenerates to its forward axis ray. 

## Constructing a cone 

A cone can be constructed from a vertex, an axis, and a slope. 

```julia 
using SearchableGeometries.GeometricPrimitives 

cone = Cone([0.0, 0.0], [1.0, 0.0], 1.0) 
``` 

In this example: 
- the vertex is the origin; 
- the cone opens in the positive ``x`` direction; 
- the radius increases by one unit for every unit traveled along the axis. 

Although `[1.0, 0.0]` is already a unit vector, the constructor accepts any nonzero axis vector and normalizes it. For example: 

```julia 
diagonal_cone = Cone([0.0, 0.0], [1.0, 1.0], 0.5) 
diagonal_cone.axis 
``` 
The normalized axis is approximately 

```julia 
[0.70710678, 0.70710678] 
``` 

The constructor throws an error when:
- the vertex and axis have different dimensions; 
- the axis is the zero vector; 
- the slope is negative. 

## Point containment 

Use [`is_contained`](@ref) to determine whether a point lies inside a cone. 

```julia 
cone = Cone([0.0, 0.0], [1.0, 0.0], 1.0) 
is_contained(cone, [2.0, 1.0]) 
``` 

This returns `true` because the point is in front of the vertex and its perpendicular distance from the cone axis is less than the cone radius at that axial position. 

A point outside the cone returns `false`. 
```julia 
is_contained(cone, [2.0, 3.0]) 
``` 

A point behind the cone vertex also returns `false`. 
```julia 
is_contained(cone, [-1.0, 0.0]) 
``` 

### Boundary behavior 

The `include_boundary` keyword controls whether points on the cone surface are included. 

```julia 
boundary_point = [2.0, 2.0] 
is_contained(cone, boundary_point; include_boundary=true) 
``` 

This returns `true`. 

```julia 
is_contained(cone, boundary_point; include_boundary=false) 
``` 

This returns `false`. 

The vertex is also considered part of the boundary. 

```julia 
is_contained(cone, cone.vertex; include_boundary=true) 
``` 

This returns `true`. 

```julia 
is_contained(cone, cone.vertex; include_boundary=false) 
``` 
This returns `false`. 

## Bounding-volume containment 

A cone may also contain an entire [`BoundingVolume`](@ref). 

```julia 
cone = Cone( [0.0, 0.0], [1.0, 0.0], 1.0) 
bv_inside = BoundingVolume( [1.0, -0.5], [2.0, 0.5]) 
is_contained(cone, bv_inside) 
``` 

This returns `true` because every point in the bounding volume lies inside the cone. If any part of the bounding volume lies outside the cone, the result is `false`. 

```julia 
bv_not_contained = BoundingVolume( [1.0, -2.0], [2.0, 2.0]) 
is_contained(cone, bv_not_contained)
``` 

## Intersecting a cone and a bounding volume 

Use [`intersects`](@ref) to determine whether a cone and a bounding volume have a nonempty intersection. 

```julia 
cone = Cone([0.0, 0.0], [1.0, 0.0], 1.0) 
bv = BoundingVolume([1.0, 0.5], [2.0, 1.5]) 
intersects(cone, bv) 
``` 

This returns `true`. 

### Boundary-only intersection 

Boundary-only contact can be included or excluded using `include_boundary`. 

```julia 
touching_bv = BoundingVolume([1.0, 2.0], [2.0, 3.0]) 
intersects(cone, touching_bv; include_boundary=true) 
``` 

This returns `true` when the bounding volume touches the cone boundary. 

```julia 
intersects(cone, touching_bv; include_boundary=false) 
``` 

This returns `false` when the only contact occurs on the boundary. 

## Bounding the intersection 

Use [`get_intersection`](@ref) to obtain an axis-aligned bounding volume that encloses the intersection between a cone and a bounding volume. 

```julia 
intersection_bv = get_intersection(cone, bv) 
``` 

The returned object is a conservative bounding representation. The exact intersection between a cone and a bounding volume may have curved or sloped boundaries. Therefore, the exact intersection is not generally another axis-aligned bounding volume. Instead, `get_intersection` returns a bounding volume that contains the exact intersection. 

If the cone and bounding volume do not intersect, the method returns an empty bounding volume. 
```julia 
outside_bv = BoundingVolume([1.0, 3.0], [2.0, 4.0]) 
intersection_bv = get_intersection(cone, outside_bv) 
intersection_bv.is_empty 
``` 

This returns `true`. 

## Bounding a truncated cone 
A finite section of a cone can be enclosed by a bounding volume using its minimum and maximum axial radii. 

```julia 
cone = Cone([0.0, 0.0], [1.0, 0.0], 1.0) 
cone_bv = BoundingVolume(cone, 0.0, 2.0) 
``` 

The resulting bounds are 
```julia 
cone_bv.lb 
``` 

which gives 

```julia 
[0.0, -2.0] 
``` 

and 

```julia 
cone_bv.ub 
``` 

which gives 

```julia 
[2.0, 2.0] 
``` 

The radii must satisfy 

```math 
0 \leq R_{\min} \leq R_{\max}. 
``` 

An error is thrown when: 
- `R_min` is negative; 
- `R_max` is smaller than `R_min`. 

## See also 
- [`Cone`](@ref) 
- [`BoundingVolume`](@ref) 
- [`is_contained`](@ref) 
- [`intersects`](@ref) 
- [`get_intersection`](@ref)