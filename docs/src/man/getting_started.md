# [Getting Started](@id getting_started)

```@meta
CurrentModule = SearchableGeometries
```

This page gives a quick introduction to using `SearchableGeometries.jl`.

The package provides geometric objects and query functions that are useful for building searchable geometric structures. At the moment, the main focus is on the geometry layer: bounding volumes, balls, and hyperplanes.

## Installation

Since `SearchableGeometries.jl` is a registered Julia package, install it from the Julia package manager with:

```julia
using Pkg

Pkg.add("SearchableGeometries")
```

Then load the package with:

```julia
using SearchableGeometries
```

If you are developing the package locally, activate the package environment from the repository root:

```julia
using Pkg

Pkg.activate(".")
Pkg.instantiate()
```

## Main geometric objects

SearchableGeometries.jl currently provides three main geometry types:

- `BoundingVolume`: an axis-aligned box or hyperrectangle;
- `Ball`: a ``p``-norm ball centered at a point;
- `Hyperplane`: a flat affine surface represented by a point and a normal vector.

These objects support common geometric queries such as containment, intersection, closest-point queries, and furthest-point queries.

## Bounding volumes

A bounding volume is an axis-aligned region defined by lower and upper coordinate bounds.

```julia
using SearchableGeometries

bv = BoundingVolume([0.0, 0.0], [2.0, 1.0])
```

This creates the rectangle

```math
[0,2] \times [0,1].
```

You can test whether a point lies inside the bounding volume:

```julia
is_contained(bv, [1.0, 0.5])   # true
is_contained(bv, [3.0, 0.5])   # false
```

You can also find the closest point in the bounding volume to a query point:

```julia
get_closest_point(bv, [3.0, 0.5])
```

The closest point is `[2.0, 0.5]`, because the query point lies to the right of the bounding volume.

Similarly, you can find a furthest point:

```julia
get_furthest_point(bv, [3.0, 0.5])
```

See [Bounding Volumes](@ref bounding_volumes_manual) for more details.

## Balls

A ball represents all points within a fixed distance of a center point.

```julia
ball = Ball([0.0, 0.0], 1.0; p=2)
```

This creates the Euclidean unit ball centered at the origin:

```math
\{x \in \mathbb{R}^2 : \|x\|_2 \le 1\}.
```

You can test point containment:

```julia
is_contained(ball, [0.5, 0.5])   # true
is_contained(ball, [2.0, 0.0])   # false
```

You can also construct a bounding volume that encloses the ball:

```julia
bv_ball = BoundingVolume(ball)
```

This is useful because bounding volumes are cheap to query and can be used as fast approximations of more complicated geometric objects.

Balls can also be checked against bounding volumes:

```julia
bv = BoundingVolume([0.5, -0.5], [2.0, 0.5])

intersects(bv, ball)
get_intersection(bv, ball)
```

See [Balls](@ref balls_manual) for more details.

## Hyperplanes

A hyperplane is a flat affine surface.

In two dimensions, a hyperplane is a line. In three dimensions, it is a plane. In higher dimensions, it is a codimension-one affine subspace.

For example, the line

```math
x - y = 0
```

can be represented as:

```julia
plane = Hyperplane([0.0, 0.0], [1.0, -1.0])
```

The first vector is a point on the hyperplane. The second vector is a normal vector.

You can test whether a point lies on the hyperplane:

```julia
is_contained(plane, [2.0, 2.0])   # true
is_contained(plane, [2.0, 0.0])   # false
```

You can also project a point onto the hyperplane:

```julia
get_closest_point([2.0, 0.0], plane)
```

For this example, the closest point is `[1.0, 1.0]`.

Hyperplanes can also be checked against bounding volumes:

```julia
bv = BoundingVolume([0.0, 0.0], [2.0, 2.0])

intersects(bv, plane)
get_intersection(bv, plane)
```

See [Hyperplanes](@ref hyperplanes_manual) for more details.

## Common workflow

A typical workflow is:

1. Construct a geometric object.
2. Query containment or intersection.
3. Use closest-point, furthest-point, or intersection routines to support a search algorithm.

For example:

```julia
using SearchableGeometries

bv = BoundingVolume([0.0, 0.0], [2.0, 1.0])
query_pt = [3.0, 0.5]

if is_contained(bv, query_pt)
    println("The query point is inside the bounding volume.")
else
    closest_pt = get_closest_point(bv, query_pt)
    println("The query point is outside. The closest point is ", closest_pt)
end
```

## Next steps

After this page, a good order is:

- [Geometric Search](@ref geometric_search)
- [Bounding Volumes](@ref bounding_volumes_manual)
- [Balls](@ref balls_manual)
- [Hyperplanes](@ref hyperplanes_manual)