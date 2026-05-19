# [Ball](@id def_Ball)

```@meta
CurrentModule = SearchableGeometries
```

```@docs
Ball
```

## Constructors

```julia
Ball(center, radius; p=2, active_indices=true, indices=eachindex(center))
```

## Conversion

```@docs
BoundingVolume(::Ball)
```

## Operations

```@docs
is_contained(::Ball, ::Vector{<:Real})
is_contained(::BoundingVolume, ::Ball)
is_contained(::Ball, ::BoundingVolume)
intersects(::BoundingVolume, ::Ball)
get_intersection(::BoundingVolume, ::Ball)
get_reduced_dim_ball
tighten_bv_bounds!
```