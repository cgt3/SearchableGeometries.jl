```@meta
CurrentModule = SearchableGeometries
```

# [Ball](@id def_Ball)

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

## Related operations

-[`BoundingVolume`](@ref def_BoundingVolume)
-[`Containment`](@ref)
-[`Intersection`](@ref)