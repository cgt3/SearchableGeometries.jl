# [Containment](@id containment)

```@meta
CurrentModule = SearchableGeometries.GeometricPrimitives
```

Containment asks whether one object lies inside another object.

The generic function is [`is_contained`](@ref).

## Generic function

```@docs
is_contained
```

## Bounding volumes

```@docs; canonical=false
is_contained(::BoundingVolume, ::Vector{<:Real})
is_contained(::BoundingVolume, ::BoundingVolume)
```

## Balls

```@docs; canonical=false
is_contained(::Ball, ::Vector{<:Real})
is_contained(::BoundingVolume, ::Ball)
is_contained(::Ball, ::BoundingVolume)
```

## Hyperplanes

```@docs; canonical=false
is_contained(::Hyperplane, ::Vector{<:Real})
```