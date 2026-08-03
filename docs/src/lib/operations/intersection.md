# [Intersection](@id intersection)

```@meta
CurrentModule = SearchableGeometries.GeometricPrimitives
```

Intersection asks whether two geometric objects overlap and, when appropriate, returns a bounding representation of their overlap.

## Generic functions

```@docs
intersects
get_intersection
```

## Bounding volume with bounding volume

```@docs; canonical=false
intersects(::BoundingVolume, ::BoundingVolume)
get_intersection(::BoundingVolume, ::BoundingVolume)
```

## Bounding volume with ball

```@docs; canonical=false
intersects(::BoundingVolume, ::Ball)
get_intersection(::BoundingVolume, ::Ball)
```

## Bounding volume with hyperplane

```@docs; canonical=false
intersects(::BoundingVolume, ::Hyperplane)
get_intersection(::BoundingVolume, ::Hyperplane)
```

## Cone with bounding volume 

```@docs; canonical=false
intersects(::Cone, ::BoundingVolume)
get_intersection(::Cone, ::BoundingVolume)
```