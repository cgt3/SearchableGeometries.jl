# [Cone](@id def_Cone) 

```@meta 
CurrentModule = SearchableGeometries.GeometricPrimitives 
```

```@docs 
Cone 
```

## Bounding-volume construction 

```@docs 
BoundingVolume(::Cone, ::Real, ::Real) 
```

## Point operations 
```@docs 
is_contained(::Cone, ::Vector{<:Real}) 
```

## Bounding-volume operations 
```@docs 
is_contained(::Cone, ::BoundingVolume) 
intersects(::Cone, ::BoundingVolume) 
get_intersection(::Cone, ::BoundingVolume) 
```