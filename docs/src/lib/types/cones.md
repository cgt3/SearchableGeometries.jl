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

``` The internal helper functions should not be placed on the public API page unless they are intended to become part of the public interface. These internal helpers include functions such as: 

```julia 
get_R get_radii 
get_alpha 
get_bound_lines 
get_bounding_radii 
```