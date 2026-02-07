using SearchableGeometries
using SafeTestsets

@safetestset "SearchableGeometries.jl" begin
    @safetestset "BoundingVolumes:" begin include("BoundingVolume_test.jl") end
end
