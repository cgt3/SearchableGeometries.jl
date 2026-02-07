using SearchableGeometries
using SafeTestsets

@safetestset "SearchableGeometries.jl" begin
    @safetestset "BoundingVolumes:" begin include("BoundingVolume_test.jl") end
    @safetestset "Balls:" begin include("Ball_test.jl") end
end