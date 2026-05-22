using SafeTestsets

@safetestset "SearchableGeometries.jl" begin
    @safetestset "BoundingVolumes:" begin
        include("GeometricPrimitives/BoundingVolume_test.jl")
    end
    @safetestset "Balls:" begin
        include("GeometricPrimitives/Ball_test.jl")
    end
    @safetestset "Hyperplanes:" begin
        include("GeometricPrimitives/Hyperplane_test.jl")
    end
end