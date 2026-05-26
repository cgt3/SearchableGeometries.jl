using SearchableGeometries.GeometricPrimitives
using LinearAlgebra: norm

# Cone constructors ----------------------------------------------------------
@testset "Constructing cones: Invalid slope" begin
    @test_throws "SearchableGeometries.GeometricPrimitives.Cone: slope must be non-negative" Cone([0.0, 0.0], [1.0, 1.0], -1.0)
end

@testset "Constructing cones: Mismatched vertex and axis dimensions" begin
    @test_throws "SearchableGeometries.GeometricPrimitives.Cone: direction vector and source point dimension must match" Cone([0.0, 0.0], [1.0, 1.0, 1.0], 1.0)
end

@testset "Constructing cones: Zero axis vector" begin
    @test_throws "SearchableGeometries.GeometricPrimitives.Cone: Cannot construct cone with zero axis vector." Cone([0.0, 0.0], [0.0, 0.0], 1.0)
end

@testset "Constructing cones: Valid cone" begin
    cone = Cone([0.0, 0.0], [1.0, 1.0], 1.0)
    @test all(cone.vertex .== [0.0, 0.0])
    @test all(cone.axis .== [1.0, 1.0] ./ norm([1.0, 1.0]))
    @test cone.slope == 1.0
end

# `get_R` -----------------------------------------------------------
@testset "get_R: Point along axis" begin
    cone = Cone([0.0, 0.0], [1.0, 0.0], 1.0)
    p = [2.0, 0.0]
    @test get_R(cone, p) == 2.0
end

@testset "get_R: Point opposite axis direction" begin
    cone = Cone([0.0, 0.0], [1.0, 0.0], 1.0)
    p = [-2.0, 0.0]
    @test get_R(cone, p) == -2.0
end

@testset "get_R: Point off axis" begin
    cone = Cone([0.0, 0.0], [1.0, 0.0], 1.0)
    p = [2.0, 1.0]
    @test get_R(cone, p) == 2.0
end

# `get_radii` ----------------------------------------------------------------
@testset "get_radii: Point along axis" begin
    cone = Cone([0.0, 0.0], [1.0, 0.0], 1.0)
    p = [3.0, 4.0]
    R, r = get_radii(cone, p)
    @test R == 3.0
    @test r == 4.0
end

@testset "get_radii: Point opposite axis direction" begin
    cone = Cone([0.0, 0.0], [1.0, 0.0], 1.0)
    p = [-3.0, 4.0]
    R, r = get_radii(cone, p)
    @test R == -3.0
    @test r == 4.0
end

@testset "get_radii: Point off axis" begin
    cone = Cone([0.0, 0.0], [1.0, 0.0], 1.0)
    p = [2.0, 1.0]
    R, r = get_radii(cone, p)
    @test R == 2.0
    @test r == 1.0
end