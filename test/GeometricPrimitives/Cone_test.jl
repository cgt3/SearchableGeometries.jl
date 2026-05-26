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
@testset "get_R(cone, pt): Point along axis" begin
    cone = Cone([0.0, 0.0], [1.0, 0.0], 1.0)
    p = [2.0, 0.0]
    @test get_R(cone, p) == 2.0
end

@testset "get_R(cone, pt): Point opposite axis direction" begin
    cone = Cone([0.0, 0.0], [1.0, 0.0], 1.0)
    p = [-2.0, 0.0]
    @test get_R(cone, p) == -2.0
end

@testset "get_R(cone, pt): Point off axis" begin
    cone = Cone([0.0, 0.0], [1.0, 0.0], 1.0)
    p = [2.0, 1.0]
    @test get_R(cone, p) == 2.0
end

# `get_radii` ----------------------------------------------------------------
@testset "get_radii(cone, pt): Point along axis" begin
    cone = Cone([0.0, 0.0], [1.0, 0.0], 1.0)
    p = [3.0, 4.0]
    R, r = get_radii(cone, p)
    @test R == 3.0
    @test r == 4.0
end

@testset "get_radii(cone, pt): Point opposite axis direction" begin
    cone = Cone([0.0, 0.0], [1.0, 0.0], 1.0)
    p = [-3.0, 4.0]
    R, r = get_radii(cone, p)
    @test R == -3.0
    @test r == 4.0
end

@testset "get_radii(cone, pt): Point off axis" begin
    cone = Cone([0.0, 0.0], [1.0, 0.0], 1.0)
    p = [2.0, 1.0]
    R, r = get_radii(cone, p)
    @test R == 2.0
    @test r == 1.0
end

# `is_contained(cone, point)` ----------------------------------------------------------------
@testset "is_contained(cone, pt): Point inside cone" begin
    cone = Cone([0.0, 0.0], [1.0, 0.0], 1.0)
    p = [2.0, 1.0]
    @test is_contained(cone, p; include_boundary=true)
    @test is_contained(cone, p; include_boundary=false)
end

@testset "is_contained(cone, pt): Point on boundary" begin
    cone = Cone([0.0, 0.0], [1.0, 0.0], 1.0)
    p = [2.0, 2.0]
    @test is_contained(cone, p; include_boundary=true)
    @test !is_contained(cone, p; include_boundary=false)
end

@testset "is_contained(cone, pt): Point outside cone" begin
    cone = Cone([0.0, 0.0], [1.0, 0.0], 1.0)
    p = [2.0, 3.0]
    @test !is_contained(cone, p; include_boundary=true)
    @test !is_contained(cone, p; include_boundary=false)
end

@testset "is_contained(cone, pt): Point behind vertex" begin
    cone = Cone([0.0, 0.0], [1.0, 0.0], 1.0)
    p = [-1.0, 0.5]
    @test !is_contained(cone, p; include_boundary=true)
    @test !is_contained(cone, p; include_boundary=false)
end

@testset "is_contained(cone, pt): Point at vertex" begin
    cone = Cone([0.0, 0.0], [1.0, 0.0], 1.0)
    p = [0.0, 0.0]
    @test is_contained(cone, p; include_boundary=true)
    @test !is_contained(cone, p; include_boundary=false)
end

@testset "is_contained(cone, pt): Point on axis" begin
    cone = Cone([0.0, 0.0], [1.0, 0.0], 1.0)
    p = [2.0, 0.0]
    @test is_contained(cone, p; include_boundary=true)
    @test is_contained(cone, p; include_boundary=false)
end

@testset "is_contained(cone, pt): Zero-slope cone only contains points on the axis" begin
    cone = Cone([0.0, 0.0], [1.0, 0.0], 0.0)
    p_on_axis = [3.0, 0.0]
    p_off_axis = [3.0, 0.1]

    @test is_contained(cone, p_on_axis; include_boundary=true)
    @test !is_contained(cone, p_on_axis; include_boundary=false)

    @test !is_contained(cone, p_off_axis; include_boundary=true)
    @test !is_contained(cone, p_off_axis; include_boundary=false)
end

# `get_bounding_radii(cone, bv)` ----------------------------------------------------------------
@testset "get_bounding_radii(cone, bv): Cone and bounding volume dimensions do not match" begin
    cone = Cone([0.0, 0.0], [1.0, 0.0], 1.0)
    bv = BoundingVolume([-1.0, -1.0, -1.0], [1.0, 1.0, 1.0])
    @test_throws "SearchableGeometries.GeometricPrimitives.get_antiextreme_point: cone dimension(2) does not match bounding volume dimension(3)" get_bounding_radii(cone, bv)
end

@testset "get_bounding_radii(cone, bv): 3D cone aligned with z-axis and centered at origin" begin
    cone = Cone([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 1.0)
    bv = BoundingVolume([-1.0, -1.0, -1.0], [1.0, 1.0, 2.0])
    
    R_min, R_max = get_bounding_radii(cone, bv)
    @test R_min == -1.0
    @test R_max == 2.0
end

@testset "get_bounding_radii(cone, bv): 3D cone with arbitrary axis and vertex" begin
    cone = Cone([1.0, 2.0, 3.0], [1.0, 1.0, 1.0], 1.0)
    bv = BoundingVolume([0.0, 0.0, 0.0], [4.0, 4.0, 4.0])
    
    R_min, R_max = get_bounding_radii(cone, bv)
    @test isapprox(R_min, -2sqrt(3), atol=1e-6)
    @test isapprox(R_max, 2sqrt(3), atol=1e-6)
end

@testset "get_bounding_radii(cone, bv): 2D cone with zero slope" begin
    cone = Cone([0.0, 0.0], [1.0, 0.0], 0.0)
    bv = BoundingVolume([-1.0, -1.0], [1.0, 1.0])
    
    R_min, R_max = get_bounding_radii(cone, bv)
    @test R_min == -1.0
    @test R_max == 1.0
end