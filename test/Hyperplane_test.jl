using SearchableGeometries
using LinearAlgebra: norm

# Hyperplane constructors --------------------------------------------------------
@testset "Hyperplane(point, n): Point and normal must have the same dimension" begin
    @test_throws "point and normal vector must have the same dimension" Hyperplane([0.0, 1.0], [1.0])
    @test_throws "point and normal vector must have the same dimension" Hyperplane([0.0], [1.0, 2.0])
end

@testset "Hyperplane(point, n): Normal vector must be nonzero" begin
    @test_throws "SearchableGeometries.Hyperplane: normal vector must be nonzero" Hyperplane([0.0, 0.0], [0.0, 0.0])
    @test_throws "SearchableGeometries.Hyperplane: normal vector must be nonzero" Hyperplane([1.0, 2.0, 3.0], [0.0, 0.0, 0.0])
end

@testset "Constructing Hyperplanes: Normal vector is normalized" begin
    plane1 = Hyperplane([0.0, 0.0], [1.0, 0.0])
    plane2 = Hyperplane([0.0, 0.0], [10.0, 0.0])

    @test plane1.n == [1.0, 0.0]
    @test plane2.n == [1.0, 0.0]
    @test all(plane1.active_dim .== plane2.active_dim)
    @test all(plane1.inactive_dim .== plane2.inactive_dim)
    @test plane1.dim == plane2.dim
end

@testset "Constructing Hyperplanes: Valid 2D constructor" begin
    plane = Hyperplane([1.0, 2.0], [3.0, 4.0])

    @test plane.point == [1.0, 2.0]
    @test plane.n == [3.0, 4.0] / norm([3.0, 4.0])
    @test plane.dim == 2
    @test plane.embedding_dim == 2
    @test all(plane.active_dim .== [1, 2])
    @test all(plane.inactive_dim .== [])
    @test all(plane.is_active .== [true, true])
end

@testset "Constructing Hyperplanes: Lower-dimensional Hyperplane" begin
    plane = Hyperplane([1.0, 2.0, 3.0, 4.0], [0.0, 5.0, 0.0, -12.0])

    @test plane.point == [1.0, 2.0, 3.0, 4.0]
    @test plane.n == [0.0, 5.0, 0.0, -12.0] / norm([0.0, 5.0, 0.0, -12.0])
    @test plane.dim == 2
    @test plane.embedding_dim == 4
    @test all(plane.active_dim .== [2, 4])
    @test all(plane.inactive_dim .== [1, 3])
    @test all(plane.is_active .== [false, true, false, true])
end

@testset "Constructing Hyperplanes: Constructor with axis-aligned normal" begin
    plane = Hyperplane([4.0, -1.0, 7.0], [0.0, 1.0, 0.0])

    @test plane.point == [4.0, -1.0, 7.0]
    @test plane.n == [0.0, 1.0, 0.0]
    @test plane.dim == 1
    @test plane.embedding_dim == 3
    @test all(plane.active_dim .== [2])
    @test all(plane.inactive_dim .== [1, 3])
    @test all(plane.is_active .== [false, true, false])
end

# Hyperplane -> BVs ------------------------------------------------------------

# `isContained(plane, query_pt)` ------------------------------------------------
@testset "isContained(plane, query_pt): Point dimension must match hyperplane embedding dimension" begin
    plane = Hyperplane([0.0, 0.0, 0.0], [1.0, 0.0, 0.0])

    @test_throws "point dimension(2) does not match hyperplane embedding dimension(3)" isContained(plane, [0.0, 0.0])
    @test_throws "point dimension(4) does not match hyperplane embedding dimension(3)" isContained(plane, [0.0, 0.0, 0.0, 0.0])
end

@testset "isContained(plane, query_pt): Point exactly on the hyperplane" begin
    plane = Hyperplane([0.0, 0.0], [1.0, 1.0])

    query_pt1 = [1.0, -1.0]
    @test isContained(plane, query_pt1)

    query_pt2 = [2.5, -2.5]
    @test isContained(plane, query_pt2)

    query_pt3 = [0.0, 0.0]
    @test isContained(plane, query_pt3)
end

@testset "isContained(plane, query_pt): Point not on the hyperplane" begin
    plane = Hyperplane([0.0, 0.0], [1.0, 1.0])

    query_pt1 = [1.0, 0.0]
    @test !isContained(plane, query_pt1)

    query_pt2 = [0.0, 1.0]
    @test !isContained(plane, query_pt2)

    query_pt3 = [2.0, -1.0]
    @test !isContained(plane, query_pt3)
end

@testset "isContained(plane, query_pt): Hyperplane with inactive dimensions" begin
    plane = Hyperplane([1.0, 2.0, 3.0], [0.0, 1.0, 0.0])

    query_pt1 = [10.0, 2.0, -7.0]
    @test isContained(plane, query_pt1)

    query_pt2 = [1.0, 2.0, 3.0]
    @test isContained(plane, query_pt2)

    query_pt3 = [10.0, 2.5, -7.0]
    @test !isContained(plane, query_pt3)
end
