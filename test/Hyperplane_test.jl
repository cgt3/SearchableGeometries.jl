using SearchableGeometries
using LinearAlgebra: norm, dot

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

# Hyperplane -> Points ----------------------------------------------------------

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

# Hyperplane -> BVs --------------------------------------------------------------
# `intersects(bv, plane)` --------------------------------------------------------
@testset "intersects(bv, plane): Dimension mismatch throws" begin
    bv = BoundingVolume([0.0, 0.0], [1.0, 1.0])
    plane = Hyperplane([0.0, 0.0, 0.0], [1.0, 0.0, 0.0])

    @test_throws "SearchableGeometries.Hyperplane: bounding volume dimension(2) does not match hyperplane embedding dimension(3)" intersects(bv, plane)
end

@testset "intersects(bv, plane): Plane crosses the interior of a full-dimensional BV" begin
    bv = BoundingVolume([0.0, 0.0], [1.0, 1.0])
    plane = Hyperplane([0.0, 1.0], [1.0, 1.0])  # x + y = 1

    @test intersects(bv, plane; include_boundary=true)
    @test intersects(bv, plane; include_boundary=false)
end

@testset "intersects(bv, plane): Plane touching only a boundary face" begin
    bv = BoundingVolume([0.0, 0.0], [1.0, 1.0])
    plane = Hyperplane([0.0, 0.0], [1.0, 0.0])  # x = 0

    @test intersects(bv, plane; include_boundary=true)
    @test !intersects(bv, plane; include_boundary=false)
end

@testset "intersects(bv, plane): Plane touching only at a corner" begin
    bv = BoundingVolume([0.0, 0.0], [1.0, 1.0])
    plane = Hyperplane([0.0, 0.0], [1.0, 1.0])  # x + y = 0

    @test intersects(bv, plane; include_boundary=true)
    @test !intersects(bv, plane; include_boundary=false)
end

@testset "intersects(bv, plane): Plane completely outside the BV" begin
    bv = BoundingVolume([0.0, 0.0], [1.0, 1.0])
    plane = Hyperplane([2.0, 0.0], [1.0, 0.0])  # x = 2

    @test !intersects(bv, plane; include_boundary=true)
    @test !intersects(bv, plane; include_boundary=false)
end

@testset "intersects(bv, plane): Lower-dimensional BV lying on the plane" begin
    bv = BoundingVolume([0.0, -1.0], [0.0, 1.0])  # segment x = 0
    plane = Hyperplane([0.0, 0.0], [1.0, 0.0])    # plane x = 0

    @test intersects(bv, plane; include_boundary=true)
    @test !intersects(bv, plane; include_boundary=false)
end

# # `getClosestPoint(bv, plane)` ------------------------------------------------
# @testset "getClosestPoint(bv, plane): Dimension mismatch throws" begin
#     bv = BoundingVolume([0.0, 0.0], [1.0, 1.0])
#     plane = Hyperplane([0.0, 0.0, 0.0], [1.0, 0.0, 0.0])

#     @test_throws "SearchableGeometries.Hyperplane: bounding volume dimension(2) does not match hyperplane embedding dimension(3)" getClosestPoint(bv, plane)
# end

# @testset "getClosestPoint(bv, plane): BV completely on the positive side of the plane" begin
#     bv = BoundingVolume([2.0, -1.0], [3.0, 4.0])
#     plane = Hyperplane([0.0, 0.0], [1.0, 0.0])   # plane x = 0

#     closest_pt = getClosestPoint(bv, plane)

#     # The whole BV is to the right of x = 0, so the closest face is x = 2.
#     # Among all points on that face, we choose the lexicographically smallest one.
#     @test closest_pt == [2.0, -1.0]
#     @test isapprox(abs(dot(plane.n, closest_pt - plane.point)), 2.0; atol=1e-12)
# end

# @testset "getClosestPoint(bv, plane): BV completely on the negative side of the plane" begin
#     bv = BoundingVolume([0.0, -1.0], [1.0, 4.0])
#     plane = Hyperplane([2.0, 0.0], [1.0, 0.0])   # plane x = 2

#     closest_pt = getClosestPoint(bv, plane)

#     # The whole BV is to the left of x = 2, so the closest face is x = 1.
#     # Among all points on that face, we choose the lexicographically smallest one.
#     @test closest_pt == [1.0, -1.0]
#     @test isapprox(abs(dot(plane.n, closest_pt - plane.point)), 1.0; atol=1e-12)
# end

# @testset "getClosestPoint(bv, plane): Lower-dimensional BV lying on the plane" begin
#     bv = BoundingVolume([0.0, -1.0], [0.0, 1.0])   # segment x = 0
#     plane = Hyperplane([0.0, 0.0], [1.0, 0.0])     # plane x = 0

#     closest_pt = getClosestPoint(bv, plane)

#     # Every point of the segment lies on the plane.
#     # The lexicographically smallest point is (0,-1).
#     @test closest_pt == [0.0, -1.0]
#     @test isContained(plane, closest_pt)
# end

# @testset "getClosestPoint(bv, plane): Plane with inactive dimensions in the normal" begin
#     bv = BoundingVolume([0.0, 0.0, 0.0], [3.0, 4.0, 5.0])
#     plane = Hyperplane([0.0, 2.0, 0.0], [0.0, 1.0, 0.0])   # plane y = 2

#     closest_pt = getClosestPoint(bv, plane)

#     # The intersection is all points with y = 2 inside the box.
#     # The lexicographically smallest such point is (0,2,0).
#     @test closest_pt == [0.0, 2.0, 0.0]
#     @test isContained(plane, closest_pt)
# end

# @testset "getClosestPoint(bv, plane): Boundary-only intersection" begin
#     bv = BoundingVolume([0.0, 0.0], [1.0, 1.0])
#     plane = Hyperplane([0.0, 0.0], [1.0, 0.0])   # plane x = 0

#     closest_pt = getClosestPoint(bv, plane)

#     # The plane meets the BV on the left face x = 0.
#     # Lexicographically smallest point on that face is (0,0).
#     @test closest_pt == [0.0, 0.0]
#     @test isContained(plane, closest_pt)
# end

# @testset "getClosestPoint(bv, plane): Corner-only intersection" begin
#     bv = BoundingVolume([0.0, 0.0], [1.0, 1.0])
#     plane = Hyperplane([0.0, 0.0], [1.0, 1.0])   # plane x + y = 0

#     closest_pt = getClosestPoint(bv, plane)

#     # The only intersection point is the corner (0,0).
#     @test closest_pt == [0.0, 0.0]
#     @test isContained(plane, closest_pt)
# end
