using SearchableGeometries

# Constructing BVs: ------------------------------------------------------
@testset "Constructing BVs: Invalid bounds" begin
    lb = [1, 2]
    ub = [-1, -2]
    @test_throws "Cannot construct bounding volume" BoundingVolume(lb, ub)

    ub = [3]
    @test_throws "Points have different dimensions" BoundingVolume(lb, ub)
end

@testset "Empty BV" begin
    bv = BoundingVolume()

    @test bv.is_empty
    @test bv.dim == 0
    @test bv.lb[1] == Inf
    @test bv.ub[1] == -Inf
    @test length(bv.active_dim) == 0
    @test length(bv.inactive_dim) == 0
    @test length(bv.is_active) == 0
end

@testset "Point BV" begin
    bv = BoundingVolume([1, 2, 3], [1, 2, 3])

    @test bv.is_empty == false
    @test bv.lb == bv.ub
    @test bv.dim == 0
    @test bv.active_dim == []
    @test bv.inactive_dim == [1, 2, 3]
    @test bv.is_active == [false, false, false]
end

@testset "Low-dimension BV" begin
    bv = BoundingVolume([1, 2, 3], [4, 2, 5])

    @test bv.is_empty == false
    @test bv.dim == 2
    @test bv.active_dim == [1, 3]
    @test bv.inactive_dim == [2]
    @test bv.is_active == [true, false, true]
    @test bv.lb[bv.inactive_dim] == bv.ub[bv.inactive_dim]
end

@testset "Full-dimension BV" begin
    bv = BoundingVolume([1, 2, 3], [4, 5, 6])

    @test bv.dim == 3
    @test bv.is_empty == false
    @test bv.active_dim == [1, 2, 3]
    @test length(bv.inactive_dim) == 0
    @test bv.is_active == ones(Bool, 3)
end

# `get_closest_point`: ------------------------------------------------------
@testset "get_closest_point(BV, pt): Interior Point (pt in BV)" begin
    bv = BoundingVolume([0, 0], [1, 1])
    pt = [0.5, 0.5]

    @test all(get_closest_point(bv, pt) .== pt)
end

@testset "get_closest_point(BV, pt): Boundary Point (pt on boundary of BV)" begin
    bv = BoundingVolume([0, 0], [1, 1])
    pt = [1, 0.5]

    @test all(get_closest_point(bv, pt) .== pt)
end

@testset "get_closest_point(BV, pt): Exterior Point (pt not in BV)" begin
    bv = BoundingVolume([0, 0], [1, 1])
    pt1 = [2, 2]
    pt2 = [-1, -0.5]

    @test all(get_closest_point(bv, pt1) .== bv.ub)
    @test all(get_closest_point(bv, pt2) .== bv.lb)
end

# `get_furthest_point`: ------------------------------------------------------
@testset "get_furthest_point(BV, pt): Interior Point (pt in BV)" begin
    bv = BoundingVolume([0, 0], [1, 1])
    pt1 = [0.5, 0.5]
    pt2 = [0.25, 0.25]

    @test all(get_furthest_point(bv, pt1) .== bv.lb) # Note: ties goes to lb
    @test all(get_furthest_point(bv, pt2) .== bv.ub)
end

@testset "get_furthest_point(BV, pt): Boundary Point (pt on boundary of BV)" begin
    bv = BoundingVolume([0, 0], [1, 1])
    pt = [0.25, 1]

    @test all(get_furthest_point(bv, pt) .== [1, 0])
end

@testset "get_furthest_point(BV, pt): Exterior Point (pt not in BV)" begin
    bv = BoundingVolume([0, 0], [1, 1])
    pt1 = [1.5, 1.5]
    pt2 = [-1, -0.5]

    @test all(get_furthest_point(bv, pt1) .== bv.lb)
    @test all(get_furthest_point(bv, pt2) .== bv.ub)
end

# `is_contained`: ----------------------------------------------------------
@testset "is_contained(BV, pt)" begin
    bv = BoundingVolume([0, 0], [1, 1])
    interior_pt = [0.5, 0.5]
    boundary_pt = [1, 0]
    exterior_pt = [2, 2]

    @test is_contained(bv, interior_pt, include_boundary=true) == true
    @test is_contained(bv, interior_pt, include_boundary=false) == true

    @test is_contained(bv, boundary_pt, include_boundary=true) == true
    @test is_contained(bv, boundary_pt, include_boundary=false) == false

    @test is_contained(bv, exterior_pt, include_boundary=true) == false
    @test is_contained(bv, exterior_pt, include_boundary=false) == false
end

@testset "is_contained(BV, BV): Empty Intersection" begin
    bv = BoundingVolume([0, 0], [1, 1])
    bv_query = BoundingVolume([-2, -2], [-1, -1])

    @test is_contained(bv, bv_query, include_boundary=true) == false
    @test is_contained(bv, bv_query, include_boundary=false) == false
end

@testset "is_contained(BV, BV): Partial Intersection" begin
    bv = BoundingVolume([0, 0], [1, 1])
    bv_query = BoundingVolume([-1, -1], [0.5, 0.5])

    # Full-dim intersection
    @test is_contained(bv, bv_query, include_boundary=true) == false
    @test is_contained(bv, bv_query, include_boundary=false) == false

    # Low-dim intersection
    bv_query = BoundingVolume([-1, 0], [0, 1])
    @test is_contained(bv, bv_query, include_boundary=true) == false
    @test is_contained(bv, bv_query, include_boundary=false) == false
end

@testset "is_contained(BV, BV): Contained" begin
    bv = BoundingVolume([0, 0], [1, 1])
    bv_query = BoundingVolume([0.25, 0.25], [0.75, 0.75])

    @test is_contained(bv, bv_query, include_boundary=true) == true
    @test is_contained(bv, bv_query, include_boundary=false) == true
end

@testset "is_contained(BV, BV): Strictly/Fully Contained" begin
    bv = BoundingVolume([0, 0], [1, 1])
    bv_query = BoundingVolume([0.25, 0.25], [1, 1])

    @test is_contained(bv, bv_query, include_boundary=true) == true
    @test is_contained(bv, bv_query, include_boundary=false) == false
end

# `intersects`: ----------------------------------------------------------
@testset "intersects(BV, BV): Empty Intersection" begin
    bv1 = BoundingVolume([1, 2], [3, 4])
    bv2 = BoundingVolume([-3, -4], [-1, -2])

    @test intersects(bv1, bv2, include_boundary=true) == false
    @test intersects(bv1, bv2, include_boundary=false) == false
end

@testset "intersects(BV, BV): Low-Dim Intersection" begin
    bv1 = BoundingVolume([1, 2], [3, 4])
    bv2 = BoundingVolume([0, 0], [1, 2])

    @test intersects(bv1, bv2, include_boundary=true) == true
    @test intersects(bv1, bv2, include_boundary=false) == false
end

@testset "intersects(BV, BV): Full-Dim Intersection" begin
    bv1 = BoundingVolume([1, 2], [3, 4])
    bv2 = BoundingVolume([0, 0], [2, 3])

    @test intersects(bv1, bv2, include_boundary=true) == true
    @test intersects(bv1, bv2, include_boundary=false) == true
end

# `get_intersection`: ------------------------------------------------------
@testset "get_intersection(BV, BV): Empty Intersection" begin
    bv1 = BoundingVolume([1, 2], [3, 4])
    bv2 = BoundingVolume([-3, -4], [-1, -2])

    intersection = get_intersection(bv1, bv2)
    @test intersection.is_empty == true
end

@testset "get_intersection(BV, BV): Low-Dim Intersection" begin
    bv = BoundingVolume([0, 0], [1, 1])

    bv_pt = BoundingVolume([-1, -1], [0, 0])
    pt_intersection = get_intersection(bv, bv_pt)
    @test pt_intersection.dim == 0
    @test pt_intersection.is_empty == false
    @test pt_intersection.lb == pt_intersection.ub

    bv_line = BoundingVolume([-1, 0], [0, 1])
    line_intersection = get_intersection(bv, bv_line)
    @test line_intersection.dim == 1
    @test line_intersection.is_empty == false
    @test line_intersection.lb == [0, 0]
    @test line_intersection.ub == [0, 1]
end

@testset "get_intersection(BV, BV): Full-Dim Intersection" begin
    bv = BoundingVolume([0, 0], [1, 1])

    bv_interior = BoundingVolume([0.25, 0.25], [0.75, 0.75])
    interior_intersection = get_intersection(bv, bv_interior)
    @test interior_intersection == bv_interior

    bv_overlapping = BoundingVolume([-1, -1], [0.25, 0.5])
    intersection_true = BoundingVolume([0, 0], [0.25, 0.5])
    intersection = get_intersection(bv, bv_overlapping)
    @test intersection == intersection_true
end

# `face_index_to_spatial_index`: ------------------------------------------------------
@testset "face_index_to_spatial_index(face_index, num_dim)" begin
    num_dim = 3
    @test face_index_to_spatial_index(1, num_dim) == 1
    @test face_index_to_spatial_index(2, num_dim) == 2
    @test face_index_to_spatial_index(3, num_dim) == 3
    @test face_index_to_spatial_index(4, num_dim) == 1
    @test face_index_to_spatial_index(5, num_dim) == 2
    @test face_index_to_spatial_index(6, num_dim) == 3
end

# `get_face_bounding_volume`: ------------------------------------------------------
@testset "get_face_bounding_volume(face_index, bv)" begin
    bv = BoundingVolume([0, 0, 0], [1, 1, 1])

    left = BoundingVolume([0, 0, 0], [0, 1, 1])
    right = BoundingVolume([1, 0, 0], [1, 1, 1])

    front = BoundingVolume([0, 0, 0], [1, 0, 1])
    back = BoundingVolume([0, 1, 0], [1, 1, 1])

    bottom = BoundingVolume([0, 0, 0], [1, 1, 0])
    top = BoundingVolume([0, 0, 1], [1, 1, 1])

    # Lower bound faces
    @test get_face_bounding_volume(1, bv) == left
    @test get_face_bounding_volume(2, bv) == front
    @test get_face_bounding_volume(3, bv) == bottom

    # Upper bound faces
    @test get_face_bounding_volume(4, bv) == right
    @test get_face_bounding_volume(5, bv) == back
    @test get_face_bounding_volume(6, bv) == top
end
