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
    @test bv.active_dim == [ ]
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

# `getClosestPoint`: ------------------------------------------------------
@testset "getClosestPoint(BV, pt): Interior Point (pt in BV)" begin
    bv = BoundingVolume([0, 0], [1, 1])
    pt = [0.5, 0.5]
    
    @test all(getClosestPoint(bv, pt) .== pt)
end

@testset "getClosestPoint(BV, pt): Boundary Point (pt on boundary of BV)" begin
    bv = BoundingVolume([0, 0], [1, 1])
    pt = [1, 0.5]
    
    @test all(getClosestPoint(bv, pt) .== pt)
end

@testset "getClosestPoint(BV, pt): Exterior Point (pt not in BV)" begin
    bv = BoundingVolume([0, 0], [1, 1])
    pt1 = [2, 2]
    pt2 = [-1, -0.5]
    
    @test all(getClosestPoint(bv, pt1) .== bv.ub)
    @test all(getClosestPoint(bv, pt2) .== bv.lb)
end

# `getFurthestPoint`: ------------------------------------------------------
@testset "getFurthestPoint(BV, pt): Interior Point (pt in BV)" begin
    bv = BoundingVolume([0, 0], [1, 1])
    pt1 = [0.5, 0.5]
    pt2 = [0.25, 0.25]
    
    @test all(getFurthestPoint(bv, pt1) .== bv.lb) # Note: ties goes to lb
    @test all(getFurthestPoint(bv, pt2) .== bv.ub)
end

@testset "getFurthestPoint(BV, pt): Boundary Point (pt on boundary of BV)" begin
    bv = BoundingVolume([0, 0], [1, 1])
    pt = [0.25, 1]
    
    @test all(getFurthestPoint(bv, pt) .== [1, 0])
end

@testset "getFurthestPoint(BV, pt): Exterior Point (pt not in BV)" begin
    bv = BoundingVolume([0, 0], [1, 1])
    pt1 = [1.5, 1.5]
    pt2 = [-1, -0.5]
    
    @test all(getFurthestPoint(bv, pt1) .== bv.lb)
    @test all(getFurthestPoint(bv, pt2) .== bv.ub)
end