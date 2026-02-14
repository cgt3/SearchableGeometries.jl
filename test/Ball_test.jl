using SearchableGeometries

# Ball Constructors --------------------------------------------------------
@testset "Constructing Balls: Invalid radius" begin
    @test_throws "Cannot construct ball with negative radius" Ball([0, 0], -1)
    
end

@testset "Constructing Balls: Zero Radius" begin
    ball = Ball([0, 0], 0)

    @test ball.radius == 0
    @test ball.dim == 0
    @test ball.embedding_dim == 2
    @test all(ball.active_dim .== [])
    @test all(ball.inactive_dim .== [1, 2])
    @test all(ball.is_active .== [false, false])
end

@testset "Constructing Balls: Low-Dimension Ball" begin
    ball = Ball([0, 0, 0], 1, p=2, active_indices=true, indices=[1, 3])
    @test ball.radius == 1
    @test ball.p == 2
    @test ball.dim == 2
    @test ball.embedding_dim == 3
    @test all(ball.active_dim .== [1, 3])
    @test all(ball.inactive_dim .== [2])
    @test all(ball.is_active .== [true, false, true])

    ball = Ball([0, 0, 0], 1, p=2, active_indices=false, indices=[2])
    @test ball.radius == 1
    @test ball.p == 2
    @test ball.dim == 2
    @test ball.embedding_dim == 3
    @test all(ball.active_dim .== [1, 3])
    @test all(ball.inactive_dim .== [2])
    @test all(ball.is_active .== [true, false, true])
end

@testset "Constructing Balls: Full-Dim, Euclidean Distance" begin
    ball = Ball([0, 0, 0], 1)

    @test ball.radius == 1
    @test ball.p == 2
    @test ball.dim == 3
    @test ball.embedding_dim == 3
    @test all(ball.active_dim .== [1, 2, 3])
    @test all(ball.inactive_dim .== [])
    @test all(ball.is_active .== [true, true, true])
end

@testset "Constructing Balls: Full-Dim, arbitrary p" begin
    ball = Ball([0, 0, 0], 1, p=Inf)

    @test ball.radius == 1
    @test ball.p == Inf
    @test ball.dim == 3
    @test ball.embedding_dim == 3
    @test all(ball.active_dim .== [1, 2, 3])
    @test all(ball.inactive_dim .== [])
    @test all(ball.is_active .== [true, true, true])
end


# Ball -> BV: --------------------------------------------------------------

# `isContained(Ball, pt)` --------------------------------------------------
@testset "isContained(Ball, pt):" begin
    ball = Ball([0, 0, 0], 1)
    pt = [0, 0]

    @test_throws "Point dimension($(length(pt))) does not match ball embedding dimension($(ball.embedding_dim))" isContained(ball, pt)
end

@testset "isContained(Ball, Pt): Zero-Dim Ball" begin
    ball_p1 = Ball([0, 0], 0, p=1)
    ball_p2 = Ball([0, 0], 0, p=2)
    ball_pInf = Ball([0, 0], 0, p=Inf)

    pt = [0, 0]

    @test isContained(ball_p1, pt, include_boundary=true)
    @test isContained(ball_p2, pt, include_boundary=true)
    @test isContained(ball_pInf, pt, include_boundary=true)

    @test !isContained(ball_p1, pt, include_boundary=false)
    @test !isContained(ball_p2, pt, include_boundary=false)
    @test !isContained(ball_pInf, pt, include_boundary=false)
end

@testset "isContained(Ball, pt): Low-Dim Ball - Boundary point" begin
    ball_p1 = Ball([0, 0, 0], 1, p=1, active_indices=true, indices=[1, 2])
    ball_p2 = Ball([0, 0, 0], 1, p=2, active_indices=true, indices=[1, 2])
    ball_pInf = Ball([0, 0, 0], 1, p=Inf, active_indices=true, indices=[1, 2])
    
    pt = [1, 0, 0]

    @test isContained(ball_p1, pt, include_boundary=true)
    @test isContained(ball_p2, pt, include_boundary=true)
    @test isContained(ball_pInf, pt, include_boundary=true)

    @test !isContained(ball_p1, pt, include_boundary=false)
    @test !isContained(ball_p2, pt, include_boundary=false)
    @test !isContained(ball_pInf, pt, include_boundary=false)
end

@testset "isContained(BV, pt): Low-Dim Ball - Exterior point" begin
    ball_p1 = Ball([0, 0, 0], 1, p=1, active_indices=true, indices=[1, 2])
    ball_p2 = Ball([0, 0, 0], 1, p=2, active_indices=true, indices=[1, 2])
    ball_pInf = Ball([0, 0, 0], 1, p=Inf, active_indices=true, indices=[1, 2])
    
    # Out-of-plane point
    pt = [0, 0, 1]

    @test !isContained(ball_p1, pt, include_boundary=true)
    @test !isContained(ball_p2, pt, include_boundary=true)
    @test !isContained(ball_pInf, pt, include_boundary=true)

    @test !isContained(ball_p1, pt, include_boundary=false)
    @test !isContained(ball_p2, pt, include_boundary=false)
    @test !isContained(ball_pInf, pt, include_boundary=false)

    # In-plane exterior point
    pt = [2, 2, 0]

    @test !isContained(ball_p1, pt, include_boundary=true)
    @test !isContained(ball_p2, pt, include_boundary=true)
    @test !isContained(ball_pInf, pt, include_boundary=true)

    @test !isContained(ball_p1, pt, include_boundary=false)
    @test !isContained(ball_p2, pt, include_boundary=false)
    @test !isContained(ball_pInf, pt, include_boundary=false)
end

@testset "isContained(BV, pt): Full-Dim - Interior point" begin
    ball_p1 = Ball([0, 0], 1, p=1)
    ball_p2 = Ball([0, 0], 1, p=2)
    ball_pInf = Ball([0, 0], 1, p=Inf)
    
    pt = [0.25, 0.25]

    @test isContained(ball_p1, pt, include_boundary=true)
    @test isContained(ball_p2, pt, include_boundary=true)
    @test isContained(ball_pInf, pt, include_boundary=true)

    @test isContained(ball_p1, pt, include_boundary=false)
    @test isContained(ball_p2, pt, include_boundary=false)
    @test isContained(ball_pInf, pt, include_boundary=false)
end

@testset "isContained(BV, pt): Full-Dim - Boundary point" begin
    ball_p1 = Ball([0, 0], 1, p=1)
    ball_p2 = Ball([0, 0], 1, p=2)
    ball_pInf = Ball([0, 0], 1, p=Inf)
    
    pt = [1, 0]

    @test isContained(ball_p1, pt, include_boundary=true)
    @test isContained(ball_p2, pt, include_boundary=true)
    @test isContained(ball_pInf, pt, include_boundary=true)

    @test !isContained(ball_p1, pt, include_boundary=false)
    @test !isContained(ball_p2, pt, include_boundary=false)
    @test !isContained(ball_pInf, pt, include_boundary=false)
end
    
@testset "isContained(Ball, pt): Full-Dim Ball - Exterior point" begin
    ball_p1 = Ball([0, 0], 1, p=1)
    ball_p2 = Ball([0, 0], 1, p=2)
    ball_pInf = Ball([0, 0], 1, p=Inf)
    
    pt = [2, 2]

    @test !isContained(ball_p1, pt, include_boundary=true)
    @test !isContained(ball_p2, pt, include_boundary=true)
    @test !isContained(ball_pInf, pt, include_boundary=true)

    @test !isContained(ball_p1, pt, include_boundary=false)
    @test !isContained(ball_p2, pt, include_boundary=false)
    @test !isContained(ball_pInf, pt, include_boundary=false)
end

# `isContained(BV, Ball)` --------------------------------------------------
@testset "isContained(BV, Ball): Interior Ball" begin
    bv = BoundingVolume([0, 0], [1, 1])
    # Ball centered at [0.5, 0.5] with radius 0.25 -> bounds [0.25, 0.25], [0.75, 0.75]. STRICTLY INSIDE.
    ball_p1 = Ball([0.5, 0.5], 0.25, p=1)
    ball_p2 = Ball([0.5, 0.5], 0.25, p=2)
    ball_pInf = Ball([0.5, 0.5], 0.25, p=Inf)

    @test isContained(bv, ball_p1, include_boundary=true)
    @test isContained(bv, ball_p2, include_boundary=true)
    @test isContained(bv, ball_pInf, include_boundary=true)

    @test isContained(bv, ball_p1, include_boundary=false)
    @test isContained(bv, ball_p2, include_boundary=false)
    @test isContained(bv, ball_pInf, include_boundary=false)
end

@testset "isContained(BV, Ball): Boundary Ball" begin
    bv = BoundingVolume([0, 0], [1, 1])
    # Ball centered at [0.5, 0.5] with radius 0.5 -> bounds [0.0, 0.0], [1.0, 1.0]. MATCHES BV BOUNDARY.
    # This ball "meets the boundary" everywhere.
    ball_p1 = Ball([0.5, 0.5], 0.5, p=1)
    ball_p2 = Ball([0.5, 0.5], 0.5, p=2)
    ball_pInf = Ball([0.5, 0.5], 0.5, p=Inf)

    @test isContained(bv, ball_p1, include_boundary=true)
    @test isContained(bv, ball_p2, include_boundary=true)
    @test isContained(bv, ball_pInf, include_boundary=true)

    @test !isContained(bv, ball_p1, include_boundary=false)
    @test !isContained(bv, ball_p2, include_boundary=false)
    @test !isContained(bv, ball_pInf, include_boundary=false)
end

@testset "isContained(BV, Ball): Partial Overlap Ball" begin
    bv = BoundingVolume([0, 0], [1, 1])
    # Ball centered at [1, 0.5] with radius 0.5 -> bounds [0.5, 0], [1.5, 1]. INTERSECTS BUT OVERFLOWS.
    ball_p1 = Ball([1, 0.5], 0.5, p=1) 
    ball_p2 = Ball([1, 0.5], 0.5, p=2)
    ball_pInf = Ball([1, 0.5], 0.5, p=Inf)

    @test !isContained(bv, ball_p1, include_boundary=true)
    @test !isContained(bv, ball_p2, include_boundary=true)
    @test !isContained(bv, ball_pInf, include_boundary=true)

    @test !isContained(bv, ball_p1, include_boundary=false)
    @test !isContained(bv, ball_p2, include_boundary=false)
    @test !isContained(bv, ball_pInf, include_boundary=false)
end

@testset "isContained(BV, Ball): Exterior Ball" begin
    bv = BoundingVolume([0, 0], [1, 1])
    # Ball centered at [2, 2] with radius 0.25 -> bounds [1.75, 1.75], [2.25, 2.25]. OUTSIDE.
    ball_p1 = Ball([2, 2], 0.25, p=1)
    ball_p2 = Ball([2, 2], 0.25, p=2)
    ball_pInf = Ball([2, 2], 0.25, p=Inf)

    @test !isContained(bv, ball_p1, include_boundary=true)
    @test !isContained(bv, ball_p2, include_boundary=true)
    @test !isContained(bv, ball_pInf, include_boundary=true)

    @test !isContained(bv, ball_p1, include_boundary=false)
    @test !isContained(bv, ball_p2, include_boundary=false)
    @test !isContained(bv, ball_pInf, include_boundary=false)
end


# # `isContained(Ball, BV)` --------------------------------------------------
# @testset "isContained(Ball, BV): Interior Ball" begin
#     bv = BoundingVolume([0, 0], [1, 1])
#     # Ball centered at 0.5 with radius 0.25 -> bounds [0.25, 0.75]. STRICTLY INSIDE.
#     ball_p1 = Ball([0.5, 0.5], 0.25, p=1)
#     ball_p2 = Ball([0.5, 0.5], 0.25, p=2)
#     ball_pInf = Ball([0.5, 0.5], 0.25, p=Inf)

#     @test isContained(ball_p1, bv, include_boundary=true)
#     @test isContained(ball_p2, bv, include_boundary=true)
#     @test isContained(ball_pInf, bv, include_boundary=true)

#     @test isContained(ball_p1, bv, include_boundary=false)
#     @test isContained(ball_p2, bv, include_boundary=false)
#     @test isContained(ball_pInf, bv, include_boundary=false)
# end

# @testset "isContained(Ball, BV): Boundary Ball" begin
#     bv = BoundingVolume([0, 0], [1, 1])
#     # Ball centered at 0.5 with radius 0.5 -> bounds [0.0, 1.0]. MATCHES BV BOUNDARY.
#     # This ball "meets the boundary" everywhere.
#     ball_p1 = Ball([0.5, 0.5], 0.5, p=1)
#     ball_p2 = Ball([0.5, 0.5], 0.5, p=2)
#     ball_pInf = Ball([0.5, 0.5], 0.5, p=Inf)

#     @test isContained(ball_p1, bv, include_boundary=true)
#     @test isContained(ball_p2, bv, include_boundary=true)
#     @test isContained(ball_pInf, bv, include_boundary=true)

#     @test !isContained(ball_p1, bv, include_boundary=false)
#     @test !isContained(ball_p2, bv, include_boundary=false)
#     @test !isContained(ball_pInf, bv, include_boundary=false)
# end

# @testset "isContained(Ball, BV): Partial Overlap Ball" begin
#     bv = BoundingVolume([0, 0], [1, 1])
#     # Ball centered at 1.0 with radius 0.25 -> bounds [0.5, 1.5]. INTERSECTS BUT OVERFLOWS.
#     ball_p1 = Ball([1.0, 0.5], 0.5, p=1) 
#     ball_p2 = Ball([1.0, 0.5], 0.5, p=2)
#     ball_pInf = Ball([1.0, 0.5], 0.5, p=Inf)

#     @test !isContained(ball_p1, bv, include_boundary=true)
#     @test !isContained(ball_p2, bv, include_boundary=true)
#     @test !isContained(ball_pInf, bv, include_boundary=true)

#     @test !isContained(ball_p1, bv, include_boundary=false)
#     @test !isContained(ball_p2, bv, include_boundary=false)
#     @test !isContained(ball_pInf, bv, include_boundary=false)
# end

# @testset "isContained(Ball, BV): Exterior Ball" begin
#     bv = BoundingVolume([0, 0], [1, 1])
#     # Ball centered at 2.0 with radius 0.25 -> bounds [1.75, 2.25]. OUTSIDE.
#     ball_p1 = Ball([2.0, 0.5], 0.25, p=1)
#     ball_p2 = Ball([2.0, 0.5], 0.25, p=2)
#     ball_pInf = Ball([2.0, 0.5], 0.25, p=Inf)

#     @test !isContained(ball_p1, bv, include_boundary=true)
#     @test !isContained(ball_p2, bv, include_boundary=true)
#     @test !isContained(ball_pInf, bv, include_boundary=true)

#     @test !isContained(ball_p1, bv, include_boundary=false)
#     @test !isContained(ball_p2, bv, include_boundary=false)
#     @test !isContained(ball_pInf, bv, include_boundary=false)
# end