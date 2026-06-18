using SearchableGeometries
using SearchableGeometries.GeometricPrimitives: BoundingVolume
using SearchableGeometries: 
    partition_data!,
    partition_bv,
    is_base_case

# PointRange Constructors ----------------------------------------------------
@testset "Constructing PointRange: Valid" begin
    data = [
        0.0 1.0 4.0 5.0;
        0.0 2.0 3.0 1.0
    ]

    r = PointRange(data, 1; n = 4)
    @test r.data === data
    @test r.first == 1
    @test r.last == 4
    @test r.n == 4
    @test length(r) == 4
    @test !isempty(r)

    r2 = PointRange(data, 2, 3)
    @test r2.first == 2
    @test r2.last == 3
    @test r2.n == 2
    @test length(r2) == 2
end

@testset "Constructing PoinRange: Invalid" begin
    data = [
        0.0 1.0 4.0 5.0;
        0.0 2.0 3.0 1.0
    ]

    r = PointRange(data, 3, 2)
    @test r.n == 0
    @test isempty(r)

    @test_throws "SearchableGeometries.PointRange: n must be nonnegative" PointRange(data, 3; n = -1)
    @test_throws "SearchableGeometries.PointRange: last must be at least first - 1" PointRange(data, 3, 1)
end

# `BoundingVolume` -------------------------------------------------------
@testset "BoundingVolume(point_range)" begin
    data = [
        0.0  1.0  4.0  5.0;
        0.0  2.0  3.0  1.0
    ]

    r = PointRange(data, 1; n = 4)
    bv = BoundingVolume(r)

    @test bv.lb == [0.0, 0.0]
    @test bv.ub == [5.0, 3.0]

    r_sub = PointRange(data, 2, 3)
    bv_sub = BoundingVolume(r_sub)

    @test bv_sub.lb == [1.0, 2.0]
    @test bv_sub.ub == [4.0, 3.0]
end

# `get_split` ----------------------------------------------------------------
@testset "get_split(node): Matrix PointRange" begin
    data = [
        0.0  1.0  4.0  5.0;
        0.0  2.0  3.0  1.0
    ]

    r = PointRange(data, 1; n = 4)
    bv = BoundingVolume(r)

    node = Node(
        r,
        bv,
        Watertight;
        is_base_case = (val, bv, parent) -> true,
        get_split = get_split,
        partition_data! = partition_data!,
        partition_bv = partition_bv,
    )

    # Because the node was made a leaf immediately, manually test get_split
    split_dim, split_val = get_split(node)

    @test split_dim == 1
    @test isapprox(split_val, 2.5; atol=1e-6)
end

# `partition_data!` ----------------------------------------------------------
@testset "partition_data!(node): Matrix PointRange" begin
    data = [
        0.0  5.0  1.0  4.0;
        0.0  1.0  2.0  3.0
    ]

    r = PointRange(data, 1; n = 4)
    bv = BoundingVolume(r)

    node = Node(
        r,
        bv,
        Watertight;
        is_base_case = (val, bv, parent) -> true,
        get_split = get_split,
        partition_data! = partition_data!,
        partition_bv = partition_bv,
    )

    node.split_dim = 1
    node.split_val = 2.5

    val_L, val_R = partition_data!(node)

    @test val_L.first == 1
    @test val_R.last == 4
    @test val_L.n + val_R.n == 4

    for i in val_L.first:val_L.last
        @test data[node.split_dim, i] <= node.split_val
    end

    for i in val_R.first:val_R.last
        @test data[node.split_dim, i] > node.split_val
    end
end

# `Partition_bv` ---------------------------------------------------------
@testset "partition_bv(node, L_val, R_val): Watertight" begin
    data = [
        0.0  1.0  4.0  5.0;
        0.0  2.0  3.0  1.0
    ]

    r = PointRange(data, 1; n = 4)
    bv = BoundingVolume(r)

    node = Node(
        r,
        bv,
        Watertight;
        is_base_case = (val, bv, parent) -> true,
        get_split = get_split,
        partition_data! = partition_data!,
        partition_bv = partition_bv,
    )

    node.split_dim = 1
    node.split_val = 2.5

    val_L = PointRange(data, 1, 2)
    val_R = PointRange(data, 3, 4)

    bv_L, bv_R = partition_bv(node, val_L, val_R)

    @test bv_L.lb == bv.lb
    @test bv_L.ub == [2.5, 3.0]

    @test bv_R.lb == [2.5, 0.0]
    @test bv_R.ub == bv.ub
end

@testset "partition_bv(node, L_val, R_val): Arbitrary" begin
    data = [
        0.0  1.0  4.0  5.0;
        0.0  2.0  3.0  1.0
    ]

    r = PointRange(data, 1; n = 4)
    bv = BoundingVolume(r)

    node = Node(
        r,
        bv,
        Arbitrary;
        is_base_case = (val, bv, parent) -> true,
        get_split = get_split,
        partition_data! = partition_data!,
        partition_bv = partition_bv,
    )

    val_L = PointRange(data, 1, 2)
    val_R = PointRange(data, 3, 4)

    bv_L, bv_R = partition_bv(node, val_L, val_R)

    @test bv_L.lb == [0.0, 0.0]
    @test bv_L.ub == [1.0, 2.0]

    @test bv_R.lb == [4.0, 1.0]
    @test bv_R.ub == [5.0, 3.0]
end

# Node constructor --------------------------------------------------
@testset "Recursive Node construction with Watertight BVH" begin
    data = [
        0.0  1.0  4.0  5.0;
        0.0  2.0  3.0  1.0
    ]

    r = PointRange(data, 1; n = 4)
    bv = BoundingVolume(r)

    root = Node(
        r,
        bv,
        Watertight;
        is_base_case = (val, bv, parent) -> val.n <= 1,
        get_split = get_split,
        partition_data! = partition_data!,
        partition_bv = partition_bv,
    )

    @test root.parent === nothing
    @test !root.is_leaf
    @test root.left !== nothing
    @test root.right !== nothing

    @test root.left.parent === root
    @test root.right.parent === root

    @test root.left.val.n + root.right.val.n == root.val.n

    @test root.split_dim == 1
    @test isapprox(root.split_val, 2.5; atol=1e-6)

    # Watertight root split
    @test root.left.bv.lb == [0.0, 0.0]
    @test root.left.bv.ub == [2.5, 3.0]

    @test root.right.bv.lb == [2.5, 0.0]
    @test root.right.bv.ub == [5.0, 3.0]

    # Since leafsize is 1, the two children should split again.
    @test !root.left.is_leaf
    @test !root.right.is_leaf

    @test root.left.left !== nothing
    @test root.left.right !== nothing
    @test root.right.left !== nothing
    @test root.right.right !== nothing

    @test root.left.left.is_leaf
    @test root.left.right.is_leaf
    @test root.right.left.is_leaf
    @test root.right.right.is_leaf

    @test root.left.left.val.n == 1
    @test root.left.right.val.n == 1
    @test root.right.left.val.n == 1
    @test root.right.right.val.n == 1
end

@testset "Recursive Node construction with Arbitrary BVH" begin
    data = [
        0.0  1.0  4.0  5.0;
        0.0  2.0  3.0  1.0
    ]

    r = PointRange(data, 1; n = 4)
    bv = BoundingVolume(r)

    root = Node(
        r,
        bv,
        Arbitrary;
        is_base_case = (val, bv, parent) -> val.n <= 1,
        get_split = get_split,
        partition_data! = partition_data!,
        partition_bv = partition_bv,
    )

    @test root.parent === nothing
    @test !root.is_leaf
    @test root.left !== nothing
    @test root.right !== nothing

    @test root.left.parent === root
    @test root.right.parent === root

    @test root.left.val.n + root.right.val.n == root.val.n

    @test root.split_dim == 1
    @test isapprox(root.split_val, 2.5; atol=1e-6)

    # Arbitrary BVH uses tight child bounding volumes.
    @test root.left.bv.lb == [0.0, 0.0]
    @test root.left.bv.ub == [1.0, 2.0]

    @test root.right.bv.lb == [4.0, 1.0]
    @test root.right.bv.ub == [5.0, 3.0]

    # These should differ from the Watertight split.
    @test root.left.bv.ub[1] != root.split_val
    @test root.right.bv.lb[1] != root.split_val

    @test root.left.bv.ub[1] == 1.0
    @test root.right.bv.lb[1] == 4.0

    # Since leafsize is 1, the two children should split again.
    @test !root.left.is_leaf
    @test !root.right.is_leaf

    @test root.left.left !== nothing
    @test root.left.right !== nothing
    @test root.right.left !== nothing
    @test root.right.right !== nothing

    @test root.left.left.is_leaf
    @test root.left.right.is_leaf
    @test root.right.left.is_leaf
    @test root.right.right.is_leaf

    @test root.left.left.val.n == 1
    @test root.left.right.val.n == 1
    @test root.right.left.val.n == 1
    @test root.right.right.val.n == 1
end