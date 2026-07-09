using SearchableGeometries
using SearchableGeometries.GeometricPrimitives: BoundingVolume
using SearchableGeometries: 
    partition_data!,
    partition_bv,
    is_leaf

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
    @test isvalid(r)

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
    @test !isvalid(r)

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

# `partition_data!` ----------------------------------------------------------
@testset "partition_data!(node): PointRange" begin
    data = [
        0.0  5.0  1.0  4.0;
        0.0  1.0  2.0  3.0
    ]

    r = PointRange(data, 1; n = 4)
    bv = BoundingVolume(r)

    node = kdTreeNode(
        r,
        bv,
        Watertight;
        is_leaf = (val) -> true,
        partition_data! = partition_data!
    )

    node.split_dim, node.split_val, val_L, val_R = partition_data!(node)

    @test node.split_dim == 1
    @test node.split_val == 2.5

    @test val_L == PointRange(data, 1, 2)
    @test val_R == PointRange(data, 3, 4)

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

    node = kdTreeNode(
        r,
        bv,
        Watertight;
        is_leaf = (val) -> true,
        partition_data! = partition_data!
    )

    node.split_dim, node.split_val, val_L, val_R = partition_data!(node)

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

    node = kdTreeNode(
        r,
        bv,
        Arbitrary;
        is_leaf = (val) -> true,
        partition_data! = partition_data!,
        partition_bv = partition_bv,
    )

    node.split_dim, node.split_val, val_L, val_R = partition_data!(node)

    bv_L, bv_R = partition_bv(node, val_L, val_R)

    @test bv_L.lb == [0.0, 0.0]
    @test bv_L.ub == [1.0, 2.0]

    @test bv_R.lb == [4.0, 1.0]
    @test bv_R.ub == [5.0, 3.0]
end

# kdTreeNode constructor --------------------------------------------------
@testset "Recursive kdTreeNode construction with Watertight BVH" begin
    data = [
        0.0  1.0  4.0  5.0;
        0.0  2.0  3.0  1.0
    ]

    r = PointRange(data, 1; n = 4)
    bv = BoundingVolume(r)

    root = kdTreeNode(
        r,
        bv,
        Watertight;
        is_leaf = (val) -> val.n <= 1,
        partition_data! = partition_data!
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

@testset "Recursive kdTreeNode construction with Arbitrary BVH" begin
    data = [
        0.0  1.0  4.0  5.0;
        0.0  2.0  3.0  1.0
    ]

    r = PointRange(data, 1; n = 4)
    bv = BoundingVolume(r)

    root = kdTreeNode(
        r,
        bv,
        Arbitrary;
        is_leaf = (val) -> val.n <= 1,
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

# leaf_search with one function -------------------------------------------
@testset "leaf_search(node, func)" begin
    data = [
        -4.0  -3.0  -4.0  -3.0   3.0   4.0   3.0   4.0;
        -4.0  -4.0  -3.0  -3.0   3.0   3.0   4.0   4.0
    ]

    r = PointRange(data, 1; n = 8)
    bv = BoundingVolume(r)

    root = kdTreeNode(
        r,
        bv,
        Arbitrary;
        is_leaf = val -> val.n <= 1,
        partition_data! = partition_data!,
        partition_bv = partition_bv,
    )

    @testset "returns all leaves when func is always true" begin
        leaves = leaf_search(root, node -> true)

        @test length(leaves) == 8
        @test all(node -> node.is_leaf, leaves)
        @test all(node -> node.val.n == 1, leaves)
    end

    @testset "returns empty when root fails func" begin
        leaves = leaf_search(root, node -> false)

        @test isempty(leaves)
    end

    @testset "prunes left subtree" begin
        leaves = leaf_search(root, node -> node === root)

        # Root is internal, so it is not returned.
        # Its children fail, so no leaves are reached.
        @test isempty(leaves)
    end

    @testset "does not return leaves that fail func" begin
        leaves = leaf_search(root, node -> !node.is_leaf)

        # Internal nodes may pass, but leaves fail,
        # so no leaf should be returned.
        @test isempty(leaves)
    end
end

# leaf_search with two functions ----------------------------------------
@testset "leaf_search(node, internal_func, leaf_func)" begin
    data = [
        -4.0  -3.0  -4.0  -3.0   3.0   4.0   3.0   4.0;
        -4.0  -4.0  -3.0  -3.0   3.0   3.0   4.0   4.0
    ]

    r = PointRange(data, 1; n = 8)
    bv = BoundingVolume(r)

    root = kdTreeNode(
        r,
        bv,
        Arbitrary;
        is_leaf = val -> val.n <= 1,
        partition_data! = partition_data!,
        partition_bv = partition_bv,
    )    

    @testset "returns all leaves when both functions allow search" begin
        leaves = leaf_search(root, node -> true, node -> true)

        @test length(leaves) == 8
        @test all(node -> node.is_leaf, leaves)
        @test all(node -> node.val.n == 1, leaves)
    end

    @testset "leaf_func controls returned leaves" begin
        leaves = leaf_search(
            root,
            node -> true,
            node -> false,
        )

        @test isempty(leaves)
    end

    @testset "internal_func controls pruning" begin
        leaves = leaf_search(
            root,
            node -> node === root,
            node -> true,
        )

        # Root passes, but its children are internal and fail.
        # Therefore we never reach the leaves.
        @test isempty(leaves)
    end

    @testset "internal_func is not applied to leaves" begin
        internal_count = 0
        leaf_count = 0

        leaves = leaf_search(
            root,
            node -> begin
                internal_count += 1
                return true
            end,
            node -> begin
                leaf_count += 1
                return true
            end,
        )

        @test length(leaves) == 8
        @test internal_count == 7
        @test leaf_count == 8
    end
end

# search ------------------------------------------------------------
@testset "search(node, func; shortcircuit=false)" begin
    data = [
        -4.0  -3.0  -4.0  -3.0   3.0   4.0   3.0   4.0;
        -4.0  -4.0  -3.0  -3.0   3.0   3.0   4.0   4.0
    ]

    r = PointRange(data, 1; n = 8)
    bv = BoundingVolume(r)

    root = kdTreeNode(
        r,
        bv,
        Arbitrary;
        is_leaf = val -> val.n <= 1,
        partition_data! = partition_data!,
        partition_bv = partition_bv,
    )  

    @testset "returns all nodes when func is always true" begin
        nodes = search(root, node -> true)

        @test length(nodes) == 15
        @test nodes[1] === root
        @test count(node -> node.is_leaf, nodes) == 8
        @test count(node -> !node.is_leaf, nodes) == 7
    end

    @testset "returns only leaves when func checks node.is_leaf" begin
        nodes = search(root, node -> node.is_leaf)

        @test length(nodes) == 8
        @test all(node -> node.is_leaf, nodes)
        @test all(node -> node.val.n == 1, nodes)
    end

    @testset "returns only internal nodes when func checks !node.is_leaf" begin
        nodes = search(root, node -> !node.is_leaf)

        @test length(nodes) == 7
        @test all(node -> !node.is_leaf, nodes)
    end

    @testset "shortcircuit=false still searches children when root fails" begin
        nodes = search(
            root,
            node -> node !== root;
            shortcircuit = false,
        )

        # Root fails and is not returned,
        # but all descendants are still searched.
        @test length(nodes) == 14
        @test !(root in nodes)
    end

    @testset "shortcircuit=true prunes when root fails" begin
        nodes = search(
            root,
            node -> node !== root;
            shortcircuit = true,
        )

        # Root fails, so the entire tree is pruned.
        @test isempty(nodes)
    end

    @testset "shortcircuit=true prunes failed left subtree" begin
        failed_node = root.left

        nodes = search(
            root,
            node -> node !== failed_node;
            shortcircuit = true,
        )

        # Root passes.
        # root.left fails, so its whole subtree is skipped.
        # root.right subtree has 7 nodes.
        # Total returned = root + right subtree = 8.
        @test length(nodes) == 8
        @test root in nodes
        @test !(failed_node in nodes)
        @test all(node -> node !== failed_node, nodes)
    end

    @testset "shortcircuit=false does not prune failed left subtree" begin
        failed_node = root.left

        nodes = search(
            root,
            node -> node !== failed_node;
            shortcircuit = false,
        )

        # failed_node is not returned,
        # but its descendants are still searched.
        # Total nodes = 15 - 1 = 14.
        @test length(nodes) == 14
        @test root in nodes
        @test !(failed_node in nodes)
    end

    @testset "preorder means root is returned first" begin
        nodes = search(root, node -> true)

        @test nodes[1] === root
    end
end