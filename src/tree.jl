# ===========================================================================
# Tree family: you get a SpatialTreeInterface tree from a grid.
#
# Consumers accept grids, never bare trees, and ask the tree only traversal
# questions (children, extents, the positions its leaves emit). Cell count,
# polygons and manifold are always asked of the grid. A tree may hold whatever
# it needs internally, but the tree contract has no `getgrid`/`getcell`/`ncells`.
# ===========================================================================

"""
    treeify(grid)

A SpatialTreeInterface tree whose leaves emit `grid`'s positions, each exactly
once. A grid with analytic structure returns its own cursor; a grid that *is*
a tree returns itself. Generic method:

- a complete level grid or a subset of one (`system(grid) !== nothing`) gets the
  generic hierarchy cursor (O(1) to build when the system has sorted subtrees);
- anything else gets an R-tree over [`cell_extent`](@ref) (extension).
"""
treeify(grid) = _default_treeify(grid)

function _default_treeify(grid)
    if system(grid) !== nothing && level(grid) !== nothing
        return hierarchy_tree(grid)
    end
    return fallback_treeify(best_manifold(grid), grid)
end

"""
    subtree(grid, cells)

A tree over the cells at positions `cells` (sorted, unique) of `grid` whose
leaves emit **`grid`'s** positions. Generic method: the hierarchy cursor
restricted to `cells` for level grids, an R-tree over those cells otherwise.
"""
subtree(grid, cells::AbstractVector{<:Integer}) = _default_subtree(grid, cells)

function _default_subtree(grid, cells)
    if system(grid) !== nothing && level(grid) !== nothing && gridrelation(grid) isa NoParent
        return hierarchy_subtree(grid, cells)
    end
    return fallback_subtree(best_manifold(grid), grid, cells)
end

"""
    extent_is_expensive(grid) -> Bool

Traversal hint: whether [`cell_extent`](@ref) computes rather than reads.
Default `true`. (The node-level counterpart is `STI.node_extent_is_expensive`.)
"""
extent_is_expensive(grid) = true

"""
    split_weight(node) -> Int

Tree-node hint: a cheap estimate of the leaves under `node`, used to balance
dual traversals. Default `1`.
"""
split_weight(node) = 1

"""
    grid_extent(grid)

An extent covering every cell of `grid` (box or cap by manifold). Generic:
the root extent of `treeify(grid)`. Named `grid_extent` rather than `extent`
to avoid colliding with `Extents.extent`.
"""
grid_extent(grid) = STI.node_extent(treeify(grid))

"""
    query(grid, ext) -> Vector{Int}

Ascending positions of the cells whose [`cell_extent`](@ref) intersects
extent `ext`, found by descending `treeify(grid)`. Exact geometric filtering is
the caller's business.
"""
function query(grid, ext)
    out = Int[]
    tree = treeify(grid)
    STI.depth_first_search(e -> Extents.intersects(e, ext), tree) do i
        push!(out, i)
        nothing
    end
    return sort!(out)
end
