struct NotAScheme end   # an id scheme no grid declares

"""
    check_grid(grid; name, nsamples = 64, full_tree_limit = 200_000) -> ConformanceReport

Run every family's laws that `grid` claims. Never throws.
"""
function check_grid(grid; name::AbstractString=string(nameof(typeof(grid))), nsamples::Int=64,
        full_tree_limit::Int=200_000, area_rtol::Real=1e-3, total_area_limit::Int=2_000_000,
        quiet::Bool=true)
    rep = ConformanceReport(String(name))
    n = ncells(grid)
    m = best_manifold(grid)
    S = sample_positions(n, nsamples)
    append!(rep.fallbacks, fallback_verbs(grid))
    metric!(rep, "ncells", n)
    metric!(rep, "manifold", string(nameof(typeof(m))))
    metric!(rep, "cellshape", string(nameof(typeof(cellshape(grid)))))
    check_identity!(rep, grid, n, S)
    # spherical type contract first: every fallback refuses other types, so a
    # grid that breaks it would otherwise fail every family with the same error
    if m isa Spherical && n > 0 && !check_spherical_types!(rep, grid, S)
        skip!(rep, :tree, "spherical type contract broken (see geometry)")
        quiet || print_report(rep)
        return rep
    end
    tree = check_tree!(rep, grid, n, m, S, full_tree_limit)
    check_geometry!(rep, grid, m, S, area_rtol, total_area_limit)
    check_location!(rep, grid, m, S, tree)
    check_neighbors!(rep, grid, m, S)
    check_dual!(rep, grid, m, S)
    check_subsets!(rep, grid, m, S)
    check_grid_performance!(rep, grid, S)
    quiet || print_report(rep)
    return rep
end

"""
    test_grid(grid; kw...) -> ConformanceReport

`check_grid`, emitted as a `@testset` with one `@test` per law.
"""
function test_grid(grid; kw...)
    rep = check_grid(grid; kw...)
    emit_tests(rep)
    return rep
end

# ------------------------------------------------------------------ identity
function check_identity!(rep, grid, n, S)
    F = :identity
    law!(rep, F, "ncells is a non-negative Int", n isa Int && n >= 0)
    n == 0 && return
    guarded!(rep, F, "localindex ∘ cellindex == id (canonical)") do
        bad = [i for i in S if localindex(grid, cellindex(grid, i)) != i]
        isempty(bad), "positions $(first(bad, 5))"
    end
    for T in cellindextypes(grid)
        guarded!(rep, F, "bijection through scheme $T") do
            bad = [i for i in S if localindex(grid, cellindex(grid, i, T)) != i]
            isempty(bad), "positions $(first(bad, 5))"
        end
        guarded!(rep, F, "reindex round trip $T") do
            all(i -> reindex(first(cellindextypes(grid)), grid, cellindex(grid, i, T)) == cellindex(grid, i), S)
        end
    end
    guarded!(rep, F, "BoundsError at 0 and n+1") do
        e0 = try cellindex(grid, 0); nothing catch e; e end
        e1 = try cellindex(grid, n + 1); nothing catch e; e end
        e0 isa BoundsError && e1 isa BoundsError, "got $(typeof(e0)), $(typeof(e1))"
    end
    guarded!(rep, F, "undeclared scheme refused up front") do
        e = try cellindex(grid, 1, NotAScheme); nothing catch e; e end
        e isa ArgumentError, "got $(typeof(e))"
    end
    # nothing for non-members
    nonmember = _nonmember(grid, n)
    if nonmember === nothing
        info!(rep, F, "nothing for non-members", "no non-member id constructible generically")
    else
        guarded!(rep, F, "localindex(non-member) === nothing") do
            localindex(grid, nonmember) === nothing, "for $nonmember"
        end
    end
    sys = GrI.system(grid)
    if sys !== nothing
        guarded!(rep, F, "level(c) on every sampled id") do
            all(i -> GrI.level(cellindex(grid, i)) isa Integer, S)
        end
    end
    # localindex answers without a search (measured)
    c = cellindex(grid, n)
    localindex(grid, c)
    t = @elapsed for _ in 1:100; localindex(grid, c); end
    metric!(rep, "localindex(last) ns", round(t / 100 * 1e9; digits=1))
    generic = :localindex in rep.fallbacks && !(c isa Integer)
    law!(rep, F, "localindex is not the linear-scan fallback for typed ids", !generic,
        generic ? "generic scan answers ($(round(t/100*1e9))) ns at position $n" : "")
end

function _nonmember(grid, n)
    rel = GrI.gridrelation(grid)
    if rel isa GrI.ParentGrid
        P = parent(grid)
        members = Set(cellindex(grid, i) for i in 1:min(n, 100_000))
        for j in 1:min(ncells(P), 100_000)
            c = cellindex(P, j)
            c in members || return c
        end
        return nothing
    end
    sys = GrI.system(grid)
    l = GrI.level(grid)
    if sys !== nothing && l !== nothing
        for l2 in GrI.levels(sys)
            l2 == l && continue
            return cellindex(GrI.levelgrid(sys, l2), 1)
        end
    end
    c = cellindex(grid, 1)
    c isa Integer && return n + 1
    return nothing
end

# ------------------------------------------------------------------ spherical types
# Points on a spherical manifold are `UnitSphericalPoint`, extents `SphericalCap`.
function check_spherical_types!(rep, grid, S)
    F = :geometry
    ok = guarded!(rep, F, "spherical points are UnitSphericalPoint") do
        bad = Any[]
        for i in S
            pts = cellshape(grid) isa SingleRing ? cell_boundary(grid, i) : all_points(cell_geometry(grid, i))
            for p in pts
                p isa UnitSphericalPoint || (push!(bad, (i, typeof(p))); break)
            end
        end
        isempty(bad), isempty(bad) ? "" : "cell $(bad[1][1]) has $(bad[1][2]); convert with GeometryOps.UnitSphereFromGeographic()"
    end
    ok || return false
    guarded!(rep, F, "spherical centroid is UnitSphericalPoint, extent is SphericalCap") do
        bad = [i for i in S if !(cell_centroid(grid, i) isa UnitSphericalPoint && cell_extent(grid, i) isa SphericalCap)]
        isempty(bad), isempty(bad) ? "" : "cell $(bad[1]): $(typeof(cell_centroid(grid, bad[1]))), $(typeof(cell_extent(grid, bad[1])))"
    end
    return true
end

# ------------------------------------------------------------------ tree
function check_tree!(rep, grid, n, m, S, full_tree_limit)
    F = :tree
    tree = try
        treeify(grid)
    catch e
        law!(rep, F, "treeify", false, "threw: " * first(sprint(showerror, e), 300))
        return nothing
    end
    law!(rep, F, "treeify returns an STI tree", STI.isspatialtree(tree), string(typeof(tree)))
    metric!(rep, "tree type", string(nameof(typeof(tree))))
    if n <= full_tree_limit
        guarded!(rep, F, "leaves emit exactly 1:ncells once each") do
            ls = sort!(collect_leaves(tree))
            ls == collect(1:n), "got $(length(ls)) leaves, $(length(unique(ls))) unique, range $(isempty(ls) ? "∅" : extrema(ls))"
        end
    else
        info!(rep, F, "leaves emit exactly 1:ncells once each", "skipped: ncells = $n > $full_tree_limit")
    end
    # extents: node contains its leaves' geometry; node contains child extents (PRD wording)
    nested_bad = Ref(0); geom_bad = Ref(0); visited = Ref(0)
    guarded!(rep, F, "every node's extent contains its leaves' geometry") do
        _walk_extents!(grid, m, tree, nested_bad, geom_bad, visited, 400)
        geom_bad[] == 0, "$(geom_bad[]) violations over $(visited[]) nodes"
    end
    # Only R-trees promise nested extents; a cursor's node extent covers its
    # descendants' geometry (law above) but a child's cap may poke out of it.
    if nameof(typeof(tree)) === :RTree
        law!(rep, F, "every node's extent contains its children's extents (R-tree)",
            nested_bad[] == 0, "$(nested_bad[]) non-nested parent/child extent pairs over $(visited[]) nodes")
    else
        info!(rep, F, "nested child extents (promised by R-trees only)",
            "$(nested_bad[]) non-nested pairs over $(visited[]) nodes")
    end
    guarded!(rep, F, "subtree(grid, cells) emits parent positions") do
        cells = S[1:2:end]
        ls = sort!(collect_leaves(subtree(grid, cells)))
        ls == cells, "expected $(length(cells)) got $(length(ls))"
    end
    # the intents: every hint and an integer yield a tree over 1:ncells
    for hint in (GrI.Balanced(), GrI.ForSweeps(), GrI.ForQueries(), 3)
        guarded!(rep, F, "treeify(grid; leaf_capacity = $hint) is a tree over 1:ncells") do
            t = treeify(grid; leaf_capacity=hint)
            ok = STI.isspatialtree(t)
            n <= full_tree_limit && (ok &= sort!(collect_leaves(t)) == 1:n)
            ok, string(typeof(t))
        end
    end
    # GridInterface's own trees honour the capacity; a provider with one tree
    # shape may ignore it (reported, not failed)
    cap = GrI.leaf_capacity(grid, GrI.ForSweeps())
    ts = treeify(grid; leaf_capacity=GrI.ForSweeps())
    worst = _max_leaf(ts, 2000)
    metric!(rep, "ForSweeps max leaf size", worst)
    if nameof(typeof(ts)) in (:HierarchyCursor, :SelectionCursor, :RTree)
        law!(rep, F, "ForSweeps leaves hold at most leaf_capacity(grid, ForSweeps()) = $cap cells", worst <= cap,
            "a leaf holds $worst")
    else
        info!(rep, F, "ForSweeps leaves vs leaf_capacity = $cap (provider tree)", "largest leaf $worst")
    end
    # wrappers: GridAndTree and mockchunks forward every verb in WRAPPED_VERBS
    nbr_ok = cellshape(grid) isa GrI.Polygonal || !(:neighbors in rep.fallbacks)
    for (wname, w) in (("GridAndTree", GrI.GridAndTree(grid, tree)), ("mockchunks", GrI.mockchunks(grid, GrI.chunks(grid))))
        guarded!(rep, F, "$wname forwards every verb in WRAPPED_VERBS") do
            bad = String[]; k = 0
            for (f, Ts) in GrI.WRAPPED_VERBS
                args = _wrapper_args(grid, S, f, Ts, nbr_ok)
                args === nothing && continue
                k += 1
                fn = getfield(GrI, f)
                _wrapper_key(f, fn(w, args...)) == _wrapper_key(f, fn(grid, args...)) || push!(bad, "$f$Ts")
            end
            metric!(rep, "$wname verbs exercised", "$k / $(length(GrI.WRAPPED_VERBS))")
            isempty(bad), "differ: " * join(first(bad, 6), ", ")
        end
    end
    guarded!(rep, F, "GridAndTree returns its tree for every hint and an integer") do
        g = GrI.GridAndTree(grid, tree)
        GrI.treeify(g) === tree && all(h -> treeify(g; leaf_capacity=h) === tree, (GrI.Balanced(), GrI.ForSweeps(), GrI.ForQueries(), 3))
    end
    guarded!(rep, F, "split_weight and subtree_indices survive GridAndTree") do
        g = GrI.GridAndTree(grid, tree)
        GrI.split_weight(g) == GrI.split_weight(tree) && collect(GrI.subtree_indices(g)) == collect(GrI.subtree_indices(tree))
    end
    guarded!(rep, F, "mockchunks answers the given chunks and chunkat agrees") do
        ch = GrI.chunks(grid)
        cells = S[1:2:end]
        mine = [GrI.subset(grid, cells), GrI.subset(grid, setdiff(1:n, cells))]
        w = GrI.mockchunks(grid, mine)
        GrI.chunks(w) === mine && GrI.chunkat(w, cells[1]) == 1 && GrI.chunkat(w, n) == (n in cells ? 1 : 2) &&
            (ch === nothing ? GrI.chunks(GrI.mockchunks(grid, nothing)) === nothing : true)
    end
    return tree
end

function _walk_extents!(grid, m, node, nested_bad, geom_bad, visited, budget)
    visited[] >= budget && return
    visited[] += 1
    ext = STI.node_extent(node)
    if STI.isleaf(node)
        for (i, e) in STI.child_indices_extents(node)
            pts = _cellpts(grid, i)
            all(p -> extent_contains_point(m, ext, p), pts) || (geom_bad[] += 1)
            all(p -> extent_contains_point(m, e, p), pts) || (geom_bad[] += 1)
        end
        return
    end
    for child in STI.getchild(node)
        ce = STI.node_extent(child)
        extent_contains_extent(m, ext, ce) || (nested_bad[] += 1)
        # leaf geometry under child must be inside this node's extent too
        if STI.isleaf(child)
            for (i, _) in STI.child_indices_extents(child)
                all(p -> extent_contains_point(m, ext, p), _cellpts(grid, i)) || (geom_bad[] += 1)
            end
        end
        _walk_extents!(grid, m, child, nested_bad, geom_bad, visited, budget)
        visited[] >= budget && return
    end
end

_cellpts(grid, i) = cellshape(grid) isa SingleRing ? collect(cell_boundary(grid, i)) : all_points(cell_geometry(grid, i))

# the largest leaf (entries) over the first `budget` nodes
function _max_leaf(tree, budget)
    worst = 0; seen = 0
    stack = Any[tree]
    while !isempty(stack) && seen < budget
        node = pop!(stack); seen += 1
        if STI.isleaf(node)
            worst = max(worst, length(STI.child_indices_extents(node)))
        else
            for c in STI.getchild(node); push!(stack, c); end
        end
    end
    return worst
end

# Sample arguments for a forwarded verb shape, or `nothing` when the suite does
# not exercise it (O(ncells) dense vertex verbs, verbs the grid cannot answer).
function _wrapper_args(grid, S, f, Ts, nbr_ok)
    i = S[end]
    nbr = (:neighbors, :ring, :halo, :border, :interior)
    f in (:dual, :nvertices, :vertex_table, :vertex_numbering, :vertexat, :cell_vertices, :cell_vertex_keys,
          :vertex_cells, :vertex_point, :vertex_isopen, :vertex_key, :vertex_index) && return nothing
    f in nbr && !nbr_ok && return nothing
    f in (:cell_boundary, :cell_corners) && !(cellshape(grid) isa SingleRing) && return nothing
    f === :parentpositions && !(GrI.gridrelation(grid) isa GrI.ParentGrid) && return nothing
    f === :globalindex && GrI.gridrelation(grid) isa GrI.ParentSystem && return nothing
    f === :chunkat && GrI.chunks(grid) === nothing && return nothing
    f === :chunks && Ts == (Integer,) && return nothing
    f === :localindex && Ts == (Integer,) && return cellindex(grid, i) isa Integer ? (i,) : nothing
    Ts == () && return ()
    Ts == (Integer,) && return f in (:halo, :maxring, :maxneighbors) ? (1,) : f === :treeify ? (3,) : (i,)
    Ts == (Integer, Type) && return (i, first(cellindextypes(grid)))
    Ts == (Integer, Integer) && return f in nbr ? (i, 1) : nothing
    Ts == (Integer, Any) && return (i, cell_centroid(grid, i))
    Ts == (Any,) && return f in (:localindex, :globalindex) ? (cellindex(grid, i),) :
                           f === :query ? (cell_extent(grid, i),) : (cell_centroid(grid, i),)
    Ts == (GrI.QueryPredicate,) && return (GrI.CentroidIn(cell_extent(grid, i)),)
    Ts == (GrI.Connectivity,) && return (Edge(),)
    Ts == (Integer, GrI.Connectivity) && return (1, Edge())
    Ts == (GrI.LeafHint,) && return f === :treeify ? (GrI.ForQueries(),) : (GrI.ForSweeps(),)
    Ts == (GrI.Balanced,) && return (GrI.Balanced(),)
    Ts == (GrI.ForSweeps,) && return (GrI.ForSweeps(),)
    Ts == (GrI.ForQueries,) && return (GrI.ForQueries(),)
    Ts == (Integer, GrI.Index{GrI.Local}) && return (i, GrI.Index(GrI.Local()))
    Ts == (Integer, GrI.Index{GrI.Global}) && return (i, GrI.Index(GrI.Global()))
    Ts == (Integer, GrI.Index{Type{T}} where {T}) && return (i, GrI.Index(first(cellindextypes(grid))))
    Ts == (AbstractVector{<:Integer},) && return (S[1:2:end],)
    Ts == (AbstractVector{<:Integer}, Integer) && return (S[1:2:end], 1)
    Ts == (GrI.Subset,) && return (GrI.subset(grid, S[1:2:end]),)
    return nothing
end
# a comparable key for a verb's answer (lazy vectors, trees, polygons)
function _wrapper_key(f, x)
    f === :chunks && return x === nothing ? nothing : length(x)
    STI.isspatialtree(x) && return (typeof(x), GrI.split_weight(x))
    f === :cell_geometry && return collect(GI.getpoint(x))
    x isa AbstractVector && return collect(x)
    return x
end

# ------------------------------------------------------------------ geometry
function check_geometry!(rep, grid, m, S, area_rtol, total_area_limit=2_000_000)
    F = :geometry
    sh = cellshape(grid)
    n = ncells(grid)
    guarded!(rep, F, "cell_geometry is polygonal") do
        all(i -> GI.trait(cell_geometry(grid, i)) isa Union{GI.PolygonTrait,GI.MultiPolygonTrait}, S)
    end
    if m isa Spherical
        guarded!(rep, F, "vertices lie on the unit sphere") do
            bad = count(i -> any(p -> abs(_norm(_v(p)) - 1) > 1e-10, _cellpts(grid, i)), S)
            bad == 0, "$bad cells"
        end
    end
    if sh isa SingleRing
        guarded!(rep, F, "ring is unclosed, ≥ 3 points, single ring") do
            bad = Int[]
            for i in S
                b = collect(cell_boundary(grid, i))
                ok = length(b) >= 3 && b[1] != b[end] && length(collect(GI.getring(cell_geometry(grid, i)))) == 1
                ok || push!(bad, i)
            end
            isempty(bad), "cells $(first(bad, 5))"
        end
        guarded!(rep, F, "ring is counter-clockwise") do
            bad = Int[]
            for i in S
                b = collect(cell_boundary(grid, i))
                ring_is_ccw(m, b, ccw_centre(m, b, cell_centroid(grid, i))) || push!(bad, i)
            end
            isempty(bad), "$(length(bad)) of $(length(S)) cells, e.g. $(first(bad, 5))"
        end
        guarded!(rep, F, "cell_corners: a subsequence of the ring, no consecutive repeats") do
            bad = Int[]
            for i in S
                cs = collect(GrI.cell_corners(grid, i)); b = collect(cell_boundary(grid, i))
                ok = length(cs) >= 3 && is_cyclic_subsequence(m, cs, b) &&
                     all(k -> !_close(m, cs[k], cs[mod1(k + 1, length(cs))]; tol=1e-12), eachindex(cs))
                ok || push!(bad, i)
            end
            isempty(bad), "cells $(first(bad, 5))"
        end
        guarded!(rep, F, "cell_boundary! fills a buffer with the ring") do
            buf = similar(collect(cell_boundary(grid, S[1])), 0)
            all(i -> collect(GrI.cell_boundary!(buf, grid, i)) == collect(cell_boundary(grid, i)), S)
        end
    else
        skip!(rep, F, "ring laws: grid declares $(nameof(typeof(sh)))")
    end
    if sh isa ConvexCells
        guarded!(rep, F, "ConvexCells: convex and inside an open hemisphere") do
            bad = [i for i in S if !ring_is_convex(m, collect(cell_boundary(grid, i)))]
            isempty(bad), "$(length(bad)) of $(length(S)) cells, e.g. $(first(bad, 5))"
        end
    end
    guarded!(rep, F, "centroid inside its cell") do
        bad = [i for i in S if !GrI.cell_contains(grid, i, cell_centroid(grid, i))]
        isempty(bad), "$(length(bad)) of $(length(S)), e.g. $(first(bad, 5))"
    end
    guarded!(rep, F, "area positive") do
        all(i -> cell_area(grid, i) > 0, S)
    end
    if !(:cell_area in rep.fallbacks) && GrI.edgekind(grid) isa GrI.CurvedEdges
        # the boundary only approximates curved edges: report the gap, do not fail
        worst = maximum(i -> (g = GrI.fallback_area(m, grid, i); abs(cell_area(grid, i) - g) / abs(g)), S)
        metric!(rep, "area max rel diff vs geometry", worst)
        info!(rep, F, "overridden area vs the densified boundary (CurvedEdges)", "max rel diff $worst")
    elseif !(:cell_area in rep.fallbacks)
        guarded!(rep, F, "overridden area equals the geometric answer (rtol $area_rtol)") do
            worst = 0.0
            for i in S
                a = cell_area(grid, i)
                g = GrI.fallback_area(m, grid, i)
                worst = max(worst, abs(a - g) / abs(g))
            end
            metric!(rep, "area max rel diff vs geometry", worst)
            worst <= area_rtol, "max rel diff $worst"
        end
    end
    guarded!(rep, F, "cellsize is positive and finite") do
        h = GrI.cellsize(grid)
        metric!(rep, "cellsize", h)
        h > 0 && isfinite(h), "got $h"
    end
    _check_total_area!(rep, grid, m, n, total_area_limit)
    guarded!(rep, F, "cell_extent contains the cell") do
        bad = [i for i in S if !all(p -> extent_contains_point(m, cell_extent(grid, i), p), _cellpts(grid, i))]
        isempty(bad), "cells $(first(bad, 5))"
    end
end

# A global tessellating grid on the sphere covers it: areas sum to 4πr².
function _check_total_area!(rep, grid, m, n, limit)
    F = :geometry
    m isa Spherical || return
    want = 4π * Float64(m.radius)^2
    if !(tessellation(grid) isa Tessellating)
        return
    elseif n > limit
        return info!(rep, F, "total area 4πr² (global grids)", "skipped: ncells = $n > $limit")
    end
    guarded!(rep, F, GrI.isglobal(grid) ? "total area equals 4πr² (global grid)" : "total area at most 4πr²") do
        A = 0.0
        for i in 1:n
            A += cell_area(grid, i)
        end
        metric!(rep, "total area / 4πr²", A / want)
        GrI.isglobal(grid) ? (isapprox(A, want; rtol=1e-6), "sum = $(A / want) × 4πr²") :
            (A <= want * (1 + 1e-6), "sum = $(A / want) × 4πr²")
    end
end

# ------------------------------------------------------------------ location
function check_location!(rep, grid, m, S, tree)
    F = :location
    tree === nothing && return skip!(rep, F, "no tree")
    analytic = GrI.location_is_analytic(grid)
    if analytic && :cellat in rep.fallbacks
        law!(rep, F, "declared analytic location is provided (not the fallback)", false,
            "location_is_analytic is true but cellat resolves to the generic fallback")
    end
    tess = tessellation(grid)
    pts = Tuple{Int,Any}[]
    for i in S
        c = cell_centroid(grid, i)
        push!(pts, (i, c))
        b = _cellpts(grid, i)
        for k in 1:min(length(b), 4)
            push!(pts, (i, _toward(m, c, b[k], 0.5)))
        end
    end
    multi = 0; agree_geo = 0; agree_n = 0; ties_bad = 0; centroid_bad = Int[]
    for (i, p) in pts
        cs = containing_cells(grid, tree, m, p)
        length(cs) > 1 && (multi += 1)
        geo = isempty(cs) ? nothing : first(cs)
        a = cellat(grid, p)
        agree_n += 1
        # PRD: on an overlapping grid an analytic answer wins over the lowest-position
        # tie rule, so there it only has to be one of the containing cells
        (a == geo || (analytic && tess isa Overlapping && a !== nothing && a in cs)) && (agree_geo += 1)
        if tess isa Overlapping && length(cs) > 1 && !analytic && a != first(cs)
            ties_bad += 1
        end
        if p === cell_centroid(grid, i) || p == cell_centroid(grid, i)
            tess isa Tessellating && a != i && push!(centroid_bad, i)
        end
    end
    metric!(rep, "points with >1 containing cell", "$multi / $(length(pts))")
    if tess isa Tessellating
        law!(rep, F, "at most one containing cell (tessellation law)", multi == 0,
            "$multi of $(length(pts)) interior sample points lie in more than one cell")
        law!(rep, F, "cellat(centroid of i) == i", isempty(centroid_bad), "cells $(first(centroid_bad, 5))")
    else
        info!(rep, F, "tessellation law", "grid declares Overlapping; $multi multiply-covered sample points")
        law!(rep, F, "tie rule (lowest position) on overlapping grid", ties_bad == 0, "$ties_bad violations")
    end
    if cellshape(grid) isa SingleRing
        # Points 1 % and 10 % of the centroid-to-edge distance inside each edge
        # midpoint (the tolerance); the 0.01 % rate is reported, not failed.
        guarded!(rep, F, "cellat(p) contains p, sampled near edges (within 1 %)") do
            bad = 0; tot = 0; close_bad = 0; close_tot = 0
            for i in S
                c = cell_centroid(grid, i)
                b = collect(GrI.cell_corners(grid, i))
                for k in eachindex(b)
                    mid = _toward(m, b[k], b[mod1(k + 1, length(b))], 0.5)
                    for t in (0.9, 0.99, 0.9999)
                        p = _toward(m, c, mid, t)
                        j = cellat(grid, p)
                        ok = j !== nothing && GrI.cell_contains(grid, j, p)
                        if t == 0.9999
                            close_tot += 1; ok || (close_bad += 1)
                        else
                            tot += 1; ok || (bad += 1)
                        end
                    end
                end
            end
            metric!(rep, "near-edge cellat disagreement (1 %, 10 %)", "$bad / $tot")
            metric!(rep, "near-edge cellat disagreement (0.01 %)", "$close_bad / $close_tot")
            bad == 0, "$bad of $tot near-edge points ($(round(100bad/tot; digits=2)) %)"
        end
    end
    if tess isa Overlapping
        guarded!(rep, F, "cellsat returns every containing cell") do
            all(((i, p),) -> collect(GrI.cellsat(grid, p)) == containing_cells(grid, tree, m, p), pts[1:min(end, 64)])
        end
    end
    law!(rep, F, analytic ? "analytic cellat agrees with the geometric answer (a containing cell when Overlapping)" : "cellat agrees with the geometric answer",
        agree_geo == agree_n, "$(agree_n - agree_geo) of $agree_n disagree")
end

# ------------------------------------------------------------------ neighbours
function check_neighbors!(rep, grid, m, S)
    F = :neighbors
    sh = cellshape(grid)
    provided = !(:neighbors in rep.fallbacks)
    if !(sh isa GrI.Polygonal) && !provided
        return skip!(rep, F, "non-polygonal grid without its own neighbors")
    end
    if GrI.neighbors_are_analytic(grid) && !provided
        law!(rep, F, "declared analytic neighbours are provided (not the fallback)", false,
            "neighbors_are_analytic is true but neighbors resolves to the generic fallback")
    end
    SS = S[1:min(end, 24)]
    N = Dict{Tuple{Int,Symbol},Vector{Int}}()
    nb(i, conn) = get!(() -> collect(neighbors(grid, i; connectivity=conn)), N, (i, nameof(typeof(conn))))
    guarded!(rep, F, "symmetry (Vertex and Edge)") do
        bad = 0
        # (a polar cell may touch tens of thousands of cells: check its first 64)
        for i in SS, conn in (Vertex(), Edge()), j in first(nb(i, conn), 64)
            i in nb(j, conn) || (bad += 1)
        end
        bad == 0, "$bad asymmetric pairs"
    end
    guarded!(rep, F, "Edge neighbours ⊆ Vertex neighbours") do
        all(i -> issubset(nb(i, Edge()), nb(i, Vertex())), SS)
    end
    guarded!(rep, F, "no self, no duplicates") do
        all(i -> !(i in nb(i, Vertex())) && allunique(nb(i, Vertex())), SS)
    end
    for conn in (Vertex(), Edge())
        mx = maxneighbors(grid, conn)
        if mx === nothing
            info!(rep, F, "maxneighbors($(nameof(typeof(conn))))", "nothing (unbounded/undeclared)")
        else
            guarded!(rep, F, "maxneighbors($(nameof(typeof(conn)))) = $mx never exceeded") do
                worst = maximum(i -> length(nb(i, conn)), SS)
                worst <= mx, "max seen $worst"
            end
        end
        w = winding(grid, conn)
        if w isa CounterClockwise || w isa Clockwise
            guarded!(rep, F, "declared winding $(nameof(typeof(w))) ($(nameof(typeof(conn))))") do
                bad = Int[]
                for i in SS
                    js = nb(i, conn)
                    c = cell_centroid(grid, i)
                    az = [_azimuth(m, c, cell_centroid(grid, j)) for j in js]
                    w isa Clockwise && reverse!(az)
                    is_ccw_sequence(az) || push!(bad, i)
                end
                isempty(bad), "cells $(first(bad, 5))"
            end
        else
            info!(rep, F, "winding($(nameof(typeof(conn))))", string(nameof(typeof(w))))
        end
    end
    guarded!(rep, F, "disc(i, k) is rings 1:k concatenated; maxneighbors(grid, k) bounds it") do
        ok = true; detail = ""
        for i in filter(i -> length(nb(i, Vertex())) <= 64, SS)[1:min(end, 8)], k in 1:2
            d = collect(neighbors(grid, i, k))
            ok &= allunique(d) && !(i in d)
            ok &= Set(d) == Set(vcat((collect(ring(grid, i, j)) for j in 1:k)...))
            mx = GrI.maxneighbors(grid, k, Vertex())
            if mx !== nothing && length(d) > mx
                ok = false; detail = "disc $k of cell $i has $(length(d)) > maxneighbors $mx"
            end
        end
        ok, detail
    end
    if cellshape(grid) isa SingleRing && conformity(grid) isa Conforming && winding(grid, Edge()) isa CounterClockwise
        guarded!(rep, F, "alignment: Edge neighbour k lies across edge (corner k, corner k+1)") do
            # each Edge neighbour holds both ends of one edge of i (corners are
            # shared exactly); edge indices increase, and equal 1:n on a full cell
            bad = Int[]
            for i in SS
                js = nb(i, Edge()); cs = collect(GrI.cell_corners(grid, i)); n = length(cs)
                es = [something(findfirst(e -> GrI.cell_contains(grid, j, cs[e]) && GrI.cell_contains(grid, j, cs[mod1(e + 1, n)]), 1:n), 0)
                      for j in js]
                ok = all(>(0), es) && issorted(es) && allunique(es) && (length(js) != n || es == 1:n)
                ok || push!(bad, i)
            end
            isempty(bad), "$(length(bad)) of $(length(SS)) cells, e.g. $(first(bad, 5))"
        end
    end
    if Structured.isstructured(grid) && sh isa SingleRing
        check_seams!(rep, grid, S)
    end
    guarded!(rep, F, "rings partition the disc (k ≤ 2)") do
        ok = true
        for i in filter(i -> length(nb(i, Vertex())) <= 64, SS)[1:min(end, 8)]
            r0 = ring(grid, i, 0); r1 = ring(grid, i, 1); r2 = ring(grid, i, 2)
            r0 = collect(r0); r1 = collect(r1); r2 = collect(r2)
            ok &= r0 == [i]
            ok &= Set(r1) == Set(nb(i, Vertex()))
            ok &= isempty(intersect(r1, r2)) && !(i in r2) && allunique(r2)
            ok &= all(j -> any(k -> k in r1, neighbors(grid, j)), r2)
        end
        ok
    end
    metric!(rep, "neighbors answered by", provided ? "provider" : "geometric fallback")
end

# --- structured seams: the arithmetic hooks against a shared-vertex scan ---------
# The oracle knows nothing of `facecross`: a cell touches another when they
# hold a common lattice vertex point, and lies across a side when it holds
# both of the side's points. Candidates are the in-face block around the cell
# and the border cells of every face (the only cells a seam can lead to).
struct SeamOracle{T}
    border::Vector{Tuple{Int,T}}       # (position, lattice corners) of every face's border cells
end
function SeamOracle(grid)
    c1 = Structured._lattice_corners(grid, 1, 1, 1)
    border = Tuple{Int,typeof(c1)}[]
    for f in 1:Structured.nfaces(grid)
        ni, nj = Structured.facesize(grid, f)
        for i in 1:ni, j in (1, nj)
            push!(border, (Structured.faceposition(grid, f, i, j), Structured._lattice_corners(grid, f, i, j)))
        end
        for j in 2:nj-1, i in (1, ni)
            push!(border, (Structured.faceposition(grid, f, i, j), Structured._lattice_corners(grid, f, i, j)))
        end
    end
    return SeamOracle(border)
end
# every candidate (position, corners) for cell (f, i, j): its in-face block, then the borders
function _candidates(o::SeamOracle, grid, f, i, j)
    ni, nj = Structured.facesize(grid, f)
    block = [(Structured.faceposition(grid, f, i2, j2), Structured._lattice_corners(grid, f, i2, j2))
             for i2 in max(1, i - 1):min(ni, i + 1), j2 in max(1, j - 1):min(nj, j + 1)]
    return Iterators.flatten((vec(block), o.border))
end
_holds(c, p) = any(q -> GrI._near(q, p), c)
# the cell's sides as (start, end), skipping collapsed ones
_sides(c) = ((c[s], c[mod1(s + 1, 4)]) for s in 1:4 if !GrI._near(c[s], c[mod1(s + 1, 4)]))
_hasside(c, a, b) = any(((x, y),) -> (GrI._near(x, a) && GrI._near(y, b)) || (GrI._near(x, b) && GrI._near(y, a)), _sides(c))
function _oracle_edge(o, grid, k, f, i, j, c)
    out = Set{Int}()
    for (a, b) in _sides(c), (q, cq) in _candidates(o, grid, f, i, j)
        q != k && _hasside(cq, a, b) && push!(out, q)
    end
    return out
end
function _oracle_vertex(o, grid, k, f, i, j, c)
    out = Set{Int}()
    for (q, cq) in _candidates(o, grid, f, i, j)
        q != k && any(p -> _holds(cq, p), c) && push!(out, q)
    end
    return out
end
_oracle_fan(o, grid, f, i, j, p) = Set{Int}(q for (q, cq) in _candidates(o, grid, f, i, j) if _holds(cq, p))

function check_seams!(rep, grid, S)
    F = :neighbors
    o = try
        SeamOracle(grid)
    catch e
        return law!(rep, F, "structured seams: shared-vertex oracle", false, "threw: " * first(sprint(showerror, e), 200))
    end
    # every cell on a small grid, else the sampled cells plus border cells spread along the borders
    cells = ncells(grid) <= 2_000 ? collect(1:ncells(grid)) :
        unique!(vcat(S[1:min(end, 16)], [o.border[k][1] for k in round.(Int, range(1, length(o.border); length=32))]))
    guarded!(rep, F, "structured seams: Edge and Vertex neighbours equal the shared-vertex oracle") do
        bad = Int[]
        for k in cells
            f, i, j = Structured.facecell(grid, k)
            c = Structured._lattice_corners(grid, f, i, j)
            Set(_nbrs(grid, k, Edge())) == _oracle_edge(o, grid, k, f, i, j, c) || push!(bad, k)
            Set(_nbrs(grid, k, Vertex())) == _oracle_vertex(o, grid, k, f, i, j, c) || push!(bad, k)
        end
        isempty(bad), "$(length(bad)) of $(length(cells)) cells, e.g. $(first(bad, 5))"
    end
    guarded!(rep, F, "structured seams: facecross lands on the side holding both points of the side, and crosses back") do
        bad = Int[]; seams = 0
        for k in cells
            f, i, j = Structured.facecell(grid, k)
            ni, nj = Structured.facesize(grid, f)
            c = Structured._lattice_corners(grid, f, i, j)
            for s in 1:4
                ((s == 1 && j == 1) || (s == 2 && i == ni) || (s == 3 && j == nj) || (s == 4 && i == 1)) || continue
                x = Structured.facecross(grid, f, i, j, s)
                x === nothing && continue
                seams += 1
                f2, i2, j2, s2 = x
                cq = Structured._lattice_corners(grid, f2, i2, j2)
                # the landing side runs the seam the other way: corner s is corner s2+1 across, corner s+1 is corner s2
                ok = 1 <= s2 <= 4 && _hasside(cq, c[s], c[mod1(s + 1, 4)]) &&
                     GrI._near(cq[mod1(s2 + 1, 4)], c[s]) && GrI._near(cq[s2], c[mod1(s + 1, 4)]) &&
                     Structured.facecross(grid, f2, i2, j2, s2) == (f, i, j, s)
                ok || push!(bad, k)
            end
        end
        metric!(rep, "seam crossings checked", seams)
        isempty(bad), "$(length(bad)) cells, e.g. $(first(bad, 5))"
    end
    guarded!(rep, F, "structured seams: vertex_cells at every corner equals the shared-vertex oracle") do
        bad = Int[]
        for k in cells
            f, i, j = Structured.facecell(grid, k)
            cs = GrI.cell_corners(grid, k)
            for t in 1:length(cs)
                key = GrI.vertex_key(grid, k, t)
                fan = collect(GrI.vertex_cells(grid, key))
                ok = Set(fan) == _oracle_fan(o, grid, f, i, j, cs[t]) && allunique(fan) &&
                     key.cell == minimum(fan) && GrI.vertex_key(grid, key.cell, key.corner) == key
                ok || push!(bad, k)
            end
        end
        isempty(bad), "$(length(bad)) cells, e.g. $(first(bad, 5))"
    end
end
_nbrs(grid, k, conn) = collect(neighbors(grid, k; connectivity=conn))

# ------------------------------------------------------------------ dual
function check_dual!(rep, grid, m, S)
    F = :dual
    if !(conformity(grid) isa Conforming)
        guarded!(rep, F, "dual refused on a non-conforming grid") do
            e = try GrI.dual(grid); nothing catch e; e end
            e isa ArgumentError, "got $(typeof(e))"
        end
        return skip!(rep, F, "grid is $(nameof(typeof(conformity(grid)))) / $(nameof(typeof(cellshape(grid))))")
    end
    if ncells(grid) > 50_000 && :nvertices in rep.fallbacks
        return skip!(rep, F, "geometric vertex merge skipped above 50k cells (ncells = $(ncells(grid)))")
    end
    if ncells(grid) > 2_000_000 && !GrI.vertex_numbering_is_analytic(grid)
        # the dense numbering is derived in O(n): check the keyed layer only
        GrI.vertices_are_analytic(grid) || return skip!(rep, F, "ncells = $(ncells(grid)) > 2M")
        guarded!(rep, F, "keyed vertices: each cell is in the fan of each of its corner keys") do
            all(i -> all(k -> i in GrI.vertex_cells(grid, GrI.vertex_key(grid, i, k)), 1:length(GrI.cell_corners(grid, i))), S)
        end
        return skip!(rep, F, "dense-numbering laws skipped: ncells = $(ncells(grid)) > 2M")
    end
    # The dense verbs are O(ncells) per call on a plain grid: the suite holds the
    # precomputation once, as a caller would, and checks the grid's own verbs
    # against it on a few samples.
    V = try
        GrI.GridAndVertices(grid)
    catch e
        return law!(rep, F, "vertex precomputation (GridAndVertices)", false, "threw: " * first(sprint(showerror, e), 200))
    end
    nv = try
        GrI.nvertices(V)
    catch e
        return law!(rep, F, "nvertices", false, "threw: " * first(sprint(showerror, e), 200))
    end
    metric!(rep, "nvertices", nv)
    metric!(rep, "vertex identity answered by", :nvertices in rep.fallbacks ? "geometric fallback" : "provider")
    if GrI.vertex_numbering_is_analytic(grid)
        guarded!(rep, F, "analytic numbering: no vertex precomputation is built") do
            V.table === nothing && GrI.dual(grid).table === nothing
        end
        guarded!(rep, F, "analytic numbering: vertex_key and vertex_index are inverse on the plain grid") do
            all(v -> GrI.vertex_index(grid, GrI.vertex_key(grid, v)) == v, sample_positions(nv, length(S))) &&
            all(i -> all(k -> GrI.vertex_key(grid, GrI.vertex_index(grid, k)) == k, GrI.cell_vertex_keys(grid, i)), S)
        end
        if !GrI.vertices_are_analytic(grid)
            guarded!(rep, F, "analytic numbering: the derived keyed verbs agree with the provider's dense verbs") do
                all(sample_positions(nv, length(S))) do v
                    k = GrI.vertex_key(grid, v)
                    collect(GrI.vertex_cells(grid, k)) == collect(GrI.vertex_cells(grid, v)) &&
                    _close(m, GrI.vertex_point(grid, k), GrI.vertex_point(grid, v)) &&
                    GrI.vertex_isopen(grid, k) == GrI.vertex_isopen(grid, v)
                end
            end
        end
    end
    guarded!(rep, F, "the grid's own dense vertex verbs agree with the held table (sampled)") do
        GrI.nvertices(grid) == nv &&
        all(i -> collect(GrI.cell_vertices(grid, i)) == collect(GrI.cell_vertices(V, i)), S[1:min(end, 3)]) &&
        all(v -> collect(GrI.vertex_cells(grid, v)) == collect(GrI.vertex_cells(V, v)) &&
                 GrI.vertex_point(grid, v) == GrI.vertex_point(V, v), sample_positions(nv, 3))
    end
    grid = V
    guarded!(rep, F, "cell_vertices and vertex_cells are inverse incidences") do
        ok = true
        for i in S, v in GrI.cell_vertices(grid, i)
            ok &= i in GrI.vertex_cells(grid, v)
        end
        for v in sample_positions(nv, length(S)), c in GrI.vertex_cells(grid, v)
            ok &= v in GrI.cell_vertices(grid, c)
        end
        ok
    end
    guarded!(rep, F, "cell_vertices are a subsequence of cell_corners") do
        bad = [i for i in S if !is_cyclic_subsequence(m, [GrI.vertex_point(grid, v) for v in GrI.cell_vertices(grid, i)],
                                                      collect(GrI.cell_corners(grid, i)))]
        isempty(bad), "cells $(first(bad, 5))"
    end
    guarded!(rep, F, "vertex keys round-trip; dense numbering derived from them") do
        all(v -> GrI.vertex_index(grid, GrI.vertex_key(grid, v)) == v, sample_positions(nv, length(S))) &&
        all(i -> [GrI.vertex_index(grid, k) for k in GrI.cell_vertex_keys(grid, i)] == collect(GrI.cell_vertices(grid, i)), S)
    end
    guarded!(rep, F, "vertex_isopen: open iff the fan does not close; none on a global grid") do
        Vs = sample_positions(nv, 4length(S))
        nopen = count(v -> GrI.vertex_isopen(grid, v), Vs)
        metric!(rep, "open vertices (sampled)", "$nopen / $(length(Vs))")
        GrI.isglobal(grid) ? (nopen == 0, "$nopen open vertices on a global grid") : true
    end
    guarded!(rep, F, "dual ring point k corresponds to fan cell k") do
        D = GrI.dual(grid); bad = 0
        for v in sample_positions(nv, length(S))
            GrI.vertex_isopen(grid, v) && continue
            cs = collect(GrI.vertex_cells(grid, v)); b = collect(cell_boundary(D, v))
            length(cs) == length(b) || (bad += 1; continue)
            all(k -> _close(m, b[k], cell_centroid(grid, cs[k]); tol=1e-6) || GrI.cell_contains(grid, cs[k], b[k]), eachindex(cs)) || (bad += 1)
        end
        bad == 0, "$bad vertices"
    end
    if m isa Spherical && GrI.isglobal(grid) && nv <= 200_000
        guarded!(rep, F, "the dual of a global grid tiles the sphere (areas sum to 4πr²)") do
            D = GrI.dual(grid)
            A = sum(v -> cell_area(D, v), 1:nv)
            want = 4π * Float64(m.radius)^2
            metric!(rep, "dual total area / 4πr²", A / want)
            isapprox(A, want; rtol=1e-6), "sum = $(A / want) × 4πr²"
        end
    end
    D = GrI.dual(grid)
    guarded!(rep, F, "dual cells are valid CCW rings (≥ 3 points)") do
        bad = Int[]; degenerate = 0
        for v in sample_positions(nv, length(S))
            b = cell_boundary(D, v)
            if length(b) < 3
                degenerate += 1; push!(bad, v); continue
            end
            # open fans put the vertex itself on the ring, so wind about the ring's own mean
            ring_is_ccw(m, b, ccw_centre(m, b, GrI.vertex_point(grid, v))) || push!(bad, v)
        end
        metric!(rep, "dual degenerate cells (sampled)", degenerate)
        isempty(bad), "$(length(bad)) bad dual cells (degenerate: $degenerate), e.g. $(first(bad, 5))"
    end
end

# ------------------------------------------------------------------ subsets and chunks
function check_subsets!(rep, grid, m, S)
    F = :subsets
    rel = GrI.gridrelation(grid)
    if rel isa GrI.ParentGrid
        P = parent(grid)
        guarded!(rep, F, "globalindex round-trips through the parent") do
            all(i -> (c = cellindex(grid, i); g = globalindex(grid, c); g isa Int && cellindex(P, g) == c), S)
        end
        guarded!(rep, F, "geometry equals the parent's") do
            all(i -> (c = cellindex(grid, i); j = localindex(P, c); cell_centroid(grid, i) == cell_centroid(P, j) &&
                cell_extent(grid, i) == cell_extent(P, j)), S)
        end
    elseif rel isa GrI.ParentSystem
        guarded!(rep, F, "a multi-level set refuses globalindex") do
            e = try globalindex(grid, cellindex(grid, 1)); nothing catch e; e end
            e isa ArgumentError
        end
        law!(rep, F, "a multi-level set answers system", GrI.system(grid) !== nothing)
    end
    # the generic Subset over this grid
    guarded!(rep, F, "generic subset: own order, globalindex, forwarded geometry") do
        cells = S[1:2:end]
        s = GrI.subset(grid, cells)
        ok = ncells(s) == length(cells)
        for (k, p) in enumerate(cells)
            c = cellindex(s, k)
            ok &= c == cellindex(grid, p) && localindex(s, c) == k && cell_centroid(s, k) == cell_centroid(grid, p)
            rel isa GrI.ParentSystem || (ok &= globalindex(s, c) == globalindex(grid, c))
        end
        ok
    end
    guarded!(rep, F, "subset algebra agrees with set algebra on parent positions") do
        n = ncells(grid)
        a = GrI.subset(grid, S[1:2:end]); b = GrI.subset(grid, S[1:3:end])
        pa, pb = collect(GrI.parentpositions(a)), collect(GrI.parentpositions(b))
        collect(GrI.parentpositions(union(a, b))) == sort!(union(pa, pb)) &&
        collect(GrI.parentpositions(intersect(a, b))) == sort!(intersect(pa, pb)) &&
        collect(GrI.parentpositions(setdiff(a, b))) == sort!(setdiff(pa, pb))
    end
    if cellshape(grid) isa SingleRing || !(:neighbors in rep.fallbacks)
        guarded!(rep, F, "subset neighbours: clipped plain form, unclipped Index forms") do
            s = GrI.subset(grid, S[1:2:end]); ok = true
            for k in 1:min(ncells(s), 6)
                p = GrI.parentpositions(s)[k]
                full = collect(neighbors(grid, p))
                loc = collect(neighbors(s, k, GrI.Index(GrI.Local())))
                ok &= length(loc) == length(full)
                ok &= all(((j, l),) -> l == something(localindex(s, cellindex(grid, j)), GrI.OUTSIDE), zip(full, loc))
                ok &= collect(neighbors(s, k)) == filter(!=(GrI.OUTSIDE), loc)
                rel isa GrI.NoParent && (ok &= collect(neighbors(s, k, GrI.Index(GrI.Global()))) == full)
            end
            ok
        end
    end
    ch = GrI.chunks(grid)
    if ch === nothing
        info!(rep, F, "chunks", "none declared")
    else
        metric!(rep, "nchunks", length(ch))
        if ncells(grid) <= 10_000_000
            guarded!(rep, F, "chunks partition the grid") do
                seen = zeros(Int, ncells(grid))
                for c in ch, k in 1:ncells(c)
                    seen[localindex(grid, cellindex(c, k))] += 1
                end
                all(==(1), seen), "counts $(extrema(seen))"
            end
        end
        guarded!(rep, F, "chunkat agrees with chunk membership (sampled)") do
            bad = 0
            for j in sample_positions(length(ch), 16)
                c = ch[j]
                for k in sample_positions(ncells(c), 8)
                    GrI.chunkat(grid, localindex(grid, cellindex(c, k))) == j || (bad += 1)
                end
            end
            bad == 0, "$bad cells"
        end
        guarded!(rep, F, "chunk extents are cheap") do
            c = ch[1]; GrI.grid_extent(c)
            t = @elapsed for j in sample_positions(length(ch), 16); GrI.grid_extent(ch[j]); end
            metric!(rep, "chunk extent μs", round(t / 16 * 1e6; digits=2))
            t / 16 < 1e-3, "$(t/16) s per chunk extent"
        end
        guarded!(rep, F, "each chunk's extent covers its cells") do
            bad = 0
            for c in ch[sample_positions(length(ch), 16)]
                e = GrI.grid_extent(c)
                for k in sample_positions(ncells(c), 8)
                    all(p -> extent_contains_point(m, e, p), _cellpts(c, k)) || (bad += 1)
                end
            end
            bad == 0, "$bad cells outside their chunk extent"
        end
    end
end

# ------------------------------------------------------------------ performance (grid side)
function _sweep_centroid(grid, S)
    s = 0.0
    for i in S
        p = cell_centroid(grid, i)
        s += Float64(GI.x(p))
    end
    return s
end
function _sweep_extent(grid, S)
    s = 0
    for i in S
        e = cell_extent(grid, i)
        s += sizeof(e)
    end
    return s
end
function _sweep_identity(grid, S)
    s = 0
    for i in S
        s += localindex(grid, cellindex(grid, i))::Int
    end
    return s
end
_alloc(f, args...) = (f(args...); @allocated f(args...))

function check_grid_performance!(rep, grid, S)
    F = :performance
    for (name, f) in (("centroid", _sweep_centroid), ("extent", _sweep_extent), ("identity", _sweep_identity))
        guarded!(rep, F, "forwarding sweep over $name allocates nothing") do
            a = _alloc(f, grid, S)
            metric!(rep, "sweep $name bytes", a)
            a == 0, "$a bytes over $(length(S)) cells"
        end
    end
    guarded!(rep, F, "cell_centroid infers") do
        T = Base.return_types(cell_centroid, Tuple{typeof(grid),Int})
        length(T) == 1 && isconcretetype(only(T)), string(T)
    end
    try
        GrI.treeify(grid)
    catch e
        return law!(rep, F, "treeify", false, "threw: " * first(sprint(showerror, e), 200))
    end
    t = @elapsed GrI.treeify(grid)
    a = _alloc(GrI.treeify, grid)
    metric!(rep, "treeify time", t)
    metric!(rep, "treeify bytes", a)
    if GrI.location_is_analytic(grid)
        p = cell_centroid(grid, S[end])
        guarded!(rep, F, "analytic cellat allocates nothing") do
            a = _alloc(cellat, grid, p)
            metric!(rep, "cellat bytes", a)
            a == 0, "$a bytes"
        end
    else
        p = cell_centroid(grid, S[end])
        metric!(rep, "cellat (fallback) bytes", _alloc(cellat, grid, p))
    end
end
