struct NotAScheme end   # an id scheme no grid declares

"""
    check_grid(grid; name, nsamples = 64, full_tree_limit = 200_000) -> ConformanceReport

Run every family's laws that `grid` claims. Never throws.
"""
function check_grid(grid; name::AbstractString=string(nameof(typeof(grid))), nsamples::Int=64,
        full_tree_limit::Int=200_000, area_rtol::Real=1e-3, quiet::Bool=true)
    rep = ConformanceReport(String(name))
    n = ncells(grid)
    m = best_manifold(grid)
    S = sample_positions(n, nsamples)
    append!(rep.fallbacks, fallback_verbs(grid))
    metric!(rep, "ncells", n)
    metric!(rep, "manifold", string(nameof(typeof(m))))
    metric!(rep, "cellshape", string(nameof(typeof(cellshape(grid)))))
    check_identity!(rep, grid, n, S)
    tree = check_tree!(rep, grid, n, m, S, full_tree_limit)
    check_geometry!(rep, grid, m, S, area_rtol)
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
    law!(rep, F, "every node's extent contains its children's extents (nested; PRD table)",
        nested_bad[] == 0, "$(nested_bad[]) non-nested parent/child extent pairs over $(visited[]) nodes")
    guarded!(rep, F, "subtree(grid, cells) emits parent positions") do
        cells = S[1:2:end]
        ls = sort!(collect_leaves(subtree(grid, cells)))
        ls == cells, "expected $(length(cells)) got $(length(ls))"
    end
    guarded!(rep, F, "GridAndTree forwards every grid verb") do
        g = GrI.GridAndTree(grid, tree)
        ok = ncells(g) == n && best_manifold(g) == m && cellshape(g) == cellshape(grid) &&
             all(i -> cellindex(g, i) == cellindex(grid, i) && cell_extent(g, i) == cell_extent(grid, i) &&
                          cell_centroid(g, i) == cell_centroid(grid, i), S)
        if cellshape(grid) isa SingleRing
            ok &= all(i -> collect(cell_boundary(g, i)) == collect(cell_boundary(grid, i)), S)
        end
        ok &= GrI.treeify(g) === tree
        ok
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

# ------------------------------------------------------------------ geometry
function check_geometry!(rep, grid, m, S, area_rtol)
    F = :geometry
    sh = cellshape(grid)
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
            bad = [i for i in S if !ring_is_ccw(m, collect(cell_boundary(grid, i)), cell_centroid(grid, i))]
            isempty(bad), "$(length(bad)) of $(length(S)) cells, e.g. $(first(bad, 5))"
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
    if !(:cell_area in rep.fallbacks)
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
    guarded!(rep, F, "cell_extent contains the cell") do
        bad = [i for i in S if !all(p -> extent_contains_point(m, cell_extent(grid, i), p), _cellpts(grid, i))]
        isempty(bad), "cells $(first(bad, 5))"
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
        a == geo && (agree_geo += 1)
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
    law!(rep, F, analytic ? "analytic cellat agrees with the geometric answer" : "cellat agrees with the geometric answer",
        agree_geo == agree_n, "$(agree_n - agree_geo) of $agree_n disagree")
end

# ------------------------------------------------------------------ neighbours
function check_neighbors!(rep, grid, m, S)
    F = :neighbors
    sh = cellshape(grid)
    provided = !(:neighbors in rep.fallbacks)
    if !(sh isa SingleRing) && !provided
        return skip!(rep, F, "GeneralCells grid without its own neighbors")
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
        for i in SS, conn in (Vertex(), Edge()), j in nb(i, conn)
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
    guarded!(rep, F, "rings partition the disc (k ≤ 2)") do
        ok = true
        for i in SS[1:min(end, 8)]
            r0 = ring(grid, i, 0); r1 = ring(grid, i, 1); r2 = ring(grid, i, 2)
            ok &= r0 == [i]
            ok &= Set(r1) == Set(nb(i, Vertex()))
            ok &= isempty(intersect(r1, r2)) && !(i in r2) && allunique(r2)
            ok &= all(j -> any(k -> k in r1, neighbors(grid, j)), r2)
        end
        ok
    end
    metric!(rep, "neighbors answered by", provided ? "provider" : "geometric fallback")
end

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
    nv = try
        GrI.nvertices(grid)
    catch e
        return law!(rep, F, "nvertices", false, "threw: " * first(sprint(showerror, e), 200))
    end
    metric!(rep, "nvertices", nv)
    metric!(rep, "vertex identity answered by", :nvertices in rep.fallbacks ? "geometric fallback" : "provider")
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
    guarded!(rep, F, "cell_vertices follow cell_corners") do
        all(i -> length(GrI.cell_vertices(grid, i)) == length(GrI.cell_corners(grid, i)), S)
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
            ring_is_ccw(m, b, m isa Spherical ? ring_mean(b) : GrI.vertex_point(grid, v)) || push!(bad, v)
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
    ch = GrI.chunks(grid)
    if ch === nothing
        info!(rep, F, "chunks", "none declared")
    else
        guarded!(rep, F, "chunks partition the grid") do
            seen = zeros(Int, ncells(grid))
            for c in ch, k in 1:ncells(c)
                seen[localindex(grid, cellindex(c, k))] += 1
            end
            all(==(1), seen), "counts $(extrema(seen))"
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
    t = @elapsed GrI.treeify(grid)
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
