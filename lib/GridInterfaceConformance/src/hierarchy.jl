"""
    check_system(sys; name, levels_checked = 3, depth = 2, nsamples = 24,
                 seed = DEFAULT_SEED, rng = nothing) -> ConformanceReport

Hierarchy laws on a system, plus the hierarchy performance laws, plus
`check_grid` on a small level grid and the multi-level-set laws on a small set
built from the system.
"""
function check_system(sys; name::AbstractString=string(nameof(typeof(sys))), levels_checked::Int=3,
        depth::Int=2, nsamples::Int=24, grid_level=nothing, seed::Integer=DEFAULT_SEED,
        rng::Union{Nothing,Random.AbstractRNG}=nothing, quiet::Bool=true)
    rng, seed = sampling_rng(seed, rng)
    rep = ConformanceReport(String(name); seed)
    F = :hierarchy
    Ls = collect(GrI.levels(sys))
    l0 = first(Ls)
    metric!(rep, "levels", "$(first(Ls)):$(last(Ls))")
    guarded!(rep, F, "levels ascending") do
        issorted(Ls) && allunique(Ls)
    end
    guarded!(rep, F, "rootcells equal the coarsest level, ascending") do
        rc = collect(GrI.rootcells(sys))
        g = GrI.levelgrid(sys, l0)
        rc == [cellindex(g, i) for i in 1:ncells(g)] && issorted(rc)
    end
    guarded!(rep, F, "level(c) on ids of every checked level") do
        all(l -> (g = GrI.levelgrid(sys, l); all(i -> GrI.level(cellindex(g, i)) == l, sample_positions(rng, ncells(g), 8))),
            Ls[1:min(end, levels_checked + 1)])
    end
    guarded!(rep, F, "level grids answer system and level") do
        all(l -> (g = GrI.levelgrid(sys, l); GrI.system(g) == sys && GrI.level(g) == l), Ls[1:min(end, levels_checked + 1)])
    end
    guarded!(rep, F, "ids are AbstractCellIndex (or declared through iscellindex)") do
        c = cellindex(GrI.levelgrid(sys, l0), 1)
        GrI.iscellindex(c), "$(typeof(c)) is neither an AbstractCellIndex nor declared with iscellindex"
    end
    guarded!(rep, F, "Extents.union of node extents is defined and covers both") do
        rc = collect(GrI.rootcells(sys))
        length(rc) < 2 && return true
        a, b = GrI.node_extent(sys, rc[1]), GrI.node_extent(sys, rc[end])
        u = Extents.union(a, b)
        m = best_manifold(GrI.levelgrid(sys, l0))
        extent_contains_extent(m, u, a) && extent_contains_extent(m, u, b), "union $(u)"
    end
    GrI.has_sorted_subtrees(sys) && guarded!(rep, F, "chunks(levelgrid, l) partitions by ancestor (sampled)") do
        l = Ls[min(end, 2)]
        lg = GrI.levelgrid(sys, l)
        ch = GrI.chunks(lg, l0)
        ok = length(ch) == ncells(GrI.levelgrid(sys, l0))
        for i in sample_positions(rng, ncells(lg), 32)
            j = GrI.chunkat(ch, i)
            ok &= i in GrI.parentpositions(ch[j])
        end
        ok
    end
    # bushy walk: roots and their descendants for `depth` levels
    walk = _bushy(sys, depth, nsamples, rng)
    guarded!(rep, F, "parent/children inverses; children ascending, one level down") do
        ok = true
        for c in walk[sample_positions(rng, length(walk), 4nsamples)]
            GrI.level(c) >= last(Ls) && continue
            ch = collect(GrI.children(sys, c))
            ok &= issorted(ch) && !isempty(ch)
            ok &= all(d -> Base.parent(sys, d) == c && GrI.level(d) == GrI.level(c) + 1, ch)
        end
        ok
    end
    guarded!(rep, F, "ancestor agrees with the parent loop") do
        all(c -> GrI.ancestor(sys, c, l0) == _rootof(sys, c), walk)
    end
    # the children of all level-l cells partition level l+1, and every level-(l+1)
    # cell has exactly one parent: exhaustive on small levels, sampled above
    guarded!(rep, F, "children of every level-l cell partition level l+1 (exhaustive below 50k cells, sampled above)") do
        ok = true; detail = ""
        for l in Ls[1:min(end - 1, levels_checked)]
            lg, fg = GrI.levelgrid(sys, l), GrI.levelgrid(sys, l + 1)
            if ncells(fg) <= 50_000
                seen = zeros(Int, ncells(fg))
                for i in 1:ncells(lg), d in GrI.children(sys, cellindex(lg, i))
                    j = localindex(fg, d)
                    j === nothing && (ok = false; detail = "child $d of $(cellindex(lg, i)) is not in level $(l + 1)"; break)
                    seen[j] += 1
                end
                ok && extrema(seen) != (1, 1) && (ok = false; detail = "level $(l + 1) is covered $(extrema(seen)) times by the children of level $l")
            else
                for j in sample_positions(rng, ncells(fg), 32)
                    d = cellindex(fg, j)
                    count(==(d), GrI.children(sys, Base.parent(sys, d))) == 1 || (ok = false; detail = "$d is not exactly once among its parent's children")
                end
            end
            ok || break
        end
        ok, detail
    end
    if GrI.has_sorted_subtrees(sys)
        guarded!(rep, F, "descendant_range is a UnitRange{Int} exact against a children-recursion oracle") do
            ok = true; detail = ""
            for c in walk[sample_positions(rng, length(walk), 12)], dl in 0:2
                l = GrI.level(c) + dl
                l > last(Ls) && continue
                r = GrI.descendant_range(sys, c, l)
                r isa UnitRange{Int} || (ok = false; detail = "descendant_range($c, $l) is a $(typeof(r))"; break)
                lg = GrI.levelgrid(sys, l)
                actual = _descend(sys, c, l)
                # a range wildly wider than the descendant set is itself the
                # violation; materialising it to say so would exhaust memory
                length(r) > 8 * length(actual) + 64 && (ok = false; detail = "descendant_range($c, $l) spans $(length(r)) positions for $(length(actual)) descendants"; break)
                pos = sort!([localindex(lg, d) for d in actual])
                if pos != collect(r)
                    ok = false; detail = "descendant_range($c, $l) = $r; the oracle gives $(length(pos)) positions $(isempty(pos) ? "" : extrema(pos))"
                end
            end
            ok, detail
        end
        guarded!(rep, F, "sibling ranges partition the parent's range") do
            ok = true
            for c in walk[1:min(end, 12)]
                l = GrI.level(c) + 2
                l > last(Ls) && continue
                rs = [GrI.descendant_range(sys, d, l) for d in GrI.children(sys, c)]
                ok &= reduce(vcat, collect.(rs)) == collect(GrI.descendant_range(sys, c, l))
            end
            ok
        end
    else
        skip!(rep, F, "descendant_range: system does not declare sorted subtrees")
    end
    guarded!(rep, F, "covering law on a bushy walk (node_extent ⊇ descendant geometry)") do
        bad = 0; m = nothing
        for c in walk[sample_positions(rng, length(walk), 12)]
            e = GrI.node_extent(sys, c)
            for dl in 0:2
                l = GrI.level(c) + dl
                l > last(Ls) && continue
                lg = GrI.levelgrid(sys, l)
                m = best_manifold(lg)
                for d in _descend(sys, c, l)[sample_positions(rng, length(_descend(sys, c, l)), 16)]
                    i = localindex(lg, d)
                    all(p -> extent_contains_point(m, e, p), _cellpts(lg, i)) || (bad += 1)
                end
            end
        end
        bad == 0, "$bad descendant cells outside their ancestor's node_extent"
    end
    # depth alone is not the law: from sampled roots, one random chain of
    # children to the deepest level, every cell inside every ancestor's extent
    guarded!(rep, F, "covering law on a chain to the deepest level (every ancestor's node_extent)") do
        bad = 0; cells = 0
        rc = collect(GrI.rootcells(sys))
        for r in rc[sample_positions(rng, length(rc), 6)]
            chain = [(r, GrI.node_extent(sys, r))]
            c = r
            while GrI.level(c) < last(Ls)
                ch = collect(GrI.children(sys, c))
                isempty(ch) && break
                c = ch[rand(rng, 1:length(ch))]
                push!(chain, (c, GrI.node_extent(sys, c)))
                lg = GrI.levelgrid(sys, GrI.level(c))
                m = best_manifold(lg)
                pts = _cellpts(lg, localindex(lg, c))
                cells += 1
                for (_, e) in chain
                    all(p -> extent_contains_point(m, e, p), pts) || (bad += 1; break)
                end
            end
        end
        metric!(rep, "covering chain cells", cells)
        bad == 0, "$bad cells escape an ancestor's node_extent"
    end
    check_system_performance!(rep, sys, walk)
    # the level grid itself, and a multi-level set
    gl = grid_level === nothing ? Ls[min(end, 3)] : grid_level
    lg = GrI.levelgrid(sys, gl)
    sub = check_grid(lg; name="$(name) level $gl", nsamples, rng)
    _merge!(rep, sub, "level $gl")
    if gl + 1 <= last(Ls)
        _merge!(rep, check_laziness(lg, GrI.levelgrid(sys, gl + 1); name="$(name) levels $gl / $(gl + 1)"), "laziness")
    end
    mls = _sample_multilevel(sys, gl)
    if mls !== nothing
        _merge!(rep, check_multilevel(mls; name="$(name) multi-level set", rng), "multi-level")
    end
    quiet || print_report(rep)
    return rep
end

"""
    test_system(sys; levels = 2:5, kw...) -> ConformanceReport

[`check_system`](@ref), emitted as a `@testset` with one `@test` per law and a
`@test_skip` per suite the system does not claim. `levels` are the level
numbers the grid laws are run on (those the system has, up to three of them;
the hierarchy laws always walk from the roots); the remaining keywords are
[`check_system`](@ref)'s.
"""
function test_system(sys; levels=2:5, kw...)
    ls = [l for l in levels if l in GrI.levels(sys)]
    isempty(ls) && (ls = [GrI.levels(sys)[min(end, 3)]])
    rep = check_system(sys; grid_level=first(ls), kw...)
    for l in ls[2:min(end, 3)]
        _merge!(rep, check_grid(GrI.levelgrid(sys, l); name="$(rep.name) level $l", nsamples=24, seed=something(rep.seed, DEFAULT_SEED)), "level $l")
    end
    emit_tests(rep)
    return rep
end

function _merge!(rep, sub, tag)
    for r in sub.results
        push!(rep.results, LawResult(r.family, "[$tag] " * r.law, r.status, r.detail))
    end
    for (k, v) in sub.metrics
        rep.metrics["[$tag] " * k] = v
    end
    for f in sub.fallbacks
        push!(rep.fallbacks, Symbol("[$tag] ", f))
    end
end

function _bushy(sys, depth, nsamples, rng)
    out = collect(GrI.rootcells(sys))
    # expand only a sample of each frontier: a system may have 10^4 roots with
    # 10^5 children each (CopernicusDEM), so the full walk is not affordable
    frontier = out[sample_positions(rng, length(out), nsamples)]
    for _ in 1:depth
        next = eltype(out)[]
        for c in frontier
            GrI.level(c) >= last(GrI.levels(sys)) && continue
            append!(next, GrI.children(sys, c))
        end
        next = next[sample_positions(rng, length(next), nsamples)]
        append!(out, next)
        frontier = next
    end
    return out
end
_rootof(sys, c) = GrI.level(c) == first(GrI.levels(sys)) ? c : _rootof(sys, Base.parent(sys, c))
function _descend(sys, c, l)
    cs = [c]
    while GrI.level(first(cs)) < l
        cs = reduce(vcat, [collect(GrI.children(sys, x)) for x in cs])
    end
    return cs
end

# a small mixed-level set: one root refined two levels, the other roots as is
function _sample_multilevel(sys, gl)
    Ls = GrI.levels(sys)
    l0 = first(Ls)
    l0 + 2 > last(Ls) && return nothing
    rc = collect(GrI.rootcells(sys))
    cells = eltype(rc)[]
    first_children = collect(GrI.children(sys, rc[1]))
    append!(cells, collect(GrI.children(sys, first_children[1])))
    append!(cells, first_children[2:end])
    append!(cells, rc[2:end])
    return GrI.MultiLevelSet(sys, cells)
end

"""
    check_multilevel(set; name) -> ConformanceReport

The multi-level-set laws (antichain, tessellation trait = refinement trait,
refused globalindex and dual, tree covers every member once) plus `check_grid`.
"""
function check_multilevel(s; name::AbstractString="multi-level set", nsamples::Int=48,
        seed::Integer=DEFAULT_SEED, rng::Union{Nothing,Random.AbstractRNG}=nothing)
    rep = check_grid(s; name, nsamples, seed, rng)
    F = :multilevel
    sys = GrI.system(s)
    law!(rep, F, "no member is an ancestor of another", GrI.check_antichain(s))
    guarded!(rep, F, "tessellation trait equals the system's refinement trait") do
        (tessellation(s) isa Tessellating) == GrI.has_congruent_refinement(sys)
    end
    guarded!(rep, F, "dual refused (hanging vertices)") do
        e = try GrI.dual(s); nothing catch e; e end
        e isa ArgumentError
    end
    guarded!(rep, F, "levelgrid infers to one concrete type (forwarding is type stable)") do
        T = Base.return_types(GrI.levelgrid, Tuple{typeof(sys),Int})
        length(T) == 1 && isconcretetype(only(T)), string(T)
    end
    return rep
end

# ------------------------------------------------------------------ performance (system side)
function check_system_performance!(rep, sys, walk)
    F = :performance
    guarded!(rep, F, "levelgrid infers to one concrete type over all levels") do
        T = Base.return_types(GrI.levelgrid, Tuple{typeof(sys),Int})
        types = unique(typeof(GrI.levelgrid(sys, l)) for l in GrI.levels(sys))
        ok = length(T) == 1 && isconcretetype(only(T)) && length(types) == 1
        ok, "return_types $(T); runtime types $(length(types))"
    end
    guarded!(rep, F, "hierarchy verbs infer to one concrete type (children, parent, rootcells, node_extent, descendant_range)") do
        c = walk[min(end, 2)]
        problems = String[]
        checks = Any[("children(sys, c)", _infer_children, (typeof(sys), typeof(c))), ("parent(sys, c)", _infer_parent, (typeof(sys), typeof(c))),
            ("rootcells(sys)", _infer_rootcells, (typeof(sys),)), ("node_extent(sys, c)", _infer_extent_sys, (typeof(sys), typeof(c)))]
        GrI.has_sorted_subtrees(sys) && push!(checks, ("descendant_range(sys, c, l)", _infer_drange, (typeof(sys), typeof(c), Int)))
        for (label, f, A) in checks
            types = Base.return_types(f, Tuple{A...})
            length(types) == 1 && _inferred_ok(only(types)) || push!(problems, "$label infers $types")
        end
        isempty(problems), join(problems, "; ")
    end
    c = walk[min(end, 2)]
    GrI.children(sys, c); GrI.rootcells(sys)
    ac = @allocated GrI.children(sys, c)
    ar = @allocated GrI.rootcells(sys)
    metric!(rep, "children bytes", ac)
    metric!(rep, "rootcells bytes", ar)
    info!(rep, F, "children/rootcells allocation (reported, not failed)", "children $ac B, rootcells $ar B")
end

# ------------------------------------------------------------------ laziness
"""
    check_laziness(small, large; name) -> ConformanceReport

The laziness law on two sizes of the same grid type (two levels of a system,
two resolutions of a raster): no verb's bytes scale with `ncells` between
them, and no verb retains state across sizes. Each verb is called on `small`,
then `large`, then `small` again: the two small answers are equal, the second
small call allocates no more than the first, and the large call allocates no
more than a constant over the small one. Point location on a grid without an
analytic `cellat` is measured through `GridAndTree` (the caller's tree);
geometric neighbours are reported, not measured (they build a tree per call).
"""
function check_laziness(small, large; name::AbstractString="$(nameof(typeof(small))) × 2", quiet::Bool=true)
    rep = ConformanceReport(String(name))
    check_laziness!(rep, small, large)
    quiet || print_report(rep)
    return rep
end

_treekey(t) = (typeof(t), GrI.split_weight(t))
# a middle position: the last one is a polar cell on lon/lat grids, whose
# neighbourhood is legitimately proportional to the resolution
_mid(g) = (ncells(g) + 1) ÷ 2
_midpt(g) = cell_centroid(g, _mid(g))

function check_laziness!(rep, small, large)
    F = :performance
    typeof(small) === typeof(large) || return info!(rep, F, "laziness law", "needs two sizes of one grid type, got $(typeof(small)) and $(typeof(large))")
    ncells(small) < ncells(large) || return info!(rep, F, "laziness law", "the second grid is not larger")
    metric!(rep, "laziness sizes", "$(ncells(small)) / $(ncells(large))")
    verbs = Pair{String,Any}[
        "cellindex" => g -> cellindex(g, _mid(g)),
        "localindex ∘ cellindex" => g -> localindex(g, cellindex(g, _mid(g))),
        "cell_centroid" => g -> cell_centroid(g, _mid(g)),
        "cell_extent" => g -> cell_extent(g, _mid(g)),
        "cell_area" => g -> cell_area(g, _mid(g)),
        "cell_boundary" => g -> collect(_cellpts(g, _mid(g))),
        "chunks" => g -> (c = GrI.chunks(g); c === nothing ? nothing : length(c)),
        "subset" => g -> ncells(GrI.subset(g, [1, ncells(g)])),
    ]
    # The R-tree fallback is O(ncells) by contract, and so is any tree over a
    # system without sorted subtrees (the selection cursor, a provider's index
    # tree); a grid with its own O(1) tree or the window-mode cursor is measured.
    sys = GrI.system(small)
    fb = fallback_verbs(small)
    if :treeify in fb
        info!(rep, F, "laziness: treeify", "the R-tree fallback is O(ncells) to build: wrap once in GridAndTree; not measured")
    elseif sys !== nothing && !GrI.has_sorted_subtrees(sys)
        info!(rep, F, "laziness: treeify", "no sorted subtrees: a tree materialises the selection, O(ncells); not measured")
    else
        push!(verbs, "treeify" => g -> _treekey(treeify(g)),
            "treeify ForSweeps" => g -> _treekey(treeify(g; leaf_capacity=GrI.ForSweeps())),
            "treeify ForQueries" => g -> _treekey(treeify(g; leaf_capacity=GrI.ForQueries())),
            "grid_extent" => g -> GrI.grid_extent(g),
            "subtree" => g -> _treekey(subtree(g, [1, ncells(g)])))
    end
    if :neighbors in fb
        info!(rep, F, "laziness: neighbors", "geometric fallback builds a tree per call; not measured")
    else
        push!(verbs, "neighbors" => g -> collect(neighbors(g, _mid(g))))
    end
    if GrI.location_is_analytic(small)
        push!(verbs, "cellat" => g -> cellat(g, _midpt(g)))
    else
        # the caller's tree, built outside the measured call
        ws, wl = GrI.GridAndTree(small; leaf_capacity=GrI.ForQueries()), GrI.GridAndTree(large; leaf_capacity=GrI.ForQueries())
        push!(verbs, "cellat through GridAndTree" => g -> cellat(g === small ? ws : wl, _midpt(g)))
    end
    for (vname, f) in verbs
        guarded!(rep, F, "laziness: $vname") do
            f(small); f(large)                       # compile
            r1 = f(small); b1 = @allocated f(small)
            bL = @allocated f(large)
            r2 = f(small); b2 = @allocated f(small)
            metric!(rep, "laziness bytes $vname", "$b1 / $bL / $b2")
            ok = r1 == r2 && b2 <= b1 && bL <= 2b1 + 256
            ok, "small $b1 B, large $bL B, small again $b2 B" * (r1 == r2 ? "" : "; answers differ")
        end
    end
end

_infer_children(s, c) = GrI.children(s, c)
_infer_parent(s, c) = Base.parent(s, c)
_infer_rootcells(s) = GrI.rootcells(s)
_infer_extent_sys(s, c) = GrI.node_extent(s, c)
_infer_drange(s, c, l) = GrI.descendant_range(s, c, l)
