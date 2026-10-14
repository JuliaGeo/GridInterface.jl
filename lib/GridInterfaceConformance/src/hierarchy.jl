"""
    check_system(sys; name, levels_checked = 3, depth = 2, nsamples = 24) -> ConformanceReport

Hierarchy laws on a system, plus the hierarchy performance laws, plus
`check_grid` on a small level grid and the multi-level-set laws on a small set
built from the system.
"""
function check_system(sys; name::AbstractString=string(nameof(typeof(sys))), levels_checked::Int=3,
        depth::Int=2, nsamples::Int=24, grid_level=nothing, quiet::Bool=true)
    rep = ConformanceReport(String(name))
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
        all(l -> (g = GrI.levelgrid(sys, l); all(i -> GrI.level(cellindex(g, i)) == l, sample_positions(ncells(g), 8))),
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
        for i in sample_positions(ncells(lg), 32)
            j = GrI.chunkat(ch, i)
            ok &= i in GrI.parentpositions(ch[j])
        end
        ok
    end
    # bushy walk: roots and their descendants for `depth` levels
    walk = _bushy(sys, depth, nsamples)
    guarded!(rep, F, "parent/children inverses; children ascending, one level down") do
        ok = true
        for c in walk[sample_positions(length(walk), 4nsamples)]
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
    if GrI.has_sorted_subtrees(sys)
        guarded!(rep, F, "descendant_range exact against a children-recursion oracle") do
            ok = true
            for c in walk[sample_positions(length(walk), 12)], dl in 0:2
                l = GrI.level(c) + dl
                l > last(Ls) && continue
                lg = GrI.levelgrid(sys, l)
                pos = sort!([localindex(lg, d) for d in _descend(sys, c, l)])
                ok &= pos == collect(GrI.descendant_range(sys, c, l))
            end
            ok
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
        for c in walk[sample_positions(length(walk), 12)]
            e = GrI.node_extent(sys, c)
            for dl in 0:2
                l = GrI.level(c) + dl
                l > last(Ls) && continue
                lg = GrI.levelgrid(sys, l)
                m = best_manifold(lg)
                for d in _descend(sys, c, l)[sample_positions(length(_descend(sys, c, l)), 16)]
                    i = localindex(lg, d)
                    all(p -> extent_contains_point(m, e, p), _cellpts(lg, i)) || (bad += 1)
                end
            end
        end
        bad == 0, "$bad descendant cells outside their ancestor's node_extent"
    end
    check_system_performance!(rep, sys, walk)
    # the level grid itself, and a multi-level set
    gl = grid_level === nothing ? Ls[min(end, 3)] : grid_level
    lg = GrI.levelgrid(sys, gl)
    sub = check_grid(lg; name="$(name) level $gl", nsamples=nsamples)
    _merge!(rep, sub, "level $gl")
    mls = _sample_multilevel(sys, gl)
    if mls !== nothing
        _merge!(rep, check_multilevel(mls; name="$(name) multi-level set"), "multi-level")
    end
    quiet || print_report(rep)
    return rep
end

"""
    test_system(sys; kw...) -> ConformanceReport

[`check_system`](@ref), emitted as a `@testset` with one `@test` per law.
"""
function test_system(sys; kw...)
    rep = check_system(sys; kw...)
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

function _bushy(sys, depth, nsamples)
    out = collect(GrI.rootcells(sys))
    # expand only a sample of each frontier: a system may have 10^4 roots with
    # 10^5 children each (CopernicusDEM), so the full walk is not affordable
    frontier = out[sample_positions(length(out), nsamples)]
    for _ in 1:depth
        next = eltype(out)[]
        for c in frontier
            GrI.level(c) >= last(GrI.levels(sys)) && continue
            append!(next, GrI.children(sys, c))
        end
        next = next[sample_positions(length(next), nsamples)]
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
function check_multilevel(s; name::AbstractString="multi-level set", nsamples::Int=48)
    rep = check_grid(s; name, nsamples)
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
    c = walk[min(end, 2)]
    GrI.children(sys, c); GrI.rootcells(sys)
    ac = @allocated GrI.children(sys, c)
    ar = @allocated GrI.rootcells(sys)
    metric!(rep, "children bytes", ac)
    metric!(rep, "rootcells bytes", ar)
    info!(rep, F, "children/rootcells allocation (reported, not failed)", "children $ac B, rootcells $ar B")
end
