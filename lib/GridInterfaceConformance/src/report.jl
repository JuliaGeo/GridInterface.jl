struct LawResult
    family::Symbol
    law::String
    status::Symbol      # :pass, :fail, :broken, :skip, :info
    detail::String
end

"""
    ConformanceReport

The result of [`check_grid`](@ref) or [`check_system`](@ref): `name`, one
`LawResult` per law in `results` (family, law, status, detail), the verbs
answered by generic fallbacks in `fallbacks`, measured `metrics`, and the
sampling `seed` (printed in the header so a failure reproduces). Summarise
with [`family_status`](@ref) and [`print_report`](@ref).

| status | meaning |
|---|---|
| `:pass`, `:fail` | the law held, or did not |
| `:broken` | a zero-allocation law that allocated on Julia < 1.12, where the guarantee is known not to hold |
| `:skip` | a family the grid does not claim |
| `:info` | a measurement, never a verdict |
"""
struct ConformanceReport
    name::String
    results::Vector{LawResult}
    fallbacks::Vector{Symbol}
    metrics::Dict{String,Any}
    seed::Union{Nothing,Int}      # the sampling seed, or nothing when the caller supplied an rng
end
ConformanceReport(name; seed=nothing) = ConformanceReport(name, LawResult[], Symbol[], Dict{String,Any}(), seed)

# The rng a check samples with: the caller's, or a fresh generator on `seed`
# (recorded in the report so the run reproduces).
function sampling_rng(seed, rng)
    rng === nothing && return _default_rng(seed), Int(seed)
    return rng, nothing
end

const FAMILIES = (:identity, :tree, :geometry, :location, :neighbors, :dual, :subsets,
    :hierarchy, :multilevel, :performance)

function law!(rep, family, law, status::Symbol, detail="")
    push!(rep.results, LawResult(family, law, status, string(detail)))
    return status === :pass
end
law!(rep, family, law, ok::Bool, detail="") = law!(rep, family, law, ok ? :pass : :fail, detail)

# Zero-allocation laws hold from Julia 1.12 on. An older compiler leaves a few
# bytes behind the same code, so there a non-zero count is `:broken`: printed,
# counted against no family, a `@test_broken` in the test forms.
const ZERO_ALLOC = VERSION >= v"1.12"
allocation_verdict(bytes::Integer) = bytes == 0 ? :pass : ZERO_ALLOC ? :fail : :broken
skip!(rep, family, why) = push!(rep.results, LawResult(family, "(family)", :skip, why))
info!(rep, family, law, detail) = push!(rep.results, LawResult(family, law, :info, string(detail)))
metric!(rep, k, v) = (rep.metrics[k] = v)

# A law that throws is a failure carrying the exception, never an abort.
function guarded!(f, rep, family, law)
    try
        r = f()
        r isa Bool && return law!(rep, family, law, r)
        r isa Tuple && return law!(rep, family, law, r[1], r[2])
        return law!(rep, family, law, true)
    catch e
        msg = sprint(showerror, e)
        return law!(rep, family, law, false, "threw: " * first(msg, 300))
    end
end

"""
    family_status(rep) -> Dict{Symbol,Symbol}

`:pass`, `:fail` or `:skip` per family present in the report. `:broken` and
`:info` laws count against no family.
"""
function family_status(rep::ConformanceReport)
    out = Dict{Symbol,Symbol}()
    for r in rep.results
        r.status in (:info, :broken) && continue
        cur = get(out, r.family, :none)
        if r.status === :fail
            out[r.family] = :fail
        elseif r.status === :pass
            cur in (:none, :skip) && (out[r.family] = :pass)
        elseif r.status === :skip
            cur === :none && (out[r.family] = :skip)
        end
    end
    return out
end

"""
    print_report([io], rep::ConformanceReport; verbose = false)

Print the family status line, every failed, broken, skipped and
informational law (every law when `verbose`), the fallbacks answering, and
the metrics.
"""
function print_report(io::IO, rep::ConformanceReport; verbose=false)
    println(io, "== ", rep.name, " ==", rep.seed === nothing ? "" : "  (seed $(rep.seed))")
    st = family_status(rep)
    println(io, "  ", join(("$f=$(st[f])" for f in FAMILIES if haskey(st, f)), "  "))
    for r in rep.results
        (verbose || r.status in (:fail, :broken, :skip, :info)) || continue
        println(io, "  [", r.status, "] ", r.family, " / ", r.law, isempty(r.detail) ? "" : ": " * r.detail)
    end
    isempty(rep.fallbacks) || println(io, "  fallbacks answering: ", join(rep.fallbacks, ", "))
    for k in sort!(collect(keys(rep.metrics)))
        println(io, "  metric ", k, " = ", rep.metrics[k])
    end
end
print_report(rep::ConformanceReport; kw...) = print_report(stdout, rep; kw...)

# Emit the report as Test results (used by test_grid / test_system): one @test
# per law, a @test_broken per known-broken law, a @test_skip per suite the grid
# does not claim, and, since `Test` never displays a skip's reason and a
# summary cannot tell a skip from a known failure, one @info naming every skip
# and every broken law at the end (DGG's reporting).
function emit_tests(rep::ConformanceReport)
    skips = String[]
    broken = String[]
    @testset "$(rep.name)" begin
        for fam in FAMILIES
            rs = filter(r -> r.family === fam, rep.results)
            isempty(rs) && continue
            @testset "$fam" begin
                for r in rs
                    if r.status === :pass
                        @test true
                    elseif r.status === :fail
                        @test r.status !== :fail
                        @info "law failed" family = fam law = r.law detail = r.detail
                    elseif r.status === :broken
                        push!(broken, "$fam / $(r.law): $(r.detail)")
                        @test_broken r.status !== :broken
                    elseif r.status === :skip
                        reason = "$fam: $(r.detail)"
                        push!(skips, reason)
                        @test_skip reason
                    end
                end
            end
        end
        isempty(skips) || @info "conformance: $(length(skips)) suite(s) skipped for $(rep.name)\n" *
            join(("  " * s for s in skips), "\n")
        isempty(broken) || @info "conformance: $(length(broken)) law(s) known-broken on Julia $(VERSION) for $(rep.name)\n" *
            join(("  " * s for s in broken), "\n")
    end
    return nothing
end
