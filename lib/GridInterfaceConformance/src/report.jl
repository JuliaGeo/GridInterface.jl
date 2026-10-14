struct LawResult
    family::Symbol
    law::String
    status::Symbol      # :pass, :fail, :skip, :info
    detail::String
end

"""
    ConformanceReport

The result of [`check_grid`](@ref) or [`check_system`](@ref): `name`, one
`LawResult` per law in `results` (family, law, `:pass`/`:fail`/`:skip`/`:info`,
detail), the verbs answered by generic fallbacks in `fallbacks`, and measured
`metrics`. Summarise with [`family_status`](@ref) and [`print_report`](@ref).
"""
struct ConformanceReport
    name::String
    results::Vector{LawResult}
    fallbacks::Vector{Symbol}
    metrics::Dict{String,Any}
end
ConformanceReport(name) = ConformanceReport(name, LawResult[], Symbol[], Dict{String,Any}())

const FAMILIES = (:identity, :tree, :geometry, :location, :neighbors, :dual, :subsets,
    :hierarchy, :multilevel, :performance)

function law!(rep, family, law, ok::Bool, detail="")
    push!(rep.results, LawResult(family, law, ok ? :pass : :fail, string(detail)))
    return ok
end
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

`:pass`, `:fail` or `:skip` per family present in the report.
"""
function family_status(rep::ConformanceReport)
    out = Dict{Symbol,Symbol}()
    for r in rep.results
        r.status === :info && continue
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

Print the family status line, every failed, skipped and informational law
(every law when `verbose`), the fallbacks answering, and the metrics.
"""
function print_report(io::IO, rep::ConformanceReport; verbose=false)
    println(io, "== ", rep.name, " ==")
    st = family_status(rep)
    println(io, "  ", join(("$f=$(st[f])" for f in FAMILIES if haskey(st, f)), "  "))
    for r in rep.results
        (verbose || r.status in (:fail, :skip, :info)) || continue
        println(io, "  [", r.status, "] ", r.family, " / ", r.law, isempty(r.detail) ? "" : ": " * r.detail)
    end
    isempty(rep.fallbacks) || println(io, "  fallbacks answering: ", join(rep.fallbacks, ", "))
    for k in sort!(collect(keys(rep.metrics)))
        println(io, "  metric ", k, " = ", rep.metrics[k])
    end
end
print_report(rep::ConformanceReport; kw...) = print_report(stdout, rep; kw...)

# Emit the report as Test results (used by test_grid / test_system).
function emit_tests(rep::ConformanceReport)
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
                    elseif r.status === :skip
                        @test_skip r.law
                    end
                end
            end
        end
    end
end
