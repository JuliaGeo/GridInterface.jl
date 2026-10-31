# Validates the O(1) dense HEALPix vertex numbering in DGGAdapter.
#   julia --project=prototype/env prototype/usecases/healpix_vertex_numbering/validate.jl
include(joinpath(@__DIR__, "..", "..", "adapters", "DGGAdapter.jl")); using .DGGAdapter
using GridInterface; import GridInterface as GrI
using Test
const VK = GrI.VertexKey
const HPS = DGG.HEALPixSystem()

function check_full(l)
    g = DGG.levelgrid(HPS, l); N = 1 << l; nc = GrI.ncells(g); nv = GrI.nvertices(g)
    @test nv == 12N^2 + 2
    keys = [GrI.vertex_key(g, v) for v in 1:nv]
    @test all(GrI.vertex_index(g, keys[v]) == v for v in 1:nv)
    nbit = Ref(0); offv = Set{Int}(); seen = falses(nv); pts = Vector{Any}(undef, nv)
    for i in 1:nc, k in 1:4
        ck = GrI.vertex_key(g, i, k)
        v = GrI.vertex_index(g, VK(i, k))
        @test v == GrI.vertex_index(g, ck)
        @test GrI.vertex_key(g, v) == ck
        ck == VK(i, k) && (@test !seen[v]; seen[v] = true)
        p = GrI.vertex_point(g, VK(i, k))
        if isassigned(pts, v)
            isequal(pts[v], p) || (nbit[] += 1; push!(offv, v))
            @test maximum(abs.(pts[v] .- p)) < 1e-12
        else
            pts[v] = p
        end
    end
    @test all(seen)                                   # canonical keys ↦ exactly 1:nv
    for v in 1:nv
        cs = GrI.vertex_cells(g, keys[v])
        @test keys[v].cell in cs && length(cs) in (3, 4)
    end
    # distinct dense numbers ↔ distinct points
    @test minimum(abs(pts[a][3] - pts[b][3]) + abs(pts[a][1] - pts[b][1]) + abs(pts[a][2] - pts[b][2]) for a in 1:min(nv, 400), b in 1:min(nv, 400) if a != b) > 1e-9
    println("level $l: $(nbit[]) corner keys at $(length(offv)) vertices not bit-identical to the owner's point (within 1e-12)")
end

for l in 0:3
    @testset "level $l" check_full(l)
end
@testset "level 6 (full)" check_full(6)

g6 = DGG.levelgrid(HPS, 6)
k6 = GrI.vertex_key(g6, 12345); GrI.vertex_index(g6, k6); GrI.vertex_index(g6, VK(777, 2))
alloc_key(g, v) = @allocated GrI.vertex_key(g, v)
alloc_idx(g, k) = @allocated GrI.vertex_index(g, k)
alloc_key(g6, 1); alloc_idx(g6, k6)
a1 = alloc_key(g6, 12345)
a2 = alloc_idx(g6, k6)
a3 = alloc_idx(g6, VK(777, 2))
@test a1 == a2 == a3 == 0
println("level 6 allocations: vertex_key $a1, vertex_index $a2 / $a3")

function bench(g, n)
    nv = GrI.nvertices(g); s = 0
    t1 = @elapsed for v in 1:max(1, nv ÷ n):nv; s += GrI.vertex_key(g, v).cell; end
    m = length(1:max(1, nv ÷ n):nv)
    ks = [GrI.vertex_key(g, v) for v in 1:max(1, nv ÷ n):nv]
    t2 = @elapsed for k in ks; s += GrI.vertex_index(g, k); end
    return t1 / m * 1e9, t2 / m * 1e9, s
end
bench(g6, 1000)
t1, t2, _ = bench(g6, 100_000)
println("level 6: vertex_key $(round(t1, digits=1)) ns/call, vertex_index $(round(t2, digits=1)) ns/call")

g12 = DGG.levelgrid(HPS, 12)
nv12 = GrI.nvertices(g12)
println("level 12: ncells $(GrI.ncells(g12)) nvertices $nv12")
@testset "level 12 spot" begin
    @test nv12 == 12 * 4096^2 + 2
    for v in (1, 2, 4096, 4097, 4096^2, 5 * 4096^2 + 3, 12 * 4096^2, nv12 - 1, nv12, 123_456_789)
        @test GrI.vertex_index(g12, GrI.vertex_key(g12, v)) == v
    end
    for i in (1, 201_326_592, 100_000_000), k in 1:4
        @test GrI.vertex_key(g12, GrI.vertex_index(g12, VK(i, k))) == GrI.vertex_key(g12, i, k)
    end
    @test alloc_key(g12, 123_456_789) == 0
    @test alloc_idx(g12, VK(100_000_000, 3)) == 0
    @test alloc_key(g12, nv12) == 0
end
t1, t2, _ = bench(g12, 100_000)
println("level 12: vertex_key $(round(t1, digits=1)) ns/call, vertex_index $(round(t2, digits=1)) ns/call")

if isdefined(GrI, :vertex_numbering_is_analytic)
    println("trait vertex_numbering_is_analytic(g12) = ", GrI.vertex_numbering_is_analytic(g12))
    GrI.dual(g6); t = @elapsed d = GrI.dual(g12)
    println("dual(level 12): $(round(t * 1e3, digits=2)) ms, table === nothing: $(d.table === nothing), ncells(dual) == nvertices: $(GrI.ncells(d) == nv12)")
else
    println("GrI.vertex_numbering_is_analytic not defined yet: dual check skipped")
end
