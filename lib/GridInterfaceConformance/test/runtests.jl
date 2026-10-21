# The suite proves itself: every mock fails exactly its family and nothing
# else, every GridInterface test grid passes, unclaimed suites are recorded as
# skipped, and the test_* forms emit one @test per law.
using Test
using GridInterface
using GridInterfaceConformance
using GridInterfaceConformance.Mocks
import GridInterface as GrI
import GridInterfaceConformance as Conf
using GridInterface.Structured
include(joinpath(@__DIR__, "..", "..", "..", "test", "grids.jl"))

failed_families(rep) = sort!(unique!([r.family for r in rep.results if r.status === :fail]); by=string)
failed_laws(rep) = [r.law for r in rep.results if r.status === :fail]
skipped_families(rep) = unique!([r.family for r in rep.results if r.status === :skip])
function expect_clean(rep)
    @test isempty(failed_families(rep))
    isempty(failed_families(rep)) || print_report(rep)
end

# A captured, silenced test set: (passes, failures, skips). A known-broken law
# (a zero-allocation law on Julia < 1.12) is a `Test.Broken` too, but not a skip.
function count_results(ts)
    passes = hasproperty(ts, :n_passed) ? ts.n_passed : 0
    fails = 0; skips = 0
    for r in ts.results
        if r isa Test.Pass
            passes += 1
        elseif r isa Test.Fail || r isa Test.Error
            fails += 1
        elseif r isa Test.Broken
            r.test_type === :skipped && (skips += 1)
        elseif r isa Test.AbstractTestSet
            p, f, k = count_results(r)
            passes += p; fails += f; skips += k
        end
    end
    return passes, fails, skips
end
# Julia 1.13 replaced the mutable testset stack with a scoped value entered
# through `Test.@with_testset`; a macro only reachable on one version has to be
# branched on at parse time.
@static if isdefined(Test, :push_testset)
    function with_testset(f, ts)
        Test.push_testset(ts)
        try
            f()
        finally
            Test.pop_testset()
        end
    end
else
    with_testset(f, ts) = Test.@with_testset ts f()
end
function capture(f)
    ts = Test.DefaultTestSet("capture")
    logger = Test.TestLogger()
    try
        Base.CoreLogging.with_logger(logger) do
            with_testset(ts) do
                redirect_stdout(devnull) do
                    f()
                end
            end
        end
    catch err
        err isa Test.TestSetException || rethrow()
    end
    return count_results(ts)..., [string(r.message) for r in logger.logs]
end

rg = RegularGrid(collect(0.0:1.0:5.0), collect(0.0:1.0:4.0))
lonlat = CellBasedGrid([GO.UnitSphereFromGeographic()((lo, la)) for lo in range(-180, 180; length=37), la in range(-90, 90; length=19)])
quad = levelgrid(QuadSystem(4), 3)
cube = CubeSystem(3)

@testset "GridInterfaceConformance" begin

@testset "every test grid passes the families it claims" begin
    for (name, g) in (("TinyPlanar", TinyPlanar(4, 3)), ("MinLonLat", MinLonLat(12, 6)), ("TinyCubeSphere", TinyCubeSphere(3)),
                      ("RegularGrid (Planar)", rg), ("CellBasedGrid lon/lat (Spherical)", lonlat),
                      ("QuadLevelGrid level 3", quad), ("CubeGrid level 2", levelgrid(cube, 2)),
                      ("GridAndTree of TinyCubeSphere", GridAndTree(TinyCubeSphere(2))),
                      ("Subset of TinyPlanar", subset(TinyPlanar(4, 3), [1, 2, 5, 6, 7])),
                      ("mockchunks of QuadLevelGrid", mockchunks(quad, LevelChunks(2))))
        expect_clean(check_grid(g; name))
    end
    # the polygon bag stores Vector{Any}: only the performance family may fail
    @test failed_families(check_grid(PolygonBag(); name="PolygonBag")) == [:performance]
    polys = [GI.Polygon([[(i + 0.2j, j), (i + 1 + 0.2j, j), (i + 1.2 + 0.2j, j + 1), (i + 0.2 + 0.2j, j + 1), (i + 0.2j, j)]]) for i in 0:7, j in 0:5]
    @test failed_families(check_grid(ExplicitPolygonGrid(polys); name="ExplicitPolygonGrid")) ⊆ [:performance]
    expect_clean(check_system(QuadSystem(5); name="QuadSystem"))
    expect_clean(check_system(cube; name="CubeSystem"))
    expect_clean(check_laziness(TinyPlanar(4, 3), TinyPlanar(16, 12)))
    expect_clean(check_laziness(levelgrid(cube, 1), levelgrid(cube, 3)))
end

@testset "unclaimed suites are recorded as skipped" begin
    rep = check_grid(PolygonBag(); name="PolygonBag")
    @test :dual in skipped_families(rep) && :geometry in skipped_families(rep)      # no rings, no dual
    @test any(r -> r.family === :dual && r.status === :skip && occursin("NonConforming", r.detail), rep.results)
    rep = check_grid(TupleLonLat(); name="TupleLonLat")
    @test :tree in skipped_families(rep)                                            # the spherical type contract failed first
    # test_grid turns every skip into a @test_skip and names it in one @info
    passes, fails, skips, msgs = capture(() -> test_grid(PolygonBag(); name="PolygonBag"))
    @test skips >= 2 && passes > 0
    @test any(m -> occursin("suite(s) skipped for PolygonBag", m) && occursin("dual:", m), msgs)
    # a grid that claims everything is told nothing
    _, _, skips, msgs = capture(() -> test_grid(TinyCubeSphere(2); name="TinyCubeSphere"))
    @test skips == 0 && !any(m -> occursin("skipped", m), msgs)
end

@testset "each mock fails exactly its family" begin
    for (mock, fam) in ((BrokenNeighbors(TinyPlanar(4, 3)), :neighbors),
                        (BrokenTree(TinyPlanar(4, 3)), :tree),
                        (BrokenIntent(quad), :tree),
                        (LeakyWrapper(quad), :tree),
                        (NarrowCapacity(rg), :tree),
                        (BrokenAreas(TinyPlanar(4, 3)), :geometry),
                        (BrokenExtent(TinyPlanar(4, 3)), :geometry),
                        (BrokenLocation(rg), :location),
                        (FalseAnalytic(TinyPlanar(4, 3)), :location),
                        (BrokenDual(rg), :dual))
        name = string(nameof(typeof(mock)))
        rep = check_grid(mock; name)
        @test failed_families(rep) == [fam]
        failed_families(rep) == [fam] || print_report(rep)
    end
    for (mock, fam) in ((BrokenCovering(cube), :hierarchy), (BrokenDescendantRange(cube), :hierarchy))
        name = string(nameof(typeof(mock)))
        rep = check_system(mock; name)
        @test failed_families(rep) == [fam]
        failed_families(rep) == [fam] || print_report(rep)
    end
    small, large = stateful_pair(TinyPlanar(4, 3), TinyPlanar(16, 12))
    rep = check_laziness(small, large)
    @test failed_families(rep) == [:performance]
    @test "laziness: cell_extent" in failed_laws(rep) && length(failed_laws(rep)) == 1
    # the stateful extent sweep allocates by design: a failure where the
    # zero-allocation laws hold, a known-broken law where they do not
    rep = check_grid(small; name="Stateful")
    @test failed_families(rep) == (Conf.ZERO_ALLOC ? [:performance] : Symbol[])
    @test any(r -> occursin("sweep over extent allocates nothing", r.law) && r.status === (Conf.ZERO_ALLOC ? :fail : :broken), rep.results)
    # the law each one trips, by name
    @test any(l -> occursin("symmetry", l), failed_laws(check_grid(BrokenNeighbors(TinyPlanar(4, 3)))))
    @test any(l -> occursin("leaves emit exactly 1:ncells once", l), failed_laws(check_grid(BrokenTree(TinyPlanar(4, 3)))))
    @test any(l -> occursin("ForSweeps", l), failed_laws(check_grid(BrokenIntent(quad))))
    @test failed_laws(check_grid(LeakyWrapper(quad))) == ["a wrapper does not demote a verb its parent provides to a fallback"]
    @test failed_laws(check_grid(NarrowCapacity(rg))) == ["subtree restricts treeify(grid, n)'s own tree (its integer method takes any Integer)"]
    @test failed_laws(check_grid(BrokenAreas(TinyPlanar(4, 3)))) == ["overridden area equals the geometric answer (rtol 0.001)"]
    @test failed_laws(check_grid(BrokenExtent(TinyPlanar(4, 3)))) == ["cell_extent contains the cell"]
    @test failed_laws(check_grid(FalseAnalytic(TinyPlanar(4, 3)))) == ["declared analytic location is provided (not the fallback)"]
    @test any(l -> occursin("fans wind CCW", l), failed_laws(check_grid(BrokenDual(rg))))
    @test any(l -> occursin("partition level", l), failed_laws(check_system(BrokenCovering(cube))))
    @test any(l -> occursin("children-recursion oracle", l), failed_laws(check_system(BrokenDescendantRange(cube))))
    # ...and through the @testset entry points
    _, fails, _, _ = capture(() -> test_grid(BrokenNeighbors(TinyPlanar(4, 3))))
    @test fails > 0
    _, fails, _, _ = capture(() -> test_system(BrokenCovering(cube)))
    @test fails > 0
    _, fails, _, _ = capture(() -> test_system(cube; levels=1:2))
    @test fails == 0
end

@testset "has_nonfallback_method" begin
    @test !has_nonfallback_method(cellat, TinyPlanar(2, 2), (0.5, 0.5))
    @test has_nonfallback_method(cellat, quad, (0.5, 0.5))
    @test has_nonfallback_method(cell_boundary, TinyPlanar(2, 2), 1)
    @test !has_nonfallback_method(cell_area, TinyPlanar(2, 2), 1)
    # wrappers forward: the question is asked of the parent
    @test has_nonfallback_method(cell_area, GridAndTree(quad), 1)
    @test !has_nonfallback_method(cell_area, GridAndTree(TinyPlanar(2, 2)), 1)
    @test has_nonfallback_method(treeify, GridAndTree(TinyPlanar(2, 2)), 3)      # the wrapper's own override
    @test !has_nonfallback_method(cellat, FalseAnalytic(TinyPlanar(2, 2)), (0.5, 0.5))
    @test has_nonfallback_method(cellat, BrokenLocation(rg), (0.5, 0.5))
    # a family supertype (AbstractStructuredGrid) is a provider method
    @test has_nonfallback_method(neighbors, rg, 1) && has_nonfallback_method(cellat, rg, (0.5, 0.5))
    @test !has_nonfallback_method(cellat, 1, (0.5, 0.5))                           # no subject
    @test has_nonfallback_method(children, cube, CubeCell(0, 0))
    @test !has_nonfallback_method(ancestor, cube, CubeCell(1, 0), 0)                # the generic parent loop
    @test :cell_area in fallback_verbs(TinyPlanar(2, 2)) && !(:cell_area in fallback_verbs(quad))
end

@testset "seeded sampling reproduces" begin
    a = check_grid(MinLonLat(36, 18); nsamples=16)
    b = check_grid(MinLonLat(36, 18); nsamples=16)
    @test [(r.law, r.status) for r in a.results] == [(r.law, r.status) for r in b.results]
    @test a.seed == DEFAULT_SEED
    @test occursin("(seed $(DEFAULT_SEED))", sprint(print_report, a))
    c = check_grid(MinLonLat(36, 18); nsamples=16, rng=Conf.Random.MersenneTwister(7))
    @test c.seed === nothing && !occursin("seed", first(split(sprint(print_report, c), '\n')))
    @test Conf.sample_positions(Conf._default_rng(), 1000, 10) == Conf.sample_positions(Conf._default_rng(), 1000, 10)
    s = Conf.sample_positions(Conf._default_rng(), 1000, 10)
    @test issorted(s) && allunique(s) && 1 in s && 1000 in s && length(s) <= 10
    @test Conf.sample_positions(Conf._default_rng(), 5, 10) == 1:5
end

end
