using Test
include("grids.jl")
using GridInterfaceConformance
const STI = GrI.STI

leaves(tree) = (out = Int[]; STI.depth_first_search(i -> (push!(out, i); nothing), _ -> true, tree); sort!(out))

@testset "GridInterface" begin

@testset "identity" begin
    g = TinyPlanar(3, 2)
    @test ncells(g) == 6
    @test cellindex(g, 4) == 4
    @test_throws BoundsError cellindex(g, 0)
    @test_throws BoundsError cellindex(g, 7)
    @test localindex(g, 5) == 5
    @test localindex(g, 9) === nothing
    @test cellindextypes(g) == (Int,)
    @test cellindex(g, 2, Int) == 2
    @test_throws ArgumentError cellindex(g, 2, Float64)
    @test globalindex(g, 3) == 3
    @test GrI.toindex(g, 2, Index(Local())) == 2
    @test GrI.toindex(g, 2, Index(Global())) == 2
    sys = QuadSystem(4); lg = levelgrid(sys, 2)
    @test cellindex(lg, 5) == QuadCell(2, 4)
    @test localindex(lg, QuadCell(2, 4)) == 5
    @test localindex(lg, QuadCell(1, 0)) === nothing
    @test level(cellindex(lg, 1)) == 2
    @test GrI.toindex(lg, 3, Index(QuadCell)) == QuadCell(2, 2)
end

@testset "geometry" begin
    g = TinyPlanar(3, 2)
    @test cellshape(g) isa SingleRing
    @test collect(cell_boundary(g, 1)) == [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
    poly = cell_geometry(g, 5)
    @test GI.trait(poly) isa GI.PolygonTrait
    @test GI.npoint(poly) == 5
    @test cell_area(g, 1) ≈ 1.0
    @test all(cell_centroid(g, 1) .≈ (0.5, 0.5))
    @test cell_extent(g, 6) == Extents.Extent(X=(2.0, 3.0), Y=(1.0, 2.0))
    s = TinyCubeSphere(2)
    @test sum(cell_area(s, i) for i in 1:ncells(s)) ≈ 4π rtol = 1e-10
    @test cell_extent(s, 1) isa SphericalCap
    bag = PolygonBag()
    @test cellshape(bag) isa GeneralCells
    @test_throws ArgumentError cell_boundary(bag, 1)
    @test cell_area(bag, 3) ≈ 16 - 4
end

@testset "location" begin
    g = TinyPlanar(3, 2)
    @test cellat(g, (0.5, 0.5)) == 1
    @test cellat(g, (2.5, 1.5)) == 6
    @test cellat(g, (1.0, 0.5)) == 1          # shared edge: lowest position
    @test cellat(g, (5.0, 5.0)) === nothing
    bag = PolygonBag()
    @test cellat(bag, (1.5, 1.5)) == 1        # overlap: lowest position
    @test cellat(bag, (2.5, 2.5)) == 2
    @test cellat(bag, (12.0, 2.0)) === nothing  # in the hole
    @test cellat(bag, (22.5, 0.5)) == 4         # second part of the multipolygon
    s = TinyCubeSphere(2)
    for i in 1:ncells(s)
        @test cellat(s, cell_centroid(s, i)) == i
    end
    lg = levelgrid(QuadSystem(4), 3)
    @test cellat(lg, (0.3, 0.7)) == GrI.fallback_cellat(GOC.Planar(), lg, (0.3, 0.7))
    @test (@allocated cellat(lg, (0.3, 0.7))) == 0
end

@testset "tree" begin
    g = TinyPlanar(4, 3)
    t = treeify(g)
    @test STI.isspatialtree(t)
    @test leaves(t) == 1:12
    @test leaves(subtree(g, [2, 5, 7])) == [2, 5, 7]
    gt = GridAndTree(g, t)
    @test treeify(gt) === t
    @test ncells(gt) == 12 && cell_area(gt, 3) == cell_area(g, 3) && cellat(gt, (1.5, 0.5)) == 2
    @test query(g, Extents.Extent(X=(0.5, 1.5), Y=(0.5, 0.6))) == [1, 2]
    sys = QuadSystem(5); lg = levelgrid(sys, 4)
    ct = treeify(lg)
    @test ct isa HierarchyCursor
    @test leaves(ct) == 1:ncells(lg)
    @test (@allocated treeify(lg)) == 0
    @test leaves(subtree(lg, [3, 40, 41, 200])) == [3, 40, 41, 200]
    s = subset(lg, 17:32)
    @test leaves(treeify(s)) == 1:16
end

@testset "subtree_indices under a root index" begin
    # STI's default: leaf indices depth-first, in traversal order
    walk(n) = STI.isleaf(n) ? Int[i for (i, _) in STI.child_indices_extents(n)] :
              reduce(vcat, (walk(c) for c in STI.getchild(n)); init=Int[])
    allnodes(n) = STI.isleaf(n) ? Any[n] : reduce(vcat, (allnodes(c) for c in STI.getchild(n)); init=Any[n])
    agrees(t) = all(n -> collect(subtree_indices(n)) == walk(n), allnodes(t))
    rq = RootedQuad(QuadSystem(5)); lg = levelgrid(rq.q, 3)
    t = GrI.HierarchyCursor(GrI.LevelFrontier(lg, rq, 3))
    groups = filter(n -> n.gl >= 1, allnodes(t))      # the top node and four groups of four roots
    @test length(groups) == 5
    for g in groups
        @test collect(subtree_indices(g)) == walk(g)
        @test length(subtree_indices(g)) == (g.gl == 2 ? 64 : 16)
    end
    @test agrees(t)
    @test sort(subtree_indices(t)) == 1:64
    keys = [3, 9, 17, 18, 40, 41, 60]
    for emit in (false, true)
        st = GrI.HierarchyCursor(GrI.SubsetFrontier(lg, rq, 3, keys, emit))
        @test agrees(st)
        @test sort(subtree_indices(st)) == (emit ? keys : 1:7)
    end
    m = MultiLevelSet(rq, [QuadCell(2, 0), QuadCell(3, 4), QuadCell(3, 5), QuadCell(3, 6), QuadCell(3, 7),
                           QuadCell(2, 2), QuadCell(2, 9), QuadCell(3, 60)])
    mt = treeify(m)
    @test mt isa HierarchyCursor && any(n -> n.gl >= 1, allnodes(mt))
    @test agrees(mt)
    @test sort(subtree_indices(mt)) == 1:ncells(m)
    @test agrees(subtree(m, [2, 5, 7]))
    gt = GridAndTree(lg, t)
    @test subtree_indices(gt) == subtree_indices(t)
    @test split_weight(gt) == split_weight(t) == 64
end

@testset "subsets and chunks" begin
    g = TinyPlanar(4, 3)
    s = subset(g, [2, 3, 6, 7])
    @test ncells(s) == 4 && parent(s) === g
    @test cellindex(s, 3) == 6
    @test globalindex(s, cellindex(s, 3)) == 6
    @test localindex(s, 7) == 4 && localindex(s, 1) === nothing
    @test cell_boundary(s, 1) == cell_boundary(g, 2)
    @test gridrelation(s) isa ParentGrid
    @test halo(s) == [1, 4, 5, 8, 9, 10, 11, 12]
    @test border(s) == [1, 2, 3, 4]
    @test isempty(interior(s))
    @test neighbors(s, 1; connectivity=Edge()) |> sort == [2, 3]
    lg = levelgrid(QuadSystem(3), 2)
    ch = chunks(lg)
    @test length(ch) == 4
    @test sort(reduce(vcat, [collect(parentpositions(c)) for c in ch])) == 1:16
    @test grid_extent(ch[2]) == Extents.Extent(X=(0.5, 1.0), Y=(0.0, 0.5))
end

@testset "neighbours" begin
    g = TinyPlanar(3, 3)
    @test sort(neighbors(g, 5)) == [1, 2, 3, 4, 6, 7, 8, 9]
    @test sort(neighbors(g, 5; connectivity=Edge())) == [2, 4, 6, 8]
    @test sort(neighbors(g, 1)) == [2, 4, 5]
    @test sort(ring(g, 1, 2)) == [3, 6, 7, 8, 9]
    @test neighbors(g, 1, Index(Local())) == neighbors(g, 1)
    s = TinyCubeSphere(2)
    @test length(neighbors(s, 1)) == 7               # cube corner: 3-valent vertex
    @test length(neighbors(s, 1; connectivity=Edge())) == 4
    lg = levelgrid(QuadSystem(3), 2)
    @test sort(neighbors(lg, 1)) == sort(GrI.fallback_neighbors(GOC.Planar(), lg, 1, Vertex()))
end

@testset "dual" begin
    g = TinyPlanar(3, 2)
    @test nvertices(g) == 12
    @test length(cell_vertices(g, 1)) == 4
    v = cell_vertices(g, 1)[3]                        # (1,1): shared by cells 1,2,4,5
    @test sort(collect(vertex_cells(g, v))) == [1, 2, 4, 5]
    D = dual(g)
    @test ncells(D) == 12
    @test length(cell_boundary(D, v)) == 4
    @test cell_centroid(D, v) == vertex_point(g, v)
    @test_throws ArgumentError dual(PolygonBag())
    s = TinyCubeSphere(2)
    @test nvertices(s) == 26                          # 6·(3²) − 12·3 + 8 = 26
    @test all(v -> length(vertex_cells(s, v)) in (3, 4), 1:26)
end

@testset "analytic vertex numbering" begin
    # explicit connectivity, no keys: nothing is built, the keys are derived
    m = ExplicitQuads()
    @test vertex_numbering_is_analytic(m) && !vertices_are_analytic(m)
    @test dual(m).table === nothing && GridAndVertices(m).table === nothing
    @test ncells(dual(m)) == 9 && nvertices(GridAndVertices(m)) == 9
    @test all(v -> vertex_index(m, vertex_key(m, v)) == v, 1:9)
    @test vertex_key(m, 5) == VertexKey(1, 3) && vertex_key(m, 2, 4) == VertexKey(1, 3) && vertex_key(m, 9) == VertexKey(4, 3)
    @test all(i -> [vertex_index(m, k) for k in cell_vertex_keys(m, i)] == cell_vertices(m, i), 1:4)
    @test all(i -> all(k -> vertex_key(m, i, k) == cell_vertex_keys(m, i)[k], 1:4), 1:4)
    @test all(v -> vertex_cells(m, vertex_key(m, v)) == vertex_cells(m, v), 1:9)
    @test all(v -> vertex_point(m, vertex_key(m, v)) == vertex_point(m, v), 1:9)
    @test all(v -> vertex_isopen(m, vertex_key(m, v)) == vertex_isopen(m, v), 1:9)
    @test count(v -> vertex_isopen(m, v), 1:9) == 8
    gv = GridAndVertices(m)
    @test cell_vertices(gv, 2) === cell_vertices(m, 2) && vertex_key(gv, 5) == VertexKey(1, 3)
    D = dual(m)
    @test length(cell_boundary(D, 5)) == 4 && cell_centroid(D, 5) == (1.0, 1.0) && cellat(D, (1.1, 0.9)) == 5
    @test dual_ring_cells(D, 5) == [4, 3, 1, 2] && dual_ring_cells(D, 2) == [0, 0, 2, 1, 0]
    rep = check_grid(m; name="ExplicitQuads")
    @test !any(r -> r.status === :fail, rep.results)
    any(r -> r.status === :fail, rep.results) && print_report(rep)
    @test any(r -> r.family === :dual && occursin("no vertex precomputation", r.law) && r.status === :pass, rep.results)
    @test any(r -> r.family === :dual && occursin("derived keyed verbs", r.law) && r.status === :pass, rep.results)
    # A planar regional dual fails the same two laws as dual(TinyPlanar(3, 2)):
    # its rings are built per call, and its own (second-order) median dual has
    # degenerate rings.
    dual_fails = ["dual cells are valid CCW rings (≥ 3 points)", "forwarding sweep over extent allocates nothing"]
    rep = check_grid(D; name="dual of ExplicitQuads")
    @test [r.law for r in rep.results if r.status === :fail] == dual_fails
    # the declaration without the methods is a provider error, not a silent fallback
    f = ForgetfulQuads()
    @test_throws ArgumentError nvertices(f)
    @test_throws ArgumentError cell_vertices(f, 1)
    @test_throws ArgumentError vertex_key(f, 1)
    e = try nvertices(f); nothing catch e; e end
    @test occursin("vertex_numbering_is_analytic", e.msg) && occursin("nvertices", e.msg)
    # analytic keys without a numbering: the numbering is derived and held by the dual
    k = KeyedPlanar(3, 2)
    @test vertices_are_analytic(k) && !vertex_numbering_is_analytic(k)
    @test nvertices(k) == 12 && dual(k).table isa GrI.VertexNumbering && GridAndVertices(k).table isa GrI.VertexNumbering
    @test all(v -> vertex_index(k, vertex_key(k, v)) == v, 1:12)
    @test collect(cell_vertices(k, 5)) == collect(cell_vertices(TinyPlanar(3, 2), 5))
    rep = check_grid(k; name="KeyedPlanar")
    @test !any(r -> r.status === :fail, rep.results)
    any(r -> r.status === :fail, rep.results) && print_report(rep)
    rep = check_grid(dual(k); name="dual of KeyedPlanar")
    @test [r.law for r in rep.results if r.status === :fail] == dual_fails
end

@testset "hierarchy and multi-level sets" begin
    sys = QuadSystem(6)
    @test ancestor(sys, QuadCell(3, 27), 1) == QuadCell(1, 27 >> 4)
    @test descendant_range(sys, QuadCell(1, 2), 3) == 33:48
    @test @inferred(levelgrid(sys, 3)) isa QuadLevelGrid
    cells = [QuadCell(1, 0), QuadCell(1, 1), QuadCell(2, 8), QuadCell(2, 9), QuadCell(3, 40), QuadCell(3, 41),
             QuadCell(3, 42), QuadCell(3, 43), QuadCell(2, 11), QuadCell(1, 3)]
    m = MultiLevelSet(sys, cells)
    @test ncells(m) == 10
    @test GrI.check_antichain(m)
    @test_throws ArgumentError MultiLevelSet(sys, [QuadCell(1, 0), QuadCell(2, 1)])
    @test localindex(m, QuadCell(3, 42)) !== nothing
    @test cellindex(m, localindex(m, QuadCell(2, 11))) == QuadCell(2, 11)
    @test system(m) === sys && gridrelation(m) isa ParentSystem
    @test_throws ArgumentError globalindex(m, QuadCell(1, 0))
    @test tessellation(m) isa Tessellating
    @test conformity(m) isa NonConforming
    @test_throws ArgumentError dual(m)
    @test leaves(treeify(m)) == 1:10
    @test location_is_analytic(m)
    @test cellat(m, (0.1, 0.1)) == localindex(m, QuadCell(1, 0))
    @test cellat(m, (0.6, 0.6)) == localindex(m, QuadCell(1, 3))
    @test cell_area(m, localindex(m, QuadCell(3, 40))) ≈ 1 / 64
    @test sum(cell_area(m, i) for i in 1:ncells(m)) ≈ 1.0
    @test @inferred(cell_centroid(m, 3)) isa Tuple{Float64,Float64}
    sweep(m) = (s = 0.0; for i in 1:ncells(m); s += cell_centroid(m, i)[1]; end; s)
    sweep(m)
    @test (@allocated sweep(m)) == 0
    @test (@allocated treeify(m)) == 0
end

@testset "round-one revisions" begin
    # required set of four, on the sphere, with collapsing pole corners
    g = MinLonLat(12, 6)
    @test sum(cell_area(g, i) for i in 1:ncells(g)) ≈ 4π rtol = 1e-10
    @test length(cell_boundary(g, 1)) == 4 && length(cell_corners(g, 1)) == 3   # south-pole cell
    @test length(cell_corners(g, 20)) == 4
    buf = UnitSphericalPoint{Float64}[]
    @test cell_boundary!(buf, g, 20) == collect(cell_boundary(g, 20)) && buf == collect(cell_boundary(g, 20))
    @test edgekind(g) isa GeodesicEdges
    @test cellsize(g) > 0
    @test cellat(g, cell_centroid(g, 30)) == 30
    @test cellat(g, GO.UnitSphereFromGeographic()((-179.0, 89.0))) == 61          # north-pole fan
    # the antipode of a cell is outside it (signed winding)
    c = cell_centroid(g, 30)
    @test !cell_contains(g, 30, UnitSphericalPoint(-c[1], -c[2], -c[3]))
    # spherical types are enforced
    @test_throws ArgumentError cellat(g, (10.0, 20.0))
    @test_throws ArgumentError cell_extent(TupleLonLat(), 1)
    rep = check_grid(TupleLonLat())
    @test any(r -> r.status === :fail && occursin("UnitSphericalPoint", r.law), rep.results)
    # dual: pole vertices are 12-cell fans, everything closes on a global grid
    @test nvertices(g) == 12 * 5 + 2
    pole = cell_vertices(g, 1)[1]
    @test length(vertex_cells(g, pole)) == 12 && !vertex_isopen(g, pole)
    @test !any(v -> vertex_isopen(g, v), 1:nvertices(g))
    D = dual(g)
    @test sum(v -> cell_area(D, v), 1:ncells(D)) ≈ 4π rtol = 1e-8
    @test cellat(D, vertex_point(g, 7)) == 7
    k = vertex_key(g, 7)
    @test vertex_index(g, k) == 7 && k isa VertexKey
    @test [vertex_index(g, kk) for kk in cell_vertex_keys(g, 20)] == collect(cell_vertices(g, 20))
    # GridAndTree, split_weight
    gt = GridAndTree(g)
    @test split_weight(gt) == split_weight(gt.tree)
    @test cellat(gt, cell_centroid(g, 3)) == 3

    # neighbours: disc, maxring, subsets clipped and unclipped
    p = TinyPlanar(5, 5)
    @test neighbors(p, 13, 1) |> sort == sort(neighbors(p, 13))
    @test length(neighbors(p, 13, 2)) == 24
    @test neighbors(p, 1; connectivity=Edge()) == [2, 6]                     # aligned CCW from edge (c1, c2)
    s = subset(p, [7, 8, 9, 12, 13, 14])
    @test collect(neighbors(s, 5)) |> sort == [1, 2, 3, 4, 6]
    loc = collect(neighbors(s, 5, Index(Local())))
    @test length(loc) == 8 && count(==(OUTSIDE), loc) == 3
    @test collect(neighbors(s, 5, Index(Global()))) == collect(neighbors(p, 13))
    @test halo(s) == halo(p, [7, 8, 9, 12, 13, 14])
    @test collect(interior(subset(p, collect(1:25)))) == 1:25
    bag = PolygonBag()
    @test neighbors(bag, 1) == [2]                                            # touching
    @test winding(bag) isa Unordered

    # collections: three spaces, algebra, cellsat, query predicates
    a = subset(p, [1, 2, 3, 4]); b = subset(p, [3, 4, 5, 6])
    @test parentpositions(union(a, b)) == 1:6
    @test parentpositions(intersect(a, b)) == [3, 4]
    @test parentpositions(setdiff(a, b)) == [1, 2]
    @test parentpositions(grow(subset(p, [13]), 1)) == [7, 8, 9, 12, 13, 14, 17, 18, 19]
    @test cellsat(bag, (1.5, 1.5)) == [1, 2]
    box = Extents.Extent(X=(0.5, 2.0), Y=(0.5, 2.0))
    @test query(p, Intersects(box)) == [1, 2, 3, 6, 7, 8, 11, 12, 13]
    @test query(p, Within(Extents.Extent(X=(0.0, 2.0), Y=(0.0, 1.0)))) == [1, 2]
    @test query(p, CentroidIn(box)) == [1, 2, 6, 7]
    @test query(p, Contains(Extents.Extent(X=(0.2, 0.4), Y=(0.2, 0.4)))) == [1]
    q = query(g, Intersects(Extents.Extent(X=(-10.0, 10.0), Y=(-10.0, 10.0))))
    @test !isempty(q) && all(i -> 2 <= (i - 1) ÷ 12 <= 3, q)
    @test query(g, CentroidIn(Extents.Extent(X=(-180.0, 180.0), Y=(60.0, 90.0)))) == 61:72

    # chunks: lazy, chunkat, subsets and multi-level sets
    sys = QuadSystem(6); lg = levelgrid(sys, 4)
    ch = chunks(lg, 1)
    @test length(ch) == 4 && chunkat(ch, 200) == 4 && 200 in parentpositions(ch[4])
    @test (@allocated chunks(lg, 2)) == 0
    sub = subset(lg, 60:140)
    sc = chunks(sub, 1)
    @test sum(ncells, sc) == 81
    m = MultiLevelSet(sys, [QuadCell(1, 0), QuadCell(2, 4), QuadCell(2, 5), QuadCell(2, 6), QuadCell(2, 7), QuadCell(1, 2), QuadCell(1, 3)])
    mc = chunks(m, 2)
    @test sum(ncells, mc) == 7 && length(mc) == 7                              # a coarse member is its own chunk
    em = expand_to_reference(m)
    @test ncells(em) == sum(c -> length(descendant_range(sys, c, m.reference_level)), m.cells) && parent(em) isa QuadLevelGrid
    @test leaves(subtree(lg, subset(lg, 17:32))) == 17:32
    @test system(lg) === sys && level(lg) == 4
end

@testset "intents and wrappers" begin
    maxleaf(t) = STI.isleaf(t) ? length(STI.child_indices_extents(t)) : maximum(maxleaf, STI.getchild(t); init=0)
    sys = QuadSystem(5); g = levelgrid(sys, 4)
    @test leaf_capacity(g, Balanced()) == 4 && leaf_capacity(g, ForSweeps()) == 8 && leaf_capacity(g, ForQueries()) == 1
    ts = treeify(g; leaf_capacity=ForSweeps())
    @test ts isa HierarchyCursor && ts.bucket == 8
    @test maxleaf(ts) <= leaf_capacity(g, ForSweeps()) && leaves(ts) == 1:ncells(g)
    @test maxleaf(treeify(g; leaf_capacity=ForQueries())) == 1
    @test maxleaf(treeify(g, 3)) <= 3 && leaves(treeify(g, 3)) == 1:ncells(g)
    @test treeify(g).bucket == 4 && treeify(g, Balanced()) == treeify(g)
    @test (@allocated treeify(g; leaf_capacity=ForSweeps())) == 0
    # the R-tree fallback takes the capacity as its node capacity (floored at 2)
    p = TinyPlanar(6, 5)
    @test maxleaf(treeify(p; leaf_capacity=ForSweeps())) <= 8 && leaves(treeify(p; leaf_capacity=ForSweeps())) == 1:30
    @test maxleaf(treeify(p; leaf_capacity=ForQueries())) <= 2 && leaves(treeify(p, 5)) == 1:30
    @test leaf_capacity(subset(p, 1:4), ForSweeps()) == 8
    # subtree of a multi-level set takes the capacity on the fallback too
    m = MultiLevelSet(sys, [QuadCell(1, 0), QuadCell(1, 1), QuadCell(1, 2), QuadCell(2, 12), QuadCell(2, 13), QuadCell(2, 14), QuadCell(2, 15)])
    @test leaves(treeify(m; leaf_capacity=ForSweeps())) == 1:7 && maxleaf(treeify(m; leaf_capacity=ForQueries())) == 1
    # GridAndTree: the stored tree for every hint and an integer
    t = treeify(g)
    gt = GridAndTree(g, t)
    for h in (Balanced(), ForSweeps(), ForQueries(), 2)
        @test treeify(gt; leaf_capacity=h) === t
        @test treeify(gt, h) === t
    end
    @test treeify(gt) === t && treeify(GridAndTree(g; leaf_capacity=ForSweeps())).bucket == 8
    @test split_weight(gt) == split_weight(t) && subtree_indices(gt) == subtree_indices(t)
    @test gt isa AbstractWrappedGrid && !(gt isa AbstractGrid) && isgrid(gt) && parent(gt) === g
    @test isgrid(g) && !isgrid(1) && !isgrid(QuadLevelGrid) && issystem(sys) && !issystem(g) &&
          iscellindex(QuadCell(1, 0)) && !iscellindex(3)
    @test grid_extent(gt) == grid_extent(g) && cellat(gt, (0.3, 0.7)) == cellat(g, (0.3, 0.7))
    # the specific methods win over the forwarding loop
    spec(f, T, args...) = Base.unwrap_unionall(which(f, Tuple{T,args...}).sig).parameters[2]
    @test spec(treeify, typeof(gt), ForQueries) === GridAndTree
    @test spec(treeify, typeof(gt), Int) === GridAndTree
    @test spec(ncells, typeof(gt)) === AbstractWrappedGrid
    # mockchunks: chunks is the mock, chunkat agrees, everything else forwards
    c = mockchunks(g, LevelChunks(2))
    @test c isa ChunkedGrid && chunks(c) == chunks(g, 2) && length(chunks(c)) == 16 && length(chunks(g)) == 4
    @test chunkat(c, 100) == chunkat(chunks(g, 2), 100) == 7 && chunkat(g, 100) == 2
    @test spec(chunks, typeof(c)) === ChunkedGrid
    for (f, args) in ((ncells, ()), (cellindex, (5,)), (cell_centroid, (7,)), (cell_extent, (9,)), (cell_area, (1,)),
                      (cellat, ((0.3, 0.7),)), (level, ()), (system, ()), (cellshape, ()), (best_manifold, ()),
                      (localindex, (QuadCell(4, 7),)), (maxneighbors, (Edge(),)), (grid_extent, ()))
        @test f(c, args...) == f(g, args...)
    end
    @test collect(neighbors(c, 6)) == collect(neighbors(g, 6)) && collect(ring(c, 6, 2)) == collect(ring(g, 6, 2))
    @test collect(cell_boundary(c, 3)) == collect(cell_boundary(g, 3)) && cellindex(c, 2, QuadCell) == cellindex(g, 2)
    @test leaves(treeify(c)) == 1:ncells(g) && leaves(subtree(c, [2, 9])) == [2, 9]
    @test chunks(mockchunks(g, nothing)) === nothing
    # composition: GridAndTree(mockchunks(subset(g, cells), LevelChunks(k)), t)
    s = subset(g, 17:48)
    cs = mockchunks(s, LevelChunks(3))
    tcs = treeify(cs; leaf_capacity=ForSweeps())     # a tree over the subset's own positions
    gts = GridAndTree(cs, tcs)
    @test treeify(gts) === tcs && leaves(treeify(gts)) == 1:32 && chunks(gts) === chunks(cs) && length(chunks(gts)) == 8
    @test ncells(gts) == 32 && cellindex(gts, 1) == cellindex(s, 1) && parentpositions(gts) == 17:48
    @test gridrelation(gts) isa ParentGrid && globalindex(gts, cellindex(gts, 3)) == 19
    @test cell_boundary(gts, 2) == cell_boundary(g, 18) && halo(gts) == halo(s)
    # the point location of a pair descends its stored tree
    pt = GridAndTree(p, treeify(p; leaf_capacity=ForQueries()))
    @test cellat(pt, (2.5, 3.5)) == cellat(p, (2.5, 3.5)) == 21 && cellsat(pt, (2.0, 3.0)) == [14, 15, 20, 21]
    @test query(pt, Extents.Extent(X=(0.5, 1.5), Y=(0.5, 0.6))) == query(p, Extents.Extent(X=(0.5, 1.5), Y=(0.5, 0.6)))
    # GridAndVertices holds the vertex table; the dual holds the same one
    q = TinyPlanar(3, 2)
    gv = GridAndVertices(q)
    @test nvertices(gv) == nvertices(q) == 12 && collect(cell_vertices(gv, 1)) == collect(cell_vertices(q, 1))
    @test (@allocated nvertices(gv)) == 0
    @test vertex_cells(gv, 5) == vertex_cells(q, 5) && vertex_key(gv, 5) == vertex_key(q, 5)
    D = dual(gv)
    @test D.table === gv.table && ncells(D) == 12 && cellat(D, vertex_point(q, 3)) == 3
    @test dual(q).table isa GrI.VertexTable
    # no global state: the greps of the outline, on the sources
    files(d) = [joinpath(r, f) for (r, _, fs) in walkdir(joinpath(@__DIR__, "..", d)) for f in fs if endswith(f, ".jl")]
    srcs = vcat(files("src"), files("ext"), files("lib"))
    @test !any(f -> occursin(r"clear_caches!|IdDict|ReentrantLock", read(f, String)), srcs)
    @test !any(f -> occursin(r"isgrid\(::Type|issystem\(::Type|iscellindex\(::Type", read(f, String)), files("src"))
end

@testset "conformance" begin
    for (name, g) in (("TinyPlanar", TinyPlanar(4, 3)), ("TinyCubeSphere", TinyCubeSphere(3)),
                      ("MinLonLat (four methods)", MinLonLat(12, 6)),
                      ("Subset of TinyPlanar", subset(TinyPlanar(4, 3), [1, 2, 5, 6, 7])))
        test_grid(g; name)
    end
    # The polygon bag stores `Vector{Any}`: the conformance suite must catch the
    # type instability (performance family) and nothing else.
    rep = check_grid(PolygonBag(); name="PolygonBag")
    st = family_status(rep)
    @test st[:performance] === :fail
    @test all(f -> st[f] !== :fail, setdiff(keys(st), [:performance]))
    @test st[:dual] === :pass && st[:neighbors] === :pass   # dual: refusal is the law that passes; touching neighbours
    # The dual of the cube builds each ring per call from stack-allocated fans:
    # even the allocation-free extent sweep passes.
    rep = check_grid(dual(TinyCubeSphere(2)); name="dual of TinyCubeSphere")
    fails = [r.law for r in rep.results if r.status === :fail]
    @test isempty(fails)
    rep = test_system(QuadSystem(5); name="QuadSystem")
    @test any(r -> r.family === :performance && occursin("laziness: cellat", r.law) && r.status === :pass, rep.results)
    rep = check_laziness(TinyPlanar(4, 3), TinyPlanar(16, 12))
    @test !any(r -> r.status === :fail, rep.results)
    any(r -> r.status === :fail, rep.results) && print_report(rep)
    test_grid(GridAndTree(TinyCubeSphere(2)); name="GridAndTree of TinyCubeSphere")
end

include("structured.jl")

end
