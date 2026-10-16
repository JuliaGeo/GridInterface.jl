# The structured family: CR's grids on the face hooks, the quadtree, tiles,
# point location, vertices and neighbours across seams.
using Random
using GridInterface.Structured

walk(n) = STI.isleaf(n) ? Int[i for (i, _) in STI.child_indices_extents(n)] :
          reduce(vcat, (walk(c) for c in STI.getchild(n)); init=Int[])
allnodes(n) = STI.isleaf(n) ? Any[n] : reduce(vcat, (allnodes(c) for c in STI.getchild(n)); init=Any[n])
maxleaf(t) = STI.isleaf(t) ? length(STI.child_indices_extents(t)) : maximum(maxleaf, STI.getchild(t); init=0)
_usp(x, y, z) = (r = sqrt(x^2 + y^2 + z^2); UnitSphericalPoint(x / r, y / r, z / r))
_sameptr(a, b) = all(abs.(a .- b) .< 1e-9)

@testset "structured" begin
    rg = RegularGrid(collect(0.0:1.0:20.0), collect(0.0:1.0:10.0))
    lons = range(-180, 180; length=37); lats = range(-90, 90; length=19)
    cb = CellBasedGrid([GO.UnitSphereFromGeographic()((lo, la)) for lo in lons, la in lats])
    cube = TinyCubeSphere(3)
    @test ncells(rg) == 200 && ncells(cb) == 648 && ncells(cube) == 54
    @test best_manifold(rg) isa GOC.Planar && best_manifold(cb) isa GOC.Spherical
    @test Structured.isstructured(rg) && Structured.isstructured(cube) && !Structured.isstructured(TinyPlanar(2, 2))
    @test cellshape(rg) isa ConvexCells && cellshape(cb) isa ConvexCells
    @test collect(cell_boundary(rg, 22)) == [(1.0, 1.0), (2.0, 1.0), (2.0, 2.0), (1.0, 2.0)]
    @test length(cell_corners(cb, 1)) == 3 && length(cell_corners(cb, 40)) == 4     # pole corner collapsed
    @test length(cell_boundary(cb, 648)) == 3
    @test cell_extent(rg, 22) == Extents.Extent(X=(1.0, 2.0), Y=(1.0, 2.0)) && !extent_is_expensive(rg)
    @test cell_extent(cb, 40) isa SphericalCap && extent_is_expensive(cb)
    @test isglobal(RegularGrid(GOC.Spherical(), collect(lons), collect(lats))) && !isglobal(rg)

    # conformance: both grids pass every family they claim
    for (name, g) in (("RegularGrid 20x10 (Planar)", rg), ("CellBasedGrid lon/lat 36x18 (Spherical)", cb))
        rep = test_grid(g; name)
        @test any(r -> r.family === :neighbors && occursin("alignment", r.law) && r.status === :pass, rep.results)
        @test any(r -> r.family === :dual && occursin("keys round-trip", r.law) && r.status === :pass, rep.results)
    end

    # cellat on the RegularGrid is the axis floor
    @test location_is_analytic(rg) && !location_is_analytic(cb) && !location_is_analytic(cube)
    rng = MersenneTwister(7)
    floor_ok = 0
    for _ in 1:1000
        x = 20rand(rng); y = 10rand(rng)
        i = min(floor(Int, x) + 1, 20); j = min(floor(Int, y) + 1, 10)
        cellat(rg, (x, y)) == i + (j - 1) * 20 && (floor_ok += 1)
    end
    @test floor_ok == 1000
    @test cellat(rg, (-1.0, 2.0)) === nothing && cellat(rg, (20.0, 10.0)) == 200 && cellat(rg, (0.0, 0.0)) == 1
    @test (@allocated cellat(rg, (3.3, 4.4))) == 0
    # a spherical lon/lat RegularGrid locates through the chart, repaired against the geodesic edges
    srg = RegularGrid(GOC.Spherical(; radius=1.0), collect(lons), collect(lats))
    @test location_is_analytic(srg) && cellat(srg, cell_centroid(srg, 100)) == 100
    @test cellat(srg, GO.UnitSphereFromGeographic()((179.9, 0.5))) == 36 + 9 * 36
    @test cellat(srg, GO.UnitSphereFromGeographic()((-180.0, -90.0))) == 1
    @test all(i -> GrI.cell_contains(srg, cellat(srg, cell_centroid(srg, i)), cell_centroid(srg, i)), 1:ncells(srg))

    # on a warped CellBasedGrid the walk agrees with the tree fallback
    warp = CellBasedGrid([GO.UnitSphereFromGeographic()((lo + 6sind(3la), la + 4sind(2lo)))
                          for lo in range(-60, 60; length=25), la in range(-40, 50; length=19)])
    gt = GridAndTree(warp, treeify(warp; leaf_capacity=ForQueries()))
    agree = 0; known = 0
    for _ in 1:1000
        k = rand(rng, 1:ncells(warp))
        w = rand(rng, 4); c = cell_boundary(warp, k)
        s = w[1] * c[1] + w[2] * c[2] + w[3] * c[3] + w[4] * c[4]
        p = _usp(s[1], s[2], s[3])
        a = cellat(warp, p); b = cellat(gt, p)
        a == b && (agree += 1)
        a == k && (known += 1)
    end
    @test agree == 1000 && known == 1000
    @test cellat(warp, UnitSphericalPoint(0.0, 0.0, -1.0)) === nothing           # far outside: the walk stalls, the tree says no
    @test cellat(cube, cell_centroid(cube, 40)) == 40                              # the walk crosses seams

    # the quadtree: leaves emit 1:ncells once; subtree_indices of a node is its leaf walk; capacities
    for g in (rg, cb, cube, warp)
        t = treeify(g)
        @test STI.isspatialtree(t) && sort!(walk(t)) == 1:ncells(g)
        @test subtree_indices(t) == walk(t)
        @test all(n -> subtree_indices(n) == walk(n) && split_weight(n) == length(walk(n)), allnodes(t))
        @test all(n -> STI.isleaf(n) || subtree_indices(n) == reduce(vcat, subtree_indices.(collect(STI.getchild(n)))), allnodes(t))
        @test all(n -> sort(collect(Structured.block_indices(n))) == sort(walk(n)), allnodes(t))
        @test maxleaf(treeify(g; leaf_capacity=ForSweeps())) <= 16 && maxleaf(treeify(g; leaf_capacity=ForQueries())) <= 4
        @test maxleaf(treeify(g; leaf_capacity=ForQueries())) == 4 && sort!(walk(treeify(g, 16))) == 1:ncells(g)
        @test (@allocated treeify(g)) == 0
    end
    t = treeify(rg)
    @test t isa Structured.TopDownQuadtreeCursor && t.leafsize == (2, 2) && STI.nchild(t) == 4
    @test treeify(rg; leaf_capacity=ForSweeps()).leafsize == (4, 4)
    @test !STI.node_extent_is_expensive(typeof(t)) && STI.node_extent_is_expensive(typeof(treeify(cb)))
    @test STI.node_extent(t) == Extents.Extent(X=(0.0, 20.0), Y=(0.0, 10.0))
    ct = treeify(cube)
    @test ct isa Structured.FaceRootTree && STI.nchild(ct) == 6 && STI.node_extent(ct).radius > π / 2
    @test all(c -> c isa Structured.TopDownQuadtreeCursor, STI.getchild(ct))
    # every node's cap covers its cells
    for n in allnodes(treeify(cb))
        e = STI.node_extent(n)
        @test all(i -> all(p -> GO.UnitSpherical.spherical_distance(e.point, p) <= e.radius, cell_boundary(cb, i)), walk(n))
    end

    # tiles: a lazy partition, O(1) extents, subtree and tree
    ch = chunks(rg, (8, 4))
    @test ch isa Structured.TileChunks && length(ch) == 9 && length(chunks(rg)) == 1 && length(chunks(cube)) == 6
    @test sort!(reduce(vcat, [collect(parentpositions(c)) for c in ch])) == 1:200
    @test all(k -> all(p -> chunkat(ch, p) == k, parentpositions(ch[k])), 1:9)
    @test chunkat(ch, 200) == 9 && chunkat(rg, 200) == 1 && chunkat(chunks(cube), 54) == 6
    c5 = ch[5]
    @test c5 isa Subset && parentpositions(c5) isa Structured.TileIndices && issorted(parentpositions(c5))
    @test grid_extent(c5) == reduce(Extents.union, (cell_extent(rg, p) for p in parentpositions(c5)))
    @test sort!(walk(subtree(rg, c5))) == collect(parentpositions(c5)) == sort!(walk(subtree(rg, parentpositions(c5))))
    @test sort!(walk(treeify(c5))) == 1:ncells(c5) && sort!(walk(treeify(c5; leaf_capacity=ForSweeps()))) == 1:ncells(c5)
    @test subtree_indices(treeify(c5)) == walk(treeify(c5)) && Structured.block_indices(treeify(c5)) == 1:ncells(c5)
    @test subtree_indices(subtree(rg, c5)) == walk(subtree(rg, c5))
    @test (@allocated grid_extent(c5)) == 0 && (@allocated subtree(rg, c5)) == 0 && (@allocated treeify(c5)) == 0
    cc = chunks(cb)[1]
    @test grid_extent(cc) isa SphericalCap && all(i -> all(p -> GO.UnitSpherical.spherical_distance(grid_extent(cc).point, p) <= grid_extent(cc).radius,
        cell_boundary(cc, i)), 1:ncells(cc))
    @test 17 in parentpositions(c5) || !(17 in parentpositions(c5))    # membership is defined
    @test all(p -> (p in parentpositions(c5)) == (p in collect(parentpositions(c5))), 1:200)

    # neighbours: aligned, counter-clockwise, across periodic seams and poles
    @test neighbors(rg, 22; connectivity=Edge()) == [2, 23, 42, 21] && neighbors(rg, 22) == [2, 3, 23, 43, 42, 41, 21, 1]
    @test neighbors(rg, 1) == [2, 22, 21] && neighbors(rg, 1; connectivity=Edge()) == [2, 21]
    @test winding(rg) isa CounterClockwise && maxneighbors(rg, Edge()) == 4 && neighbors_are_analytic(rg)
    @test neighbors(cb, 1; connectivity=Edge()) == [2, 37, 36]                 # south-pole cell: no neighbour across the pole
    @test neighbors(cb, 36; connectivity=Edge()) == [1, 72, 35]                # periodic seam at the dateline
    @test length(neighbors(cb, 1)) == 38 && 36 in neighbors(cb, 1) && 3 in neighbors(cb, 1)   # the pole fan
    @test sort(neighbors(cb, 40)) == [3, 4, 5, 39, 41, 75, 76, 77]
    # the seams are declared by arithmetic: periodic longitude, collapsed pole rows
    @test cb.periodic == (true, false) && cb.collapsed == (true, true) && rg.periodic == (false, false) && !any(rg.collapsed)
    @test Structured.facecross(cb, 1, 36, 2, 2) == (1, 1, 2, 4) && Structured.facecross(cb, 1, 1, 2, 4) == (1, 36, 2, 2)
    @test Structured.facecross(cb, 1, 5, 1, 1) === nothing && Structured.facecross(rg, 1, 20, 5, 2) === nothing
    @test neighbors(cb, 72; connectivity=Edge()) == [36, 37, 108, 71]           # dateline cell (36, 2): east is column 1
    @test neighbors(cb, 72) == [36, 1, 37, 73, 108, 107, 71, 35]                # with the two diagonal wraps
    @test Structured.facefan(cb, 1, 7, 1) == 36:-1:1 && Structured.facefan(cb, 1, 7, 19) == 613:648   # pole fans: every polar cell
    @test Structured.facefan(cb, 1, 7, 5) isa SmallVector{4,Int} && Structured.facefan(rg, 1, 1, 1) == [1]
    @test maxneighbors(cb, Vertex()) === nothing && maxneighbors(rg, Vertex()) == 8 && maxneighbors(cube, Vertex()) == 8
    @test neighbors(cb, 40; connectivity=Edge()) isa SmallVector{4,Int} && neighbors(cb, 40) isa SmallVector{8,Int} && neighbors(cb, 72) isa SmallVector{8,Int}
    @test neighbors(cb, 1) isa Vector{Int} && neighbors(cube, 1) isa SmallVector{8,Int}
    @test (@allocated neighbors(cb, 40; connectivity=Edge())) == 0 && (@allocated neighbors(cb, 72; connectivity=Edge())) == 0
    @test (@allocated neighbors(cb, 40)) == 0 && (@allocated neighbors(cb, 72)) == 0 && (@allocated neighbors(cube, 14)) == 0
    # declared seams beat detection: a wrap column off by 1e-9 reads as an open border unless the grid says otherwise
    let pts = [GO.UnitSphereFromGeographic()((lo, la)) for lo in lons, la in lats]
        pts[end, :] .= [GO.UnitSphereFromGeographic()((180 - 1e-9, la)) for la in lats]
        regional = CellBasedGrid(pts); global_ = CellBasedGrid(pts; periodic=(true, false))
        @test regional.periodic == (false, false) && global_.periodic == (true, false) && global_.collapsed == (true, true)
        @test CellBasedGrid(pts; collapsed=(false, false)).collapsed == (false, false)
        @test CellBasedGrid(GOC.Spherical(), pts; periodic=(true, false), collapsed=(true, true)).periodic == (true, false)
        @test neighbors(regional, 72; connectivity=Edge()) == [36, 108, 71] && neighbors(global_, 72; connectivity=Edge()) == [36, 37, 108, 71]
        @test neighbors(regional, 37; connectivity=Edge()) == [1, 38, 73] && neighbors(global_, 37; connectivity=Edge()) == [1, 38, 73, 72]
        @test Structured.facecross(regional, 1, 36, 2, 2) === nothing && Structured.facecross(global_, 1, 36, 2, 2) == (1, 1, 2, 4)
        @test RegularGrid(GOC.Spherical(), collect(range(-180, 180 - 1e-6; length=37)), collect(lats)).periodic == (false, false)
        @test RegularGrid(GOC.Spherical(), collect(range(-180, 180 - 1e-6; length=37)), collect(lats); periodic=(true, false)).periodic == (true, false)
        @test RegularGrid(collect(0.0:1.0:20.0), collect(0.0:1.0:10.0); periodic=(true, false), collapsed=(false, false)).periodic == (true, false)
    end
    # a regional planar grid: open borders, no wrap
    @test length(neighbors(rg, 10; connectivity=Edge())) == 3 && length(neighbors(rg, 20; connectivity=Edge())) == 2
    @test neighbors(rg, 20; connectivity=Edge()) == [40, 19] && neighbors(rg, 200; connectivity=Edge()) == [180, 199]
    @test all(k -> !(k + 1 in neighbors(rg, k)), 20:20:200)
    # the cube: every face edge has a neighbour on another face sharing the edge's two vertices
    seam = 0
    for i in 1:ncells(cube)
        js = neighbors(cube, i; connectivity=Edge()); cs = collect(cell_corners(cube, i))
        @test length(js) == 4
        for (t, j) in enumerate(js)
            cj = collect(cell_corners(cube, j))
            @test any(q -> _sameptr(q, cs[t]), cj) && any(q -> _sameptr(q, cs[mod1(t + 1, 4)]), cj)
            Structured.facecell(cube, j)[1] != Structured.facecell(cube, i)[1] && (seam += 1)
        end
    end
    @test seam == 12 * 3 * 2                                                    # 12 cube edges × 3 cells, both sides
    @test count(i -> length(neighbors(cube, i)) == 7, 1:ncells(cube)) == 24 && all(i -> length(neighbors(cube, i)) in (7, 8), 1:ncells(cube))   # 3 cells at each of the 8 cube corners
    @test neighbors_are_analytic(cube) && neighbors_are_analytic(cb)
    # the seam table pairs each cube edge with one other face, reversed; crossing back lands on the cell left
    for f in 1:6, s in 1:4
        f2, s2, rev = _CUBE_SEAMS[f][s]
        @test f2 != f && rev && _CUBE_SEAMS[f2][s2] == (f, s, true)
        for p in 1:3
            i, j = _cube_sidecell(3, s, p)
            x = Structured.facecross(cube, f, i, j, s)
            @test x[1] == f2 && x[4] == s2 && Structured.facecross(cube, x[1], x[2], x[3], x[4]) == (f, i, j, s)
        end
    end
    for i in (1, 14, 40)
        js = neighbors(cube, i); c = cell_centroid(cube, i)
        @test length(unique(js)) == length(js) && !(i in js)
        @test all(j -> i in neighbors(cube, j), js)
    end
    # vertices: analytic keys; cube corners are 3-cell fans; pole fans close
    @test vertices_are_analytic(cube) && nvertices(cube) == 56                 # V = 2 + E - F = 2 + 108 - 54
    @test count(v -> length(vertex_cells(cube, v)) == 3, 1:nvertices(cube)) == 8
    @test all(v -> !vertex_isopen(cube, v), 1:nvertices(cube))
    @test nvertices(cb) == 36 * 17 + 2 && length(vertex_cells(cb, vertex_key(cb, 1, 1))) == 36
    @test vertex_key(cb, 36, 1) == VertexKey(1, 1) && vertex_key(cb, 36, 2) == VertexKey(1, 3) && vertex_key(rg, 22, 1) == VertexKey(1, 3)
    @test (@allocated vertex_key(rg, 22, 1)) == 0 && (@allocated vertex_key(cb, 72, 2)) == 0
    # keys across the dateline and across a cube seam are canonical: they round-trip through vertex_index
    @test vertex_key(cb, 72, 2) == VertexKey(1, 3) && vertex_key(cb, 72, 3) == VertexKey(37, 4)
    @test all(c -> vertex_key(cb, vertex_index(cb, VertexKey(72, c))) == vertex_key(cb, 72, c), 1:4)
    @test all(c -> vertex_key(cube, vertex_index(cube, VertexKey(12, c))) == vertex_key(cube, 12, c), 1:4)
    @test Set(vertex_cells(cb, vertex_key(cb, 72, 2))) == Set([36, 72, 37, 1]) && vertex_cells(cb, vertex_key(cb, 72, 2)) isa SmallVector{4,Int}
    @test nvertices(rg) == 21 * 11 && count(v -> vertex_isopen(rg, v), 1:nvertices(rg)) == 2 * 21 + 2 * 9
    D = dual(cube)
    @test sum(v -> cell_area(D, v), 1:ncells(D)) ≈ 4π rtol = 1e-8

    # ExplicitPolygonGrid: positions, tiles and the tree from the lattice, cells from the polygons
    polys = [GI.Polygon([[(i + 0.2j, j), (i + 1 + 0.2j, j), (i + 1.2 + 0.2j, j + 1), (i + 0.2 + 0.2j, j + 1), (i + 0.2j, j)]])
             for i in 0:7, j in 0:5]
    ep = ExplicitPolygonGrid(polys)
    @test cellshape(ep) isa GeneralCells && ncells(ep) == 48 && cell_geometry(ep, 10) === polys[2, 2]
    @test_throws ArgumentError cell_boundary(ep, 1)
    @test_throws ArgumentError Structured.facevertex(ep, 1, 1, 1)
    @test sort!(walk(treeify(ep))) == 1:48 && cell_extent(ep, 1) == GI.extent(polys[1, 1])
    @test cellat(ep, (0.5, 0.5)) == 1 && sort(neighbors(ep, 10)) == [1, 2, 3, 9, 11, 17, 18, 19]     # touching (GeneralCells)
    rep = check_grid(ep; name="ExplicitPolygonGrid sheared 8x6")
    @test all(r -> r.status !== :fail || r.family === :performance, rep.results)
end
