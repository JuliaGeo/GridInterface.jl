module StructuredDefaultsTests
using Test, GeometryOps, GridInterfaceConformance
import GridInterface as GrI
import GridInterface.Structured as S
import GeometryOpsCore as GOC
import GeoInterface as GI
import Extents

struct SharedDefaults{P} <: S.AbstractStructuredGrid
    points::P
end
GOC.best_manifold(::SharedDefaults) = GOC.Planar()
S.facesize(g::SharedDefaults, f) = (size(g.points, 1)-1, size(g.points, 2)-1)
S.facevertex(g::SharedDefaults, f, a, b) = g.points[a,b]

struct MultiDefaults <: S.AbstractStructuredGrid
    sizes::Vector{Tuple{Int,Int}}
end
GOC.best_manifold(::MultiDefaults) = GOC.Planar()
S.nfaces(g::MultiDefaults) = length(g.sizes)
S.facesize(g::MultiDefaults, f) = g.sizes[f]
S.facevertex(g::MultiDefaults, f, a, b) = (Float64(200*f+a), Float64(b))

struct ReversedDefaults{G} <: S.AbstractStructuredGrid
    grid::G
end
GOC.best_manifold(g::ReversedDefaults) = GOC.best_manifold(g.grid)
S.facesize(g::ReversedDefaults, f) = S.facesize(g.grid, f)
S.facevertex(g::ReversedDefaults, f, a, b) = S.facevertex(g.grid, f, a, b)
S.faceposition(g::ReversedDefaults, f, i, j) = GrI.ncells(g)+1-S.faceposition(g.grid, f, i, j)
S.facecell(g::ReversedDefaults, k) = S.facecell(g.grid, GrI.ncells(g)+1-k)
S.facecolumnmajor(::ReversedDefaults, f) = false

struct StoredDefaults{M,P} <: S.AbstractLatticeGrid
    manifold::M
    polygons::P
end
GOC.best_manifold(g::StoredDefaults) = g.manifold
S.facesize(g::StoredDefaults, f) = size(g.polygons)
GrI.cell_geometry(g::StoredDefaults, k::Integer) = g.polygons[k]

struct ForeignRaster{P}
    points::P
end
GOC.best_manifold(::ForeignRaster) = GOC.Planar()
S.facesize(g::ForeignRaster, f) = (size(g.points,1)-1, size(g.points,2)-1)
S.facevertex(g::ForeignRaster, f, a, b) = g.points[a,b]

@testset "optional face numbering and storage-only grids" begin
    shared = SharedDefaults([(Float64(i),Float64(j)) for i in 0:4,j in 0:3])
    multi = MultiDefaults([(2,3),(5,1),(1,4)])
    reversed = ReversedDefaults(shared)
    @test S.nfaces(shared) == 1
    @test GrI.ncells(shared) == 12 && GrI.ncells(multi) == 15
    @test S.facecell(multi,7) == (2,1,1)
    @test S.facecell(multi,12) == (3,1,1)
    for g in (shared,multi,reversed)
        @test all(k -> S.faceposition(g,S.facecell(g,k)...) == k, 1:GrI.ncells(g))
        @test @inferred(S.facecell(g,2)) isa NTuple{3,Int}
        @test @inferred(S.faceposition(g,1,2,1)) isa Int
        @test sort!(collect(GrI.subtree_indices(GrI.treeify(g)))) == collect(1:GrI.ncells(g))
        tiles = GrI.chunks(g,(2,2))
        @test sort!(vcat([collect(GrI.parentpositions(c)) for c in tiles]...)) == collect(1:GrI.ncells(g))
        @test all(c -> issorted(GrI.parentpositions(c)),tiles)
        @test all(k -> k in GrI.parentpositions(tiles[GrI.chunkat(tiles,k)]),1:GrI.ncells(g))
        @test_throws BoundsError GrI.cell_boundary(g,0)
        @test_throws BoundsError GrI.cell_boundary(g,GrI.ncells(g)+1)
        @test isempty([r for r in check_grid(g).results if r.status === :fail])
    end
    @test S.faceposition(reversed,1,1,1) == 12
    @test S.facecell(reversed,1) == (1,4,3)
    foreign = ForeignRaster(shared.points)
    @test !S.isstructured(foreign)
    @test_throws MethodError S.nfaces(foreign)
    @test_throws MethodError S.faceposition(foreign,1,1,1)

    outer = [(0.,0.),(4.,0.),(4.,4.),(0.,4.),(0.,0.)]
    hole = [(1.,1.),(1.,3.),(3.,3.),(3.,1.),(1.,1.)]
    p = GI.Polygon([outer,hole])
    stored = StoredDefaults(GOC.Planar(),fill(p,3,3))
    @test GrI.cellshape(stored) isa GrI.GeneralCells
    @test GrI.conformity(stored) isa GrI.NonConforming
    @test !GrI.neighbors_are_analytic(stored) && !GrI.vertices_are_analytic(stored) && !GrI.location_is_analytic(stored)
    @test GrI.maxneighbors(stored,GrI.Edge()) === nothing
    @test GrI.maxneighbors(stored,GrI.Vertex()) === nothing
    @test_throws ArgumentError GrI.cell_boundary(stored,1)
    @test_throws ArgumentError GrI.cell_corners(stored,1)
    @test GrI.cellat(stored,(0.5,0.5)) == 1
    @test GrI.cellat(stored,(2.,2.)) === nothing
    @test S.blockextent(stored,1,1:3,1:3) == GrI.cell_extent(stored,1)
    @test_throws ArgumentError S.blockextent(stored,1,1:0,1:3)
    @test_throws ArgumentError S.blockextent(stored,1,1:3,1:0)
    @test_throws ArgumentError S.blockextent(StoredDefaults(GOC.Spherical(),stored.polygons),1,1:3,1:3)
    # Matrix position does not imply geometric proximity to the block perimeter.
    stored.polygons[2,2] = GI.Polygon([[(100.,100.),(101.,100.),(101.,101.),(100.,101.),(100.,100.)]])
    expected = Extents.Extent(X=(0.,101.),Y=(0.,101.))
    @test S.blockextent(stored,1,1:3,1:3) == expected
    @test GrI.grid_extent(GrI.chunks(stored,(3,3))[1]) == expected
    @test sort!(collect(GrI.subtree_indices(GrI.treeify(stored)))) == collect(1:9)
end
end
