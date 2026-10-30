module StreamingPlanarExtentRegression
using Test
import GridInterface as GrI
import GeometryOps as GO
import GeometryOpsCore as GOC
import GeoInterface as GI
import Extents
using GeometryOps.UnitSpherical: UnitSphericalPoint, SphericalCap
const Ext = Base.get_extension(GrI, :GridInterfaceGeometryOpsExt)
@assert Ext !== nothing

# Same three-method GeneralCells contract as test/grids.jl's PolygonBag.
struct ExtentBag{M,P} <: GrI.AbstractGrid
    manifold::M
    geoms::P
end
GrI.ncells(g::ExtentBag) = length(g.geoms)
GOC.best_manifold(g::ExtentBag) = g.manifold
GrI.cell_geometry(g::ExtentBag, k::Integer) = g.geoms[k]
bag(geom; manifold=GOC.Planar()) = ExtentBag(manifold, [geom])
legacy_planar(g, k) = Ext.points_box(Ext._cellpoints(g, k))
function outcome(f)
    try
        e = f()
        return (:value, typeof(e), values(e))
    catch err
        return (:error, typeof(err), sprint(showerror, err))
    end
end
function matches_legacy(g, k=1)
    actual = outcome(() -> GrI.fallback_extent(GOC.Planar(), g, k))
    old = outcome(() -> legacy_planar(g, k))
    @test isequal(actual, old)
    return actual
end
ring(pts) = GI.LinearRing(pts)
poly(pts) = GI.Polygon([ring(pts)])

# GeoInterface's own test_primitives.jl implements empty polygons this way.
struct EmptyPolygon end
GI.isgeometry(::Type{EmptyPolygon}) = true
GI.geomtrait(::EmptyPolygon) = GI.PolygonTrait()
GI.ngeom(::GI.PolygonTrait, ::EmptyPolygon) = 0
GI.getgeom(::GI.PolygonTrait, ::EmptyPolygon, i) = nothing
GI.isempty(::GI.PolygonTrait, ::EmptyPolygon) = true

# A standard custom GeoInterface polygon forwarding its ring contract to a
# wrapped polygon, while furnishing a lazy, single-pass point iterator.
mutable struct PointVisits
    requests::Int
    starts::Int
    yielded::Int
end
struct OnePassPoints{P}
    points::P
    visits::PointVisits
    started::Base.RefValue{Bool}
end
Base.IteratorSize(::Type{<:OnePassPoints}) = Base.SizeUnknown()
Base.IteratorEltype(::Type{<:OnePassPoints}) = Base.HasEltype()
Base.eltype(::Type{OnePassPoints{P}}) where {P} = eltype(P)
function Base.iterate(it::OnePassPoints)
    it.started[] && error("point iterator restarted")
    it.started[] = true
    it.visits.starts += 1
    return iterate(it, 1)
end
function Base.iterate(it::OnePassPoints, k::Int)
    k > length(it.points) && return nothing
    it.visits.yielded += 1
    return (it.points[k], k + 1)
end
struct LazyPolygon{G,P}
    polygon::G
    points::P
    visits::PointVisits
end
GI.isgeometry(::Type{<:LazyPolygon}) = true
GI.geomtrait(::LazyPolygon) = GI.PolygonTrait()
GI.ncoord(::GI.PolygonTrait, g::LazyPolygon) = GI.ncoord(g.polygon)
GI.ngeom(::GI.PolygonTrait, g::LazyPolygon) = GI.ngeom(g.polygon)
GI.getgeom(::GI.PolygonTrait, g::LazyPolygon, i) = GI.getgeom(g.polygon, i)
function GI.getpoint(::GI.PolygonTrait, g::LazyPolygon)
    g.visits.requests += 1
    return OnePassPoints(g.points, g.visits, Ref(false))
end
function lazy_polygon(points)
    p = poly(points)
    visits = PointVisits(0, 0, 0)
    return LazyPolygon(p, points, visits)
end

# Some point sources reuse a mutable coordinate buffer. A streaming bound must
# cover emitted coordinates even though eager collect retains buffer aliases.
struct ReusedPoints{P}
    points::P
    buffer::Vector{Float64}
end
Base.IteratorSize(::Type{<:ReusedPoints}) = Base.SizeUnknown()
Base.IteratorEltype(::Type{<:ReusedPoints}) = Base.HasEltype()
Base.eltype(::Type{<:ReusedPoints}) = Vector{Float64}
function Base.iterate(it::ReusedPoints, k::Int=1)
    k > length(it.points) && return nothing
    it.buffer[1], it.buffer[2] = it.points[k]
    return (it.buffer, k + 1)
end
struct ReusingPolygon{G,P}
    polygon::G
    points::P
end
GI.isgeometry(::Type{<:ReusingPolygon}) = true
GI.geomtrait(::ReusingPolygon) = GI.PolygonTrait()
GI.ncoord(::GI.PolygonTrait, g::ReusingPolygon) = GI.ncoord(g.polygon)
GI.ngeom(::GI.PolygonTrait, g::ReusingPolygon) = GI.ngeom(g.polygon)
GI.getgeom(::GI.PolygonTrait, g::ReusingPolygon, i) = GI.getgeom(g.polygon, i)
GI.getpoint(::GI.PolygonTrait, g::ReusingPolygon) = ReusedPoints(g.points, zeros(2))

# Same checked SingleRing boundary contract as TinyPlanar and TupleLonLat.
struct RingGrid{M,P} <: GrI.AbstractGrid
    manifold::M
    points::P
end
GrI.ncells(::RingGrid) = 1
GOC.best_manifold(g::RingGrid) = g.manifold
GrI.cellshape(::RingGrid) = GrI.SimpleCells()
function GrI.cell_boundary(g::RingGrid, k::Integer)
    k == 1 || throw(BoundsError(g, k))
    return g.points
end

@testset "streaming planar GeneralCells extents" begin
    @testset "holes and separate multipolygon components" begin
        outer = [(0., 0.), (8., 0.), (8., 7.), (0., 7.), (0., 0.)]
        hole = [(1., 1.), (1., 6.), (6., 6.), (6., 1.), (1., 1.)]
        holey = GI.Polygon([ring(outer), ring(hole)])
        g = bag(holey)
        @test GrI.cellshape(g) isa GrI.GeneralCells
        @test GrI.cell_extent(g, 1) == Extents.Extent(X=(0., 8.), Y=(0., 7.))
        matches_legacy(g)
        remote = poly([(-9., -3.), (-7., -3.), (-7., -1.), (-9., -1.), (-9., -3.)])
        multi = GI.MultiPolygon([holey, remote])
        g = bag(multi)
        @test GrI.cell_extent(g, 1) == Extents.Extent(X=(-9., 8.), Y=(-3., 7.))
        matches_legacy(g)
        # Public geometry indexing must still reject invalid cell positions.
        @test_throws BoundsError GrI.cell_extent(g, 0)
        @test_throws BoundsError GrI.cell_extent(g, 2)
        matches_legacy(g, 0)
        matches_legacy(g, 2)
    end
    @testset "coordinate conversions and extent endpoint types" begin
        cases = (
            [(Int16(-4), Int16(2)), (Int16(5), Int16(2)), (Int16(5), Int16(8)), (Int16(-4), Int16(2))],
            [(Float32(-1.5), Float32(-2.25)), (Float32(3.5), Float32(-2.25)), (Float32(3.5), Float32(4.75)), (Float32(-1.5), Float32(-2.25))],
            [(-7//3, 2//5), (11//4, 2//5), (11//4, 13//6), (-7//3, 2//5)],
            [(big"-1.234567890123456789", big"0.25"), (big"8.76543210987654321", big"0.25"), (big"8.76543210987654321", big"3.125"), (big"-1.234567890123456789", big"0.25")],
            Tuple{Real,Real}[(-4, Float32(1.25)), (7//2, big"1.25"), (big"3.5", 9//2), (-4, Float32(1.25))],
        )
        for pts in cases
            g = bag(poly(pts))
            e = GrI.cell_extent(g, 1)
            xs = Float64[first(p) for p in pts]; ys = Float64[last(p) for p in pts]
            @test e == Extents.Extent(X=extrema(xs), Y=extrema(ys))
            @test e.X isa Tuple{Float64,Float64}
            @test e.Y isa Tuple{Float64,Float64}
            matches_legacy(g)
        end
    end
    @testset "empty GeoInterface polygon preserves baseline outcome" begin
        g = bag(EmptyPolygon())
        result = matches_legacy(g)
        if first(result) === :value
            @test isequal(result[3], ((Inf, -Inf), (Inf, -Inf)))
        end
    end
    @testset "lazy point iteration is requested and consumed once" begin
        pts = [(-3., 2.), (6., 2.), (6., 8.), (-3., 8.), (-3., 2.)]
        for extent_fn in ((g -> GrI.cell_extent(g, 1)), (g -> legacy_planar(g, 1)))
            lazy = lazy_polygon(pts)
            e = extent_fn(bag(lazy))
            @test e == Extents.Extent(X=(-3., 6.), Y=(2., 8.))
            @test lazy.visits.requests == 1
            @test lazy.visits.starts == 1
            @test lazy.visits.yielded == length(pts)
        end
    end
    @testset "reused mutable point buffer produces conservative bounds" begin
        pts = [(-3., 2.), (6., 2.), (6., 8.), (-3., 8.), (-3., 2.)]
        g = bag(ReusingPolygon(poly(pts), pts))
        expected = Extents.Extent(X=(-3., 6.), Y=(2., 8.))
        actual = GrI.cell_extent(g, 1)
        @test actual == expected
        @test all(p -> actual.X[1] <= p[1] <= actual.X[2] && actual.Y[1] <= p[2] <= actual.Y[2], pts)
    end
    @testset "SingleRing and curved validation stay intact" begin
        ringpts = [(0., 0.), (2., 0.), (2., 1.), (0., 1.)]
        planar_ring = RingGrid(GOC.Planar(), ringpts)
        @test GrI.cell_extent(planar_ring, 1) == Extents.Extent(X=(0., 2.), Y=(0., 1.))
        matches_legacy(planar_ring)
        badring = RingGrid(GOC.Spherical(), ringpts)
        @test_throws ArgumentError GrI.cell_extent(badring, 1)
        @test_throws ArgumentError GrI.fallback_extent(GOC.Planar(), badring, 1)
        matches_legacy(badring)
        badpoly = bag(poly(ringpts); manifold=GOC.Spherical())
        @test_throws ArgumentError GrI.cell_extent(badpoly, 1)
        @test_throws ArgumentError GrI.fallback_extent(GOC.Planar(), badpoly, 1)
        matches_legacy(badpoly)
        project = GrI.UnitSphereFromGeographic()
        spherepts = project.([(0., 0.), (2., 0.), (2., 1.), (0., 1.), (0., 0.)])
        goodpoly = bag(poly(spherepts); manifold=GOC.Spherical())
        e = GrI.cell_extent(goodpoly, 1)
        old = Ext.points_cap(Ext._cellpoints(goodpoly, 1))
        @test e isa SphericalCap
        @test e.point == old.point
        @test e.radius == old.radius
        matches_legacy(goodpoly)
    end
end
end
