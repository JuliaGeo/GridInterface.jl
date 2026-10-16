# ===========================================================================
# ConservativeRegridding's three curvilinear grids on the face hooks: a corner
# matrix, a pair of axis vectors, a polygon matrix. Each is one face whose
# positions are column-major, `i + (j - 1) * ni`.
# ===========================================================================

@inline _n1(g) = facesize(g, 1)[1]
@inline _colmajor_position(g, i, j) = i + (j - 1) * _n1(g)
@inline _colmajor_cell(g, k) = (n = _n1(g); (1, (k - 1) % n + 1, (k - 1) ÷ n + 1))

# The seams of one lattice face, declared at construction or detected there:
# an axis is periodic when its first and last vertex columns (rows) coincide,
# a row is collapsed (a pole) when all its vertices do. `facecross` wraps a
# periodic axis onto the opposite side and stops at a collapsed row, where
# `facefan` is every cell of the row (counter-clockwise: westward along the
# bottom row, eastward along the top).
function _seamcross(periodic, collapsed, ni, nj, i, j, s)
    if s == 1
        collapsed[1] && return nothing
        return periodic[2] ? (1, i, nj, 3) : nothing
    elseif s == 2
        return periodic[1] ? (1, 1, j, 4) : nothing
    elseif s == 3
        collapsed[2] && return nothing
        return periodic[2] ? (1, i, 1, 1) : nothing
    end
    return periodic[1] ? (1, ni, j, 2) : nothing
end
# a seam flag given by the caller, else detected
_seamflag(given, detect) = given === nothing ? detect() : NTuple{2,Bool}(given)
function _polefan(g, collapsed, b)
    ni, nj = facesize(g, 1)
    collapsed[1] && b == 1 && return collect(ni:-1:1)
    collapsed[2] && b == nj + 1 && return Int[(nj - 1) * ni + i for i in 1:ni]
    return nothing
end
function _facefan(g, periodic, collapsed, f, a, b)
    if collapsed[1] || collapsed[2]
        fan = _polefan(g, collapsed, b)
        fan === nothing || return fan
    end
    return _fanwalk(g, f, a, b)
end
_vertexmaxneighbors(collapsed) = (collapsed[1] || collapsed[2]) ? nothing : 8

"""
    CellBasedGrid(points::AbstractMatrix; periodic = nothing, collapsed = nothing)
    CellBasedGrid(manifold, points; periodic = nothing, collapsed = nothing)

One lattice face from an `(n1+1) × (n2+1)` matrix of corner points: cell
`(i, j)` has the corners `points[i, j], points[i+1, j], points[i+1, j+1],
points[i, j+1]`. A `UnitSphericalPoint` matrix is `Spherical()`, anything
else `Planar()`, unless the manifold is given.

`periodic` and `collapsed` (each an `NTuple{2,Bool}`, per lattice axis)
declare the seams: an axis is periodic when its last vertex column (row) is
its first again, and a row is collapsed to a point (a pole) when all its
vertices coincide. When `nothing` they are detected from the points with a
relative tolerance of `1e-12`, a default that reads a wrap column off by
more than that as an open border; a provider that knows its topology should
pass them.
"""
struct CellBasedGrid{M<:Manifold,P<:AbstractMatrix} <: AbstractStructuredGrid
    manifold::M
    points::P
    periodic::NTuple{2,Bool}
    collapsed::NTuple{2,Bool}
end
CellBasedGrid(points::AbstractMatrix; kw...) =
    CellBasedGrid(eltype(points) <: UnitSphericalPoint ? Spherical() : Planar(), points; kw...)
CellBasedGrid(manifold::Manifold, points::AbstractMatrix; periodic=nothing, collapsed=nothing) =
    CellBasedGrid(manifold, points, _seamflag(periodic, () -> _periodic(points)),
        _seamflag(collapsed, () -> _collapsed(points)))
function _periodic(points)
    ni, nj = size(points)
    return (all(b -> GrI._near(points[1, b], points[ni, b]), 1:nj),
            all(a -> GrI._near(points[a, 1], points[a, nj]), 1:ni))
end
function _collapsed(points)
    ni, nj = size(points)
    return (all(a -> GrI._near(points[a, 1], points[1, 1]), 2:ni),
            all(a -> GrI._near(points[a, nj], points[1, nj]), 2:ni))
end
GOC.best_manifold(g::CellBasedGrid) = g.manifold
nfaces(::CellBasedGrid) = 1
facesize(g::CellBasedGrid, f) = (size(g.points, 1) - 1, size(g.points, 2) - 1)
Base.@propagate_inbounds facevertex(g::CellBasedGrid, f, a, b) = g.points[a, b]
@inline faceposition(g::CellBasedGrid, f, i, j) = _colmajor_position(g, i, j)
@inline facecell(g::CellBasedGrid, k) = _colmajor_cell(g, k)
@inline facecross(g::CellBasedGrid, f, i, j, s) =
    ((ni, nj) = facesize(g, 1); _seamcross(g.periodic, g.collapsed, ni, nj, i, j, s))
facefan(g::CellBasedGrid, f, a, b) = _facefan(g, g.periodic, g.collapsed, f, a, b)
GrI.maxneighbors(g::CellBasedGrid, ::Vertex) = _vertexmaxneighbors(g.collapsed)

"""
    RegularGrid(x::AbstractVector, y::AbstractVector; periodic = nothing, collapsed = nothing)
    RegularGrid(manifold, x, y; periodic = nothing, collapsed = nothing)

One lattice face with separable axes: vertex `(a, b)` is `(x[a], y[b])` on
the plane (`Planar()`, the default) or the unit-sphere point of longitude
`x[a]` and latitude `y[b]` on `Spherical()`. Axes ascend. Point location is
the axis floor; cell extents are four reads on the plane
(`extent_is_expensive` is `false` there).

`periodic` and `collapsed` (each an `NTuple{2,Bool}`, per axis) declare the
seams: a periodic axis closes on itself (longitude spanning 360°), a
collapsed row is a pole. When `nothing` they are detected from the axes: on
a curved manifold `x` is periodic when it spans 360° to within `1e-9 × 360`,
and `y` collapses at an end reaching `±90`; a plane has neither. That
default reads a span off by more than the tolerance as an open border; a
provider that knows its topology should pass them.
"""
struct RegularGrid{M<:Manifold,DX,DY} <: AbstractStructuredGrid
    manifold::M
    x::DX
    y::DY
    periodic::NTuple{2,Bool}
    collapsed::NTuple{2,Bool}
end
RegularGrid(x, y; kw...) = RegularGrid(Planar(), x, y; kw...)
RegularGrid(manifold::Manifold, x, y; periodic=nothing, collapsed=nothing) =
    RegularGrid(manifold, x, y,
        _seamflag(periodic, () -> (_spansperiod(x, _axisperiod(manifold)[1]), _spansperiod(y, _axisperiod(manifold)[2]))),
        _seamflag(collapsed, () -> (iscurved(manifold) && first(y) <= -90, iscurved(manifold) && last(y) >= 90)))
# a lon axis spanning its period closes on itself
_spansperiod(x, period) = period !== nothing && abs((last(x) - first(x)) - period) <= 1e-9 * period
GOC.best_manifold(g::RegularGrid) = g.manifold
nfaces(::RegularGrid) = 1
facesize(g::RegularGrid, f) = (length(g.x) - 1, length(g.y) - 1)
Base.@propagate_inbounds facevertex(g::RegularGrid, f, a, b) = _vertex(g.manifold, g.x[a], g.y[b])
_vertex(::Planar, x, y) = (x, y)
_vertex(m, x, y) = iscurved(m) ? UnitSphereFromGeographic()((x, y)) : (x, y)
@inline faceposition(g::RegularGrid, f, i, j) = _colmajor_position(g, i, j)
@inline facecell(g::RegularGrid, k) = _colmajor_cell(g, k)
@inline facecross(g::RegularGrid, f, i, j, s) =
    ((ni, nj) = facesize(g, 1); _seamcross(g.periodic, g.collapsed, ni, nj, i, j, s))
facefan(g::RegularGrid, f, a, b) = _facefan(g, g.periodic, g.collapsed, f, a, b)
GrI.maxneighbors(g::RegularGrid, ::Vertex) = _vertexmaxneighbors(g.collapsed)
hasfaceaxes(::RegularGrid, f) = true
faceaxes(g::RegularGrid, f) = (g.x, g.y)
facecoords(g::RegularGrid, f, p) = iscurved(g.manifold) ? GeographicFromUnitSphere()(p) :
    (Float64(GI.x(p)), Float64(GI.y(p)))
# longitude is periodic whatever the grid's span
_axisperiod(m) = iscurved(m) ? (360.0, nothing) : (nothing, nothing)
faceperiod(g::RegularGrid, f) = _axisperiod(g.manifold)
GrI.extent_is_expensive(g::RegularGrid) = iscurved(g.manifold)
# a spherical lon/lat grid spanning 360° × 180° covers the sphere
GrI.isglobal(g::RegularGrid) = iscurved(g.manifold) &&
    last(g.x) - first(g.x) >= 360 && first(g.y) <= -90 && last(g.y) >= 90
function _blockextent(::Planar, g::RegularGrid, f, ir::UnitRange{Int}, jr::UnitRange{Int})
    return Extents.Extent(X=(g.x[first(ir)], g.x[last(ir)+1]), Y=(g.y[first(jr)], g.y[last(jr)+1]))
end
_cellextent(m::Planar, g::RegularGrid, f, i, j) = _blockextent(m, g, f, i:i, j:j)

"""
    ExplicitPolygonGrid(polygons::AbstractMatrix)
    ExplicitPolygonGrid(manifold, polygons)

One face of stored polygons (`Planar()` by default). The lattice gives it
positions, tiles and the quadtree; its cells are the polygons as given
(`GeneralCells`: [`cell_geometry`](@ref GridInterface.cell_geometry), no
lattice ring, no `facevertex`), with geometric fallbacks for everything else.
"""
struct ExplicitPolygonGrid{M<:Manifold,P<:AbstractMatrix} <: AbstractStructuredGrid
    manifold::M
    polygons::P
end
ExplicitPolygonGrid(polygons::AbstractMatrix) = ExplicitPolygonGrid(Planar(), polygons)
GOC.best_manifold(g::ExplicitPolygonGrid) = g.manifold
nfaces(::ExplicitPolygonGrid) = 1
facesize(g::ExplicitPolygonGrid, f) = size(g.polygons)
facevertex(g::ExplicitPolygonGrid, f, a, b) = _throw_facevertex(g)
@inline faceposition(g::ExplicitPolygonGrid, f, i, j) = _colmajor_position(g, i, j)
@inline facecell(g::ExplicitPolygonGrid, k) = _colmajor_cell(g, k)

GrI.cellshape(::ExplicitPolygonGrid) = GeneralCells()
GrI.cell_boundary(g::ExplicitPolygonGrid, k::Integer) = GrI._refuse_boundary(g)
GrI.cell_corners(g::ExplicitPolygonGrid, k::Integer) = GrI._refuse_boundary(g)
GrI.cell_geometry(g::ExplicitPolygonGrid, k::Integer) = ((_, i, j) = _facecell(g, k); g.polygons[i, j])
GrI.cell_extent(g::ExplicitPolygonGrid, k::Integer) =
    ((_, i, j) = _facecell(g, k); _blockextent(g.manifold, g, 1, i:i, j:j))
_cellextent(m, g::ExplicitPolygonGrid, f, i, j) = _blockextent(m, g, f, i:i, j:j)
GrI.neighbors_are_analytic(::ExplicitPolygonGrid) = false
GrI.winding(::ExplicitPolygonGrid, ::Connectivity) = Unordered()
GrI.maxneighbors(::ExplicitPolygonGrid, ::Edge) = nothing
GrI.neighbors(g::ExplicitPolygonGrid, k::Integer; connectivity::Connectivity=Vertex()) =
    GrI.fallback_neighbors(g.manifold, g, Int(k), connectivity)
GrI.vertices_are_analytic(::ExplicitPolygonGrid) = false
GrI.location_is_analytic(::ExplicitPolygonGrid) = false
GrI.cellat(g::ExplicitPolygonGrid, p) = GrI.fallback_cellat(g.manifold, g, check_point(g.manifold, p))

# CR's extents for a block of stored polygons: the union of their boxes on
# the plane; on the sphere the cap about the mean of every polygon vertex.
function _blockextent(::Planar, g::ExplicitPolygonGrid, f, ir::UnitRange{Int}, jr::UnitRange{Int})
    xmin = ymin = Inf; xmax = ymax = -Inf
    for j in jr, i in ir, pt in GI.getpoint(GI.getexterior(g.polygons[i, j]))
        x, y = Float64(GI.x(pt)), Float64(GI.y(pt))
        xmin = min(xmin, x); xmax = max(xmax, x)
        ymin = min(ymin, y); ymax = max(ymax, y)
    end
    return Extents.Extent(X=(xmin, xmax), Y=(ymin, ymax))
end
function _blockextent(::Spherical, g::ExplicitPolygonGrid, f, ir::UnitRange{Int}, jr::UnitRange{Int})
    all_points = UnitSphericalPoint{Float64}[]
    for j in jr, i in ir
        for pt in GI.getpoint(GI.getexterior(g.polygons[i, j]))
            push!(all_points, pt)
        end
    end
    isempty(all_points) && return SphericalCap(UnitSphericalPoint(0.0, 0.0, 1.0), 0.0)
    center = _normalize3(sum(all_points) / length(all_points))
    radius = maximum(p -> spherical_distance(center, p), all_points)
    return SphericalCap(center, radius * 1.0001)
end
