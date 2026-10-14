"""
Geometry-backed fallbacks for GridInterface. Activates when GeometryOps is
loaded. Adds methods only to the hooks in `src/fallbacks.jl` (dispatching on the
manifold), never to a public verb, so provider methods can never be ambiguous
with these.

| hook | what it does |
|---|---|
| `fallback_treeify` / `fallback_subtree` | packed STR R-tree over `cell_extent` (caps become XYZ boxes) |
| `fallback_centroid` | planar: `GO.centroid` (or a point on surface); sphere: normalised vertex mean |
| `fallback_area` | `GO.area(manifold, cell_geometry)` (units of radius²) |
| `fallback_extent` | planar: bounding box; sphere: vertex cap (CAP_SLACK), whole sphere past π/2 |
| `fallback_contains` | half-space tests for `ConvexCells` on the sphere, winding number otherwise |
| `fallback_cellat` / `fallback_cellsat` | descend the cached tree by extent, interior test |
| `fallback_neighbors` | single ring: shared boundary points, CCW from edge (v1, v2); general: touching |
| `fallback_vertex_table` | merge corners once (quantised): ≥ 3 cells or boundary, open fans flagged |
| `fallback_query` | candidates from the tree, then an exact predicate (plane and sphere) |

Every fallback checks the spherical type contract (`UnitSphericalPoint` points,
`SphericalCap` extents) at its boundary.
"""
module GridInterfaceGeometryOpsExt

import GridInterface as GrI
using GridInterface: ncells, cellshape, cell_boundary, cell_geometry, cell_extent, cell_corners,
    cell_centroid, treeify, SingleRing, ConvexCells, GeneralCells, Edge, Vertex, Connectivity,
    check_point, check_extent, check_ring
const STI = GrI.STI
import GeometryOps as GO
import GeometryOps: FlexibleRTrees
using GeometryOps.UnitSpherical: UnitSphericalPoint, SphericalCap, spherical_distance
import GeometryOpsCore as GOC
using GeometryOpsCore: Planar, Spherical, Manifold
import Extents
import GeoInterface as GI

function __init__()
    GrI._HAS_BACKEND[] = true
end

# --- small vector helpers (no LinearAlgebra dependency) ----------------------
@inline _v(p::UnitSphericalPoint) = (Float64(p[1]), Float64(p[2]), Float64(p[3]))
@inline _v(p::NTuple{3,Real}) = (Float64(p[1]), Float64(p[2]), Float64(p[3]))
@inline _v(p) = (Float64(GI.x(p)), Float64(GI.y(p)), Float64(GI.z(p)))
@inline _dot(a, b) = a[1] * b[1] + a[2] * b[2] + a[3] * b[3]
@inline _cross(a, b) = (a[2] * b[3] - a[3] * b[2], a[3] * b[1] - a[1] * b[3], a[1] * b[2] - a[2] * b[1])
@inline _norm(a) = sqrt(_dot(a, a))
@inline _scale(a, s) = (a[1] * s, a[2] * s, a[3] * s)
@inline _sub(a, b) = (a[1] - b[1], a[2] - b[2], a[3] - b[3])
@inline _add(a, b) = (a[1] + b[1], a[2] + b[2], a[3] + b[3])
@inline _normalize(a) = _scale(a, inv(_norm(a)))
@inline _usp(a) = UnitSphericalPoint(a[1], a[2], a[3])
@inline _xy(p::NTuple{2,Real}) = (Float64(p[1]), Float64(p[2]))
@inline _xy(p) = (Float64(GI.x(p)), Float64(GI.y(p)))
# a point as a coordinate tuple on manifold m
@inline _c(::Planar, p) = _xy(p)
@inline _c(::Spherical, p) = _v(p)

# all points of a cell, unclosed ring when SingleRing; type-checked on the sphere
_cellpoints(grid, i) = check_ring(GOC.best_manifold(grid), _cellpoints(cellshape(grid), grid, i))
_cellpoints(::SingleRing, grid, i) = cell_boundary(grid, i)
_cellpoints(::GrI.CellShape, grid, i) = collect(GI.getpoint(cell_geometry(grid, i)))

const CAP_SLACK = 1.0001
_full_cap() = SphericalCap(UnitSphericalPoint(0.0, 0.0, 1.0), Float64(π))

function points_cap(pts)
    s = (0.0, 0.0, 0.0)
    for p in pts
        s = _add(s, _v(p))
    end
    ns = _norm(s)
    (ns > 0 && isfinite(ns)) || return _full_cap()
    c = _usp(_scale(s, inv(ns)))
    r = 0.0
    for p in pts
        r = max(r, spherical_distance(c, _usp(_v(p))))
    end
    r = nextfloat(r * CAP_SLACK)
    (isfinite(r) && r <= π / 2) || return _full_cap()
    return SphericalCap(c, r)
end

function points_box(pts)
    xmin = ymin = Inf; xmax = ymax = -Inf
    for p in pts
        x, y = _xy(p)
        xmin = min(xmin, x); xmax = max(xmax, x)
        ymin = min(ymin, y); ymax = max(ymax, y)
    end
    return Extents.Extent(X=(xmin, xmax), Y=(ymin, ymax))
end

# --- centroid, area, extent -------------------------------------------------
GrI.fallback_extent(::Planar, grid, i) = points_box(_cellpoints(grid, i))
GrI.fallback_extent(::Spherical, grid, i) = points_cap(_cellpoints(grid, i))

# The area centroid when it lies in the cell; otherwise (holes, multipolygons,
# strongly concave rings) a point on the surface: the midpoint of the widest
# interior span of a horizontal scanline through the largest part.
function GrI.fallback_centroid(m::Planar, grid, i)
    g = cell_geometry(grid, i)
    c = GO.centroid(g)
    p = (Float64(GI.x(c)), Float64(GI.y(c)))
    cellshape(grid) isa ConvexCells && return p
    GrI.fallback_contains(m, grid, i, p) && return p
    return _point_on_surface(g)
end

function _point_on_surface(g)
    poly = g
    if GI.trait(g) isa GI.MultiPolygonTrait
        polys = collect(GI.getgeom(g))
        poly = polys[argmax([abs(GO.area(GO.Planar(), q)) for q in polys])]
    end
    rings = [collect(GI.getpoint(r)) for r in GI.getring(poly)]
    ys = sort!(unique!([_xy(p)[2] for p in rings[1]]))
    k = argmax(diff(ys))
    y = (ys[k] + ys[k+1]) / 2
    xs = Float64[]
    for r in rings, j in 1:length(r)
        (x1, y1) = _xy(r[j]); (x2, y2) = _xy(r[mod1(j + 1, length(r))])
        (y1 > y) != (y2 > y) && push!(xs, x1 + (y - y1) * (x2 - x1) / (y2 - y1))
    end
    sort!(xs)
    best = 1
    for j in 1:2:length(xs)-1
        xs[j+1] - xs[j] > xs[best+1] - xs[best] && (best = j)
    end
    return ((xs[best] + xs[best+1]) / 2, y)
end
function GrI.fallback_centroid(::Spherical, grid, i)
    s = (0.0, 0.0, 0.0)
    for p in _cellpoints(grid, i)
        s = _add(s, _v(p))
    end
    return _usp(_normalize(s))
end

function GrI.fallback_area(m::Manifold, grid, i)
    g = cell_geometry(grid, i)
    GrI.iscurved(m) && (GI.trait(g) isa GI.PolygonTrait ? check_point(m, GI.getpoint(GI.getexterior(g), 1)) : nothing)
    return GO.area(m, g)
end

# --- interior tests ---------------------------------------------------------
# planar ring, closed (boundary counts): crossing number + on-segment check.
function _in_ring_planar(ring, q; tol=1e-12)
    x, y = q
    inside = false
    n = length(ring)
    n == 0 && return false
    last_ = _xy(ring[n]) == _xy(ring[1]) ? n - 1 : n
    for k in 1:last_
        x1, y1 = _xy(ring[k]); x2, y2 = _xy(ring[k == last_ ? 1 : k + 1])
        dx, dy = x2 - x1, y2 - y1
        cr = dx * (y - y1) - dy * (x - x1)
        if abs(cr) <= tol * max(abs(dx) + abs(dy), 1.0) &&
           min(x1, x2) - tol <= x <= max(x1, x2) + tol && min(y1, y2) - tol <= y <= max(y1, y2) + tol
            return true
        end
        if (y1 > y) != (y2 > y)
            xi = x1 + (y - y1) * dx / dy
            xi > x && (inside = !inside)
        end
    end
    return inside
end

function _on_ring_planar(ring, q; tol=1e-12)
    x, y = q
    n = length(ring)
    for k in 1:n
        x1, y1 = _xy(ring[k]); x2, y2 = _xy(ring[mod1(k + 1, n)])
        dx, dy = x2 - x1, y2 - y1
        abs(dx * (y - y1) - dy * (x - x1)) <= tol * max(abs(dx) + abs(dy), 1.0) &&
            min(x1, x2) - tol <= x <= max(x1, x2) + tol && min(y1, y2) - tol <= y <= max(y1, y2) + tol && return true
    end
    return false
end

# spherical: convex ring by half-spaces; otherwise the winding angle about p,
# which is +2π for a point inside a counter-clockwise ring and 0 or -2π (the
# antipodal side) outside, so only the sign-aware test is correct.
function _in_convex_sphere(ring, p; tol=1e-12)
    n = length(ring)
    for k in 1:n
        a = _v(ring[k]); b = _v(ring[mod1(k + 1, n)])
        nab = _cross(a, b)
        _dot(nab, p) < -tol * _norm(nab) && return false
    end
    return true
end
function _winding_sphere(ring, p; tol=1e-12)
    n = length(ring)
    total = 0.0
    for k in 1:n
        a = _v(ring[k]); b = _v(ring[mod1(k + 1, n)])
        _on_arc(a, b, p, tol) && return true
        _dot(a, p) >= 1 - 1e-15 && return true
        total += atan(_dot(p, _cross(a, b)), _dot(a, b) - _dot(a, p) * _dot(b, p))
    end
    return total > π
end
@inline function _on_arc(a, b, p, tol)
    nab = _cross(a, b)
    nn = _norm(nab)
    nn == 0 && return _norm(_sub(a, p)) <= tol
    return abs(_dot(nab, p)) <= tol * nn && _dot(_cross(a, p), nab) >= -tol * nn && _dot(_cross(p, b), nab) >= -tol * nn
end

GrI.fallback_ring_contains(::Planar, ring, p) = _in_ring_planar(ring, _xy(p))
GrI.fallback_ring_contains(::Spherical, ring, p) = _winding_sphere(ring, _v(p))

function GrI.fallback_contains(::Planar, grid, i, p)
    q = _xy(p)
    cellshape(grid) isa SingleRing && return _in_ring_planar(cell_boundary(grid, i), q)
    return _in_polys(Planar(), _polys(Planar(), cell_geometry(grid, i)), q)
end
function GrI.fallback_contains(m::Spherical, grid, i, p)
    q = _v(p)
    sh = cellshape(grid)
    sh isa ConvexCells && return _in_convex_sphere(check_ring(m, cell_boundary(grid, i)), q)
    sh isa SingleRing && return _winding_sphere(check_ring(m, cell_boundary(grid, i)), q)
    return _in_polys(m, _polys(m, cell_geometry(grid, i)), q)
end

# --- polygons as coordinate rings (the predicate kernel) ---------------------
# A polygon is a vector of unclosed rings of coordinate tuples (exterior first).
const Ring2 = Vector{NTuple{2,Float64}}
const Ring3 = Vector{NTuple{3,Float64}}
_ringtype(::Planar) = Ring2
_ringtype(::Spherical) = Ring3

function _unclosed(m, pts)
    r = _ringtype(m)(undef, 0)
    for p in pts
        GrI.iscurved(m) && check_point(m, p)
        push!(r, _c(m, p))
    end
    length(r) > 1 && r[1] == r[end] && pop!(r)
    return r
end
function _polys(m, g)
    t = GI.trait(g)
    t isa GI.MultiPolygonTrait && return [_poly(m, p) for p in GI.getgeom(g)]
    t isa GI.PolygonTrait && return [_poly(m, g)]
    throw(ArgumentError("expected a polygon or multipolygon, got $(typeof(t))"))
end
_poly(m, p) = [_unclosed(m, GI.getpoint(r)) for r in GI.getring(p)]
_cellpolys(m, grid, i) = cellshape(grid) isa SingleRing ?
    [[_unclosed(m, cell_boundary(grid, i))]] : _polys(m, cell_geometry(grid, i))

# closed containment in one ring
_in_ring(::Planar, r, q) = _in_ring_planar(r, q)
_in_ring(::Spherical, r, q) = _winding_sphere(r, q)
_on_ring(::Planar, r, q) = _on_ring_planar(r, q)
function _on_ring(::Spherical, r, q; tol=1e-12)
    n = length(r)
    for k in 1:n
        _on_arc(r[k], r[mod1(k + 1, n)], q, tol) && return true
    end
    return false
end
# closed containment in a polygon set; holes are open (their boundary is inside)
function _in_polys(m, polys, q)
    for poly in polys
        _in_ring(m, poly[1], q) || continue
        any(h -> _in_ring(m, h, q) && !_on_ring(m, h, q), @view poly[2:end]) && continue
        return true
    end
    return false
end
_on_polys(m, polys, q) = any(poly -> any(r -> _on_ring(m, r, q), poly), polys)
_strictly_in_polys(m, polys, q) = _in_polys(m, polys, q) && !_on_polys(m, polys, q)

# segment/arc intersection: :none, :touch (at an endpoint, collinear) or :cross (proper)
function _seg_relation(::Planar, a, b, c, d; tol=1e-12)
    o(p, q, r) = (q[1] - p[1]) * (r[2] - p[2]) - (q[2] - p[2]) * (r[1] - p[1])
    s = tol * max(1.0, abs(a[1]), abs(a[2]), abs(b[1]), abs(b[2]))
    d1 = o(c, d, a); d2 = o(c, d, b); d3 = o(a, b, c); d4 = o(a, b, d)
    if ((d1 > s && d2 < -s) || (d1 < -s && d2 > s)) && ((d3 > s && d4 < -s) || (d3 < -s && d4 > s))
        return :cross
    end
    on(p, q, r, dd) = abs(dd) <= s && min(p[1], q[1]) - s <= r[1] <= max(p[1], q[1]) + s &&
        min(p[2], q[2]) - s <= r[2] <= max(p[2], q[2]) + s
    (on(c, d, a, d1) || on(c, d, b, d2) || on(a, b, c, d3) || on(a, b, d, d4)) && return :touch
    return :none
end
function _seg_relation(::Spherical, a, b, c, d; tol=1e-12)
    n1 = _cross(a, b); n2 = _cross(c, d)
    (_norm(n1) == 0 || _norm(n2) == 0) && return :none
    t = _cross(n1, n2)
    if _norm(t) <= tol * _norm(n1) * _norm(n2)
        # same great circle: overlapping arcs touch
        (_on_arc(a, b, c, tol) || _on_arc(a, b, d, tol) || _on_arc(c, d, a, tol) || _on_arc(c, d, b, tol)) && return :touch
        return :none
    end
    t = _normalize(t)
    for x in (t, _scale(t, -1.0))
        _on_arc(a, b, x, tol) && _on_arc(c, d, x, tol) || continue
        endpoint = min(_norm(_sub(x, a)), _norm(_sub(x, b)), _norm(_sub(x, c)), _norm(_sub(x, d))) <= 1e-12
        return endpoint ? :touch : :cross
    end
    return :none
end

# worst relation between the edges of two polygon sets (:cross dominates :touch)
function _edge_relation(m, A, B)
    rel = :none
    for pa in A, ra in pa, pb in B, rb in pb
        na, nb = length(ra), length(rb)
        for i in 1:na, j in 1:nb
            r = _seg_relation(m, ra[i], ra[mod1(i + 1, na)], rb[j], rb[mod1(j + 1, nb)])
            r === :cross && return :cross
            r === :touch && (rel = :touch)
        end
    end
    return rel
end
_allpts(P) = (q for poly in P for r in poly for q in r)

_intersects(m, A, B) = _edge_relation(m, A, B) !== :none ||
    any(q -> _in_polys(m, B, q), _allpts(A)) || any(q -> _in_polys(m, A, q), _allpts(B))
# A within B: no proper crossing, every point of A in B, no point of B strictly inside A
_within(m, A, B) = _edge_relation(m, A, B) !== :cross &&
    all(q -> _in_polys(m, B, q), _allpts(A)) && !any(q -> _strictly_in_polys(m, A, q), _allpts(B))

# --- regions ---------------------------------------------------------------
# A lon/lat box on the sphere: its boundary sampled every ≤ 2° along parallels
# and meridians, joined by great-circle arcs (DGG's reading).
function _box_ring(::Spherical, ext)
    (x0, x1), (y0, y1) = ext.X, ext.Y
    nx = max(1, ceil(Int, (x1 - x0) / 2)); ny = max(1, ceil(Int, (y1 - y0) / 2))
    pts = NTuple{3,Float64}[]
    for k in 0:nx-1; push!(pts, _v(_ll(x0 + k * (x1 - x0) / nx, y0))); end
    for k in 0:ny-1; push!(pts, _v(_ll(x1, y0 + k * (y1 - y0) / ny))); end
    for k in 0:nx-1; push!(pts, _v(_ll(x1 - k * (x1 - x0) / nx, y1))); end
    for k in 0:ny-1; push!(pts, _v(_ll(x0, y1 - k * (y1 - y0) / ny))); end
    out = NTuple{3,Float64}[]
    for p in pts                                       # poles repeat: drop repeats
        (isempty(out) || _norm(_sub(p, out[end])) > 1e-14) && push!(out, p)
    end
    length(out) > 1 && _norm(_sub(out[1], out[end])) <= 1e-14 && pop!(out)
    return out
end
_box_ring(::Planar, ext) = [(ext.X[1], ext.Y[1]), (ext.X[2], ext.Y[1]), (ext.X[2], ext.Y[2]), (ext.X[1], ext.Y[2])]
_ll(x, y) = GO.UnitSphereFromGeographic()((Float64(x), Float64(y)))

_regionpolys(m, ext::Extents.Extent) = [[_box_ring(m, ext)]]
_regionpolys(m, g) = _polys(m, g)

GrI.fallback_region_extent(::Planar, g) = points_box(_allpts(_polys(Planar(), g)))
GrI.fallback_region_extent(m::Spherical, g) = points_cap(_allpts(_polys(m, g)))

# A lon/lat box as a covering cap: centre at the mean of its densified boundary,
# radius the largest boundary distance plus half the sample spacing.
function GrI.query_region(::Spherical, ext::Extents.Extent{(:X, :Y)})
    (x0, x1), (y0, y1) = ext.X, ext.Y
    n = 32
    pts = UnitSphericalPoint{Float64}[]
    for t in range(0, 1; length=n + 1)
        x = x0 + t * (x1 - x0); y = y0 + t * (y1 - y0)
        push!(pts, _ll(x, y0), _ll(x1, y), _ll(x, y1), _ll(x0, y))
    end
    s = (0.0, 0.0, 0.0)
    for p in pts
        s = _add(s, _v(p))
    end
    ns = _norm(s)
    ns > 1e-9 || return _full_cap()
    c = _usp(_scale(s, inv(ns)))
    r = maximum(p -> spherical_distance(c, p), pts)
    step = deg2rad(max(x1 - x0, y1 - y0)) / n
    r = r + step
    r < π || return _full_cap()
    return SphericalCap(c, r)
end

# --- queries ------------------------------------------------------------------
_candidates(m, grid, region) = _candidates_in(m, _cached_tree(grid), GrI.query_region(m, region))
function _candidates_in(m, tree, q)
    out = Int[]
    STI.depth_first_search(e -> _ext_intersects(e, q), tree) do i
        push!(out, i)
        nothing
    end
    return sort!(out)
end
@inline _ext_intersects(e, q) = Extents.intersects(e, q)

function GrI.fallback_query(m::Manifold, grid, pred::GrI.QueryPredicate)
    r = pred.region
    r isa SphericalCap && m isa Spherical && return _query_cap(m, grid, pred)
    R = _regionpolys(m, r)
    cands = _candidates(m, grid, r)
    return filter!(i -> _test(m, grid, i, pred, R), cands)
end
_test(m, grid, i, ::GrI.CentroidIn, R) = _in_polys(m, R, _c(m, cell_centroid(grid, i)))
function _test(m, grid, i, ::GrI.Intersects, R)
    _in_polys(m, R, _c(m, cell_centroid(grid, i))) && return true
    return _intersects(m, _cellpolys(m, grid, i), R)
end
_test(m, grid, i, ::GrI.Within, R) = _within(m, _cellpolys(m, grid, i), R)
_test(m, grid, i, ::GrI.Contains, R) = _within(m, R, _cellpolys(m, grid, i))

# A spherical cap region: distances, exact for geodesic edges (caps ≤ hemisphere are convex).
function _query_cap(m, grid, pred)
    cap = pred.region
    c = _v(cap.point); r = Float64(cap.radius)
    cands = _candidates_in(m, _cached_tree(grid), cap)
    return filter!(cands) do i
        if pred isa GrI.CentroidIn
            return _angle(c, _v(cell_centroid(grid, i))) <= r
        end
        P = _cellpolys(m, grid, i)
        if pred isa GrI.Within
            return all(q -> _angle(c, q) <= r, _allpts(P))
        elseif pred isa GrI.Intersects
            _in_polys(m, P, c) && return true
            return any(poly -> any(ring -> _ring_mindist(ring, c) <= r, poly), P)
        else    # Contains: the cap lies in the cell
            _in_polys(m, P, c) || return false
            return all(poly -> all(ring -> _ring_mindist(ring, c) >= r, poly), P)
        end
    end
end
_angle(a, b) = atan(_norm(_cross(a, b)), _dot(a, b))
function _ring_mindist(ring, c)
    n = length(ring); d = Inf
    for k in 1:n
        d = min(d, _arc_dist(ring[k], ring[mod1(k + 1, n)], c))
    end
    return d
end
# distance from c to the great-circle arc a-b
function _arc_dist(a, b, c)
    n = _cross(a, b); nn = _norm(n)
    nn == 0 && return _angle(a, c)
    n = _scale(n, inv(nn))
    proj = _sub(c, _scale(n, _dot(c, n)))
    if _norm(proj) > 0
        x = _normalize(proj)
        _on_arc(a, b, x, 1e-12) && return _angle(c, x)
    end
    return min(_angle(a, c), _angle(b, c))
end

# --- point extents ----------------------------------------------------------
_point_extent(::Planar, p) = (q = _xy(p); Extents.Extent(X=(q[1], q[1]), Y=(q[2], q[2])))
_point_extent(::Spherical, p) = SphericalCap(_usp(_v(p)), 0.0)
GrI.fallback_point_extent(m::Manifold, p) = _point_extent(m, p)

# --- R-tree fallback --------------------------------------------------------
_box(::Planar, e::Extents.Extent) = e
_box(m::Spherical, e) = convert(Extents.Extent, check_extent(m, e))

function _rtree(m, grid, cells)
    isempty(cells) && throw(ArgumentError("cannot build a tree over no cells"))
    exts = [_box(m, cell_extent(grid, i)) for i in cells]
    cap = GrI.leaf_capacity(grid)
    return FlexibleRTrees.RTree(FlexibleRTrees.STR(), Base.OneTo(ncells(grid));
        extents=exts, indices=collect(Int, cells), nodecapacity=cap === nothing ? 16 : Int(cap))
end
GrI.fallback_treeify(m::Manifold, grid) = _rtree(m, grid, 1:ncells(grid))
GrI.fallback_subtree(m::Manifold, grid, cells) = _rtree(m, grid, cells)

GrI.split_weight(t::FlexibleRTrees.RTree) = length(t.indices)
GrI.split_weight(n::FlexibleRTrees.RTreeNode) =
    n.tree.nodecapacity^(length(n.tree.levels) - 1 - n.level)

# --- per-grid cache for the fallbacks that query repeatedly ------------------
# Point location, neighbour and query fallbacks would otherwise rebuild an
# R-tree per call. Keyed by object identity (egal for isbits grids). A grid
# whose own tree has expensive node extents (CR's quadtree cursor on the
# sphere) gets an R-tree here instead. `GrI.clear_caches!()` empties it.
const _TREE_CACHE = IdDict{Any,Any}()
const _TREE_LOCK = ReentrantLock()
function _cached_tree(grid)
    grid isa GrI.GridAndTree && return grid.tree
    lock(_TREE_LOCK) do
        get!(() -> _lookup_tree(grid), _TREE_CACHE, grid)
    end
end
function _lookup_tree(grid)
    t = treeify(grid)
    # a lazily-extented provider cursor without analytic structure: prefer an R-tree
    if GrI.system(grid) === nothing && STI.node_extent_is_expensive(typeof(t)) && !(t isa FlexibleRTrees.RTree)
        return GrI.fallback_treeify(GOC.best_manifold(grid), grid)
    end
    return t
end

# --- point location (function barriers around the Any-typed cached tree) -----
GrI.fallback_cellat(m::Manifold, grid, p) = _cellat_in(m, grid, _cached_tree(grid), p)
function _cellat_in(m, grid, tree, p)
    q = _point_extent(m, p)
    best = Ref(typemax(Int))
    STI.depth_first_search(e -> _ext_intersects(e, q), tree) do i
        if i < best[] && GrI.cell_contains(grid, i, p)
            best[] = i
        end
        nothing
    end
    return best[] == typemax(Int) ? nothing : best[]
end
GrI.fallback_cellsat(m::Manifold, grid, p) = _cellsat_in(m, grid, _cached_tree(grid), p)
function _cellsat_in(m, grid, tree, p)
    q = _point_extent(m, p)
    out = Int[]
    STI.depth_first_search(e -> _ext_intersects(e, q), tree) do i
        GrI.cell_contains(grid, i, p) && push!(out, i)
        nothing
    end
    return sort!(out)
end

# --- geometric neighbours ---------------------------------------------------
_onseg(::Planar, a, b, q, tol) = begin
    a2, b2, q2 = _xy(a), _xy(b), _xy(q)
    dx, dy = b2[1] - a2[1], b2[2] - a2[2]
    L = abs(dx) + abs(dy)
    abs(dx * (q2[2] - a2[2]) - dy * (q2[1] - a2[1])) <= tol * max(L, 1e-300) &&
        min(a2[1], b2[1]) - tol <= q2[1] <= max(a2[1], b2[1]) + tol &&
        min(a2[2], b2[2]) - tol <= q2[2] <= max(a2[2], b2[2]) + tol
end
_onseg(::Spherical, a, b, q, tol) = begin
    a3, b3, q3 = _v(a), _v(b), _v(q)
    nab = _cross(a3, b3); nn = _norm(nab)
    nn == 0 && return _norm(_sub(a3, q3)) <= tol
    abs(_dot(nab, q3)) <= tol * nn && _dot(_cross(a3, q3), nab) >= -tol * nn && _dot(_cross(q3, b3), nab) >= -tol * nn
end
_same(::Planar, a, b, tol) = (x = _xy(a); y = _xy(b); abs(x[1] - y[1]) + abs(x[2] - y[2]) <= tol)
_same(::Spherical, a, b, tol) = _norm(_sub(_v(a), _v(b))) <= tol

function _on_boundary(m, ring, q, tol)
    n = length(ring)
    for k in 1:n
        _onseg(m, ring[k], ring[mod1(k + 1, n)], q, tol) && return true
    end
    return false
end

# number of distinct boundary points shared by two rings (2+ = an edge)
function _shared_points(m, ra, rb, tol)
    n = 0
    first_shared = nothing
    for q in ra
        _on_boundary(m, rb, q, tol) || continue
        if first_shared === nothing
            first_shared = q; n = 1
        elseif !_same(m, first_shared, q, tol)
            return 2
        end
    end
    for q in rb
        _on_boundary(m, ra, q, tol) || continue
        if first_shared === nothing
            first_shared = q; n = 1
        elseif !_same(m, first_shared, q, tol)
            return 2
        end
    end
    return n
end

_tol(::Planar, ext) = 1e-9 * max(ext.X[2] - ext.X[1], ext.Y[2] - ext.Y[1], 1e-300)
_tol(::Spherical, ext) = 1e-9

# azimuth of q about centre c, counter-clockwise seen from outside / from above
_azimuth(::Planar, c, q) = (a = _xy(c); b = _xy(q); atan(b[2] - a[2], b[1] - a[1]))
function _azimuth(::Spherical, c, q)
    n = _v(c)
    ref = abs(n[3]) < 0.9 ? (0.0, 0.0, 1.0) : (1.0, 0.0, 0.0)
    e1 = _normalize(_sub(ref, _scale(n, _dot(ref, n))))
    e2 = _cross(n, e1)
    v = _v(q)
    return atan(_dot(v, e2), _dot(v, e1))
end
_mid(::Planar, a, b) = ((_xy(a)[1] + _xy(b)[1]) / 2, (_xy(a)[2] + _xy(b)[2]) / 2)
_mid(::Spherical, a, b) = _normalize(_add(_v(a), _v(b)))

GrI.fallback_neighbors(m::Manifold, grid, i, conn::Connectivity) =
    _neighbors_in(m, grid, _cached_tree(grid), i, conn, cellshape(grid))

function _neighbors_in(m, grid, tree, i, conn, ::SingleRing)
    ra = check_ring(m, cell_boundary(grid, i))
    ext = cell_extent(grid, i)
    tol = _tol(m, m isa Planar ? ext : nothing)
    out = Int[]
    STI.depth_first_search(e -> _ext_intersects(e, ext), tree) do j
        j == i && return nothing
        k = _shared_points(m, ra, cell_boundary(grid, j), tol)
        (k == 0 || (conn isa Edge && k < 2)) || push!(out, j)
        nothing
    end
    # Counter-clockwise along cell i's ring, starting across the edge (corner 1,
    # corner 2): each neighbour sits at the middle of the stretch of i's boundary
    # it shares, measured in corner units (edge k spans [k-1, k)), so Edge
    # neighbour k lies across edge k. On non-conforming or curved-edge grids the
    # shared stretch is measured along the full (densified) boundary.
    rb = cell_boundary(grid, i)
    par = _corner_params(m, rb, cell_corners(grid, i))
    n = length(cell_corners(grid, i))
    # Ties (several cells meeting cell i at one point, e.g. around a pole) by
    # azimuth about i's centroid, measured from the side opposite that point.
    c = cell_centroid(grid, i)
    key = map(out) do j
        t = _shared_mid(m, rb, par, n, cell_boundary(grid, j), tol)
        q = _ring_point(m, rb, par, n, t)
        (round(mod(t - 0.5, n); digits=9), mod2pi(_azimuth(m, c, cell_centroid(grid, j)) - _azimuth(m, c, q) + π))
    end
    return out[sortperm(key)]
end

# corner-unit parameter of each boundary point of ring rb (corner k at k-1)
function _corner_params(m, rb, cs)
    nb, n = length(rb), length(cs)
    idx = zeros(Int, n)
    b = 1
    for k in 1:n, _ in 1:nb
        if _same(m, rb[b], cs[k], 1e-12)
            idx[k] = b; break
        end
        b = mod1(b + 1, nb)
    end
    par = zeros(nb)
    any(==(0), idx) && return [(b - 1) * n / nb for b in 1:nb]     # corners not found: uniform
    for k in 1:n
        b0 = idx[k]; b1 = idx[mod1(k + 1, n)]
        span = mod(b1 - b0, nb); span == 0 && (span = nb)
        for s in 0:span-1
            par[mod1(b0 + s, nb)] = (k - 1) + s / span
        end
    end
    return par
end
# parameter (corner units) of point q on ring rb, or nothing
function _ring_param(m, rb, par, n, q, tol)
    nb = length(rb)
    for b in 1:nb
        a = rb[b]; c = rb[mod1(b + 1, nb)]
        _same(m, a, q, tol) && return par[b]
        if _onseg(m, a, c, q, tol)
            f = _chord(m, a, q) / max(_chord(m, a, c), 1e-300)
            pb = par[b]; pc = par[mod1(b + 1, nb)]
            pc < pb && (pc += n)
            return pb + f * (pc - pb)
        end
    end
    return nothing
end
# the boundary point at corner-unit parameter t (nearest boundary point)
function _ring_point(m, rb, par, n, t)
    best = 1; d = Inf
    for b in eachindex(rb)
        e = abs(mod(par[b] - t + n / 2, n) - n / 2)
        e < d && ((best, d) = (b, e))
    end
    return rb[best]
end
_chord(::Planar, a, b) = (x = _xy(a); y = _xy(b); hypot(x[1] - y[1], x[2] - y[2]))
_chord(::Spherical, a, b) = _norm(_sub(_v(a), _v(b)))
function _shared_mid(m, rb, par, n, rj, tol)
    ts = Float64[]
    for q in rj
        t = _ring_param(m, rb, par, n, q, tol)
        t === nothing || push!(ts, mod(t, n))
    end
    for (b, q) in enumerate(rb)
        _on_boundary(m, rj, q, tol) && push!(ts, mod(par[b], n))
    end
    isempty(ts) && return 0.5
    sort!(ts)
    # the shared stretch is short: if it wraps past corner 1, unwrap it
    gaps = [ts[mod1(k + 1, length(ts))] - ts[k] + (k == length(ts) ? n : 0) for k in eachindex(ts)]
    k = argmax(gaps)                          # the largest gap is outside the stretch
    lo = ts[mod1(k + 1, length(ts))]; hi = ts[k]
    hi < lo && (hi += n)
    return mod((lo + hi) / 2, n)
end
# general cells: "touching" (shared boundary or overlapping area); unordered
function _neighbors_in(m, grid, tree, i, conn, ::GrI.CellShape)
    A = _cellpolys(m, grid, i)
    ext = cell_extent(grid, i)
    out = Int[]
    STI.depth_first_search(e -> _ext_intersects(e, ext), tree) do j
        j == i && return nothing
        _intersects(m, A, _cellpolys(m, grid, j)) && push!(out, j)
        nothing
    end
    return sort!(out)
end

# --- vertex identity --------------------------------------------------------
_qkey(::Planar, p, q) = (round(Int, GI.x(p) / q), round(Int, GI.y(p) / q), 0)
_qkey(::Spherical, p, q) = (v = _v(p); (round(Int, v[1] / q), round(Int, v[2] / q), round(Int, v[3] / q)))

function GrI.fallback_vertex_table(m::Manifold, grid)
    n = ncells(grid)
    q = 1e-9
    if m isa Planar
        e = GrI.grid_extent(grid)
        q *= max(e.X[2] - e.X[1], e.Y[2] - e.Y[1])
    end
    P = GrI._pointtype(grid)
    points = P[]
    index = Dict{NTuple{3,Int},Int}()
    cell_ptr = Vector{Int}(undef, n + 1)
    cell_verts = Int[]; cell_corner = Int[]
    cell_ptr[1] = 1
    for i in 1:n
        k0 = length(cell_verts) + 1
        for (k, p) in enumerate(check_ring(m, cell_corners(grid, i)))
            v = get!(index, _qkey(m, p, q)) do
                push!(points, p)
                length(points)
            end
            # a repeated corner (the pole of a lon/lat cell) is one vertex, listed once
            length(cell_verts) >= k0 && cell_verts[end] == v && continue
            push!(cell_verts, v); push!(cell_corner, k)
        end
        if length(cell_verts) > k0 && cell_verts[end] == cell_verts[k0]
            pop!(cell_verts); pop!(cell_corner)
        end
        cell_ptr[i+1] = length(cell_verts) + 1
    end
    nv = length(points)
    incid = [Int[] for _ in 1:nv]
    for i in 1:n, k in cell_ptr[i]:cell_ptr[i+1]-1
        push!(incid[cell_verts[k]], i)
    end
    _pn(cp, cv, c, v) = (r = cp[c]:cp[c+1]-1; k = findfirst(j -> cv[j] == v, r);
        L = length(r); (cv[r[mod1(k - 1, L)]], cv[r[mod1(k + 1, L)]]))
    # A point where only two cells meet, inside the grid (H3 class III bends,
    # densification points) is not a vertex: drop it and renumber.
    keep = trues(nv)
    for v in 1:nv
        length(incid[v]) == 2 || continue
        a, b = incid[v]
        pa = _pn(cell_ptr, cell_verts, a, v); pb = _pn(cell_ptr, cell_verts, b, v)
        (pa[1] == pb[2] && pb[1] == pa[2]) && (keep[v] = false)
    end
    if !all(keep)
        newid = cumsum(keep)
        points = points[keep]
        cv = Int[]; cc = Int[]; cp = similar(cell_ptr); cp[1] = 1
        for i in 1:n
            for k in cell_ptr[i]:cell_ptr[i+1]-1
                keep[cell_verts[k]] || continue
                push!(cv, newid[cell_verts[k]]); push!(cc, cell_corner[k])
            end
            cp[i+1] = length(cv) + 1
        end
        cell_verts, cell_corner, cell_ptr, incid = cv, cc, cp, incid[keep]
        nv = length(points)
    end
    vert_ptr = Vector{Int}(undef, nv + 1)
    vert_cells = Int[]
    isopen = falses(nv)
    vert_ptr[1] = 1
    for v in 1:nv
        cs = unique(incid[v])
        pn = [_pn(cell_ptr, cell_verts, c, v) for c in cs]
        # closed iff every cell's next edge is another cell's previous edge
        start = findfirst(k -> !any(q -> q[1] == pn[k][2], pn), eachindex(cs))
        if start === nothing && length(cs) >= 3
            az = [_azimuth(m, points[v], cell_centroid(grid, c)) for c in cs]
            append!(vert_cells, cs[sortperm(az)])
        else
            isopen[v] = true
            start === nothing && (start = 1)
            order = Int[start]
            while length(order) < length(cs)
                k = findfirst(q -> q[2] == pn[order[end]][1], pn)
                (k === nothing || k in order) && break
                push!(order, k)
            end
            length(order) < length(cs) && append!(order, setdiff(eachindex(cs), order))
            append!(vert_cells, cs[order])
        end
        vert_ptr[v+1] = length(vert_cells) + 1
    end
    return GrI.VertexTable{P}(points, cell_ptr, cell_verts, cell_corner, vert_ptr, vert_cells, isopen)
end

function GrI.clear_caches!()
    lock(_TREE_LOCK) do
        empty!(_TREE_CACHE)
    end
    lock(GrI._VERTEX_LOCK) do
        empty!(GrI._VERTEX_CACHE)
        empty!(GrI._NUMBERING_CACHE)
        GrI._LAST_VT[] = (nothing, nothing)
    end
    return nothing
end

end
