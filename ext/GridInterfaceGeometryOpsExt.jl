"""
Geometry-backed fallbacks for GridInterface. Activates when GeometryOps is
loaded. Adds methods only to the hooks in `src/fallbacks.jl` (dispatching on the
manifold), never to a public verb, so provider methods can never be ambiguous
with these.

| hook | what it does |
|---|---|
| `fallback_treeify` / `fallback_subtree` | packed STR R-tree over `cell_extent` (caps become XYZ boxes) |
| `fallback_centroid` | planar: `GO.centroid`; sphere: normalised vertex mean |
| `fallback_area` | `GO.area(manifold, cell_geometry)` |
| `fallback_extent` | planar: bounding box; sphere: vertex cap (CAP_SLACK), whole sphere past π/2 |
| `fallback_contains` | half-space tests for `ConvexCells` on the sphere, winding number otherwise |
| `fallback_cellat` | descend the (cached) tree by extent, interior test, lowest position wins |
| `fallback_neighbors` | tree query + shared-boundary classification (≥ 2 shared points = edge) |
| `fallback_vertex_table` | merge coincident boundary points once (quantised), CCW incidence |
"""
module GridInterfaceGeometryOpsExt

import GridInterface as GrI
using GridInterface: ncells, cellshape, cell_boundary, cell_geometry, cell_extent,
    cell_centroid, treeify, SingleRing, ConvexCells, GeneralCells, Edge, Vertex, Connectivity
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
@inline _v(p) = (Float64(p[1]), Float64(p[2]), Float64(p[3]))
@inline _dot(a, b) = a[1] * b[1] + a[2] * b[2] + a[3] * b[3]
@inline _cross(a, b) = (a[2] * b[3] - a[3] * b[2], a[3] * b[1] - a[1] * b[3], a[1] * b[2] - a[2] * b[1])
@inline _norm(a) = sqrt(_dot(a, a))
@inline _scale(a, s) = (a[1] * s, a[2] * s, a[3] * s)
@inline _sub(a, b) = (a[1] - b[1], a[2] - b[2], a[3] - b[3])
@inline _add(a, b) = (a[1] + b[1], a[2] + b[2], a[3] + b[3])
@inline _normalize(a) = _scale(a, inv(_norm(a)))
@inline _usp(a) = UnitSphericalPoint(a[1], a[2], a[3])
@inline _xy(p) = (Float64(GI.x(p)), Float64(GI.y(p)))

# all points of a cell, unclosed ring when SingleRing
_cellpoints(grid, i) = _cellpoints(cellshape(grid), grid, i)
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
    # scan halfway between the two most separated consecutive vertex heights
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

GrI.fallback_area(m::Manifold, grid, i) = GO.area(m, cell_geometry(grid, i))

# --- interior tests ---------------------------------------------------------
# planar ring, closed (boundary counts): crossing number + on-segment check.
function _in_ring_planar(ring, q; tol=1e-12)
    x, y = q
    inside = false
    n = length(ring)
    n == 0 && return false
    # tolerate explicitly closed rings
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

function _in_polygon_planar(poly, q)
    rings = collect(GI.getring(poly))
    _in_ring_planar(collect(GI.getpoint(rings[1])), q) || return false
    for h in Iterators.drop(rings, 1)
        hp = collect(GI.getpoint(h))
        # strictly inside a hole means outside; on the hole boundary counts as inside
        _in_ring_planar(hp, q) && !_on_ring_planar(hp, q) && return false
    end
    return true
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

function GrI.fallback_contains(::Planar, grid, i, p)
    q = _xy(p)
    cellshape(grid) isa SingleRing && return _in_ring_planar(cell_boundary(grid, i), q)
    g = cell_geometry(grid, i)
    if GI.trait(g) isa GI.MultiPolygonTrait
        return any(poly -> _in_polygon_planar(poly, q), GI.getgeom(g))
    end
    return _in_polygon_planar(g, q)
end

# spherical: convex ring by half-spaces; otherwise winding angle about p.
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
    # the tangent-plane winding sum is also ±2π about the antipode of a ring:
    # a cell smaller than a hemisphere holds p only on p's side of the sphere
    s = (0.0, 0.0, 0.0)
    for k in 1:n
        s = _add(s, _v(ring[k]))
    end
    _dot(s, p) > 0 || return false
    total = 0.0
    for k in 1:n
        a = _v(ring[k]); b = _v(ring[mod1(k + 1, n)])
        # on the arc a-b counts as inside
        nab = _cross(a, b)
        nn = _norm(nab)
        if nn > 0 && abs(_dot(nab, p)) <= tol * nn && _dot(_cross(a, p), nab) >= -tol && _dot(_cross(p, b), nab) >= -tol
            return true
        end
        _dot(a, p) >= 1 - 1e-15 && return true
        total += atan(_dot(p, _cross(a, b)), _dot(a, b) - _dot(a, p) * _dot(b, p))
    end
    return abs(total) > π
end
function GrI.fallback_contains(::Spherical, grid, i, p)
    q = _v(p)
    sh = cellshape(grid)
    sh isa ConvexCells && return _in_convex_sphere(cell_boundary(grid, i), q)
    sh isa SingleRing && return _winding_sphere(cell_boundary(grid, i), q)
    g = cell_geometry(grid, i)
    polys = GI.trait(g) isa GI.MultiPolygonTrait ? collect(GI.getgeom(g)) : [g]
    for poly in polys
        rings = collect(GI.getring(poly))
        ext = collect(GI.getpoint(rings[1]))
        _winding_sphere(_unclose(ext), q) || continue
        any(h -> _winding_sphere(_unclose(collect(GI.getpoint(h))), q), Iterators.drop(rings, 1)) && continue
        return true
    end
    return false
end
_unclose(pts) = (length(pts) > 1 && _v(pts[1]) == _v(pts[end])) ? pts[1:end-1] : pts

# --- point extents ----------------------------------------------------------
_point_extent(::Planar, p) = (q = _xy(p); Extents.Extent(X=(q[1], q[1]), Y=(q[2], q[2])))
_point_extent(::Spherical, p) = SphericalCap(_usp(_v(p)), 0.0)
GrI.fallback_point_extent(m::Manifold, p) = _point_extent(m, p)

# --- R-tree fallback --------------------------------------------------------
_box(::Planar, e::Extents.Extent) = e
_box(::Spherical, e::Extents.Extent) = e
_box(::Spherical, c::SphericalCap) = convert(Extents.Extent, c)

function _rtree(m, grid, cells)
    isempty(cells) && throw(ArgumentError("cannot build a tree over no cells"))
    exts = [_box(m, cell_extent(grid, i)) for i in cells]
    return FlexibleRTrees.RTree(FlexibleRTrees.STR(), Base.OneTo(ncells(grid));
        extents=exts, indices=collect(Int, cells))
end
GrI.fallback_treeify(m::Manifold, grid) = _rtree(m, grid, 1:ncells(grid))
GrI.fallback_subtree(m::Manifold, grid, cells) = _rtree(m, grid, cells)

GrI.split_weight(t::FlexibleRTrees.RTree) = length(t.indices)
GrI.split_weight(n::FlexibleRTrees.RTreeNode) =
    n.tree.nodecapacity^(length(n.tree.levels) - 1 - n.level)

# --- per-grid cache for the fallbacks that query repeatedly ------------------
# Point location and neighbour fallbacks would otherwise rebuild an R-tree per
# call. Keyed by object identity (egal for isbits grids). `GrI.clear_caches!()`.
const _TREE_CACHE = IdDict{Any,Any}()
const _TREE_LOCK = ReentrantLock()
function _cached_tree(grid)
    grid isa GrI.GridAndTree && return grid.tree
    lock(_TREE_LOCK) do
        get!(() -> treeify(grid), _TREE_CACHE, grid)
    end
end

# --- point location ---------------------------------------------------------
GrI.fallback_cellat(m::Manifold, grid, p) = _cellat_in(m, grid, _cached_tree(grid), p)
# function barrier: the cache is an IdDict{Any,Any}, so the tree's type is known only here
function _cellat_in(m, grid, tree, p)
    q = _point_extent(m, p)
    best = Ref(typemax(Int))
    STI.depth_first_search(e -> Extents.intersects(e, q), tree) do i
        if i < best[] && GrI.cell_contains(grid, i, p)
            best[] = i
        end
        nothing
    end
    return best[] == typemax(Int) ? nothing : best[]
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

function _shared_points(m, ra, rb, tol)
    shared = eltype(ra)[]
    for q in ra
        _on_boundary(m, rb, q, tol) && !any(s -> _same(m, s, q, tol), shared) && push!(shared, q)
    end
    for q in rb
        _on_boundary(m, ra, q, tol) && !any(s -> _same(m, s, q, tol), shared) && push!(shared, q)
    end
    return length(shared)
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

function GrI.fallback_neighbors(m::Manifold, grid, i, conn::Connectivity)
    ra = cell_boundary(grid, i)
    ext = cell_extent(grid, i)
    tol = _tol(m, m isa Planar ? ext : nothing)
    tree = _cached_tree(grid)
    cands = Int[]
    STI.depth_first_search(e -> Extents.intersects(e, ext), tree) do j
        j == i || push!(cands, j)
        nothing
    end
    out = Int[]
    for j in cands
        k = _shared_points(m, ra, cell_boundary(grid, j), tol)
        k == 0 && continue
        (conn isa Edge && k < 2) && continue
        push!(out, j)
    end
    c = cell_centroid(grid, i)
    az = [_azimuth(m, c, cell_centroid(grid, j)) for j in out]
    p = sortperm(collect(zip(az, out)))
    return out[p]
end

# --- vertex identity --------------------------------------------------------
_qkey(::Planar, p, q) = (round(Int, GI.x(p) / q), round(Int, GI.y(p) / q))
_qkey(::Spherical, p, q) = (v = _v(p); (round(Int, v[1] / q), round(Int, v[2] / q), round(Int, v[3] / q)))

function GrI.fallback_vertex_table(m::Manifold, grid)
    n = ncells(grid)
    q = 1e-9
    if m isa Planar
        e = GrI.grid_extent(grid)
        q *= max(e.X[2] - e.X[1], e.Y[2] - e.Y[1])
    end
    P = typeof(first(GrI.cell_corners(grid, 1)))
    points = P[]
    index = Dict{Any,Int}()
    cell_ptr = Vector{Int}(undef, n + 1)
    cell_verts = Int[]
    cell_ptr[1] = 1
    for i in 1:n
        k0 = length(cell_verts) + 1
        for p in GrI.cell_corners(grid, i)
            v = get!(index, _qkey(m, p, q)) do
                push!(points, p)
                length(points)
            end
            # a repeated corner (the pole of a lon/lat cell) is one vertex, listed once
            length(cell_verts) >= k0 && cell_verts[end] == v && continue
            push!(cell_verts, v)
        end
        length(cell_verts) > k0 && cell_verts[end] == cell_verts[k0] && pop!(cell_verts)
        cell_ptr[i+1] = length(cell_verts) + 1
    end
    nv = length(points)
    incid = [Int[] for _ in 1:nv]
    for i in 1:n, k in cell_ptr[i]:cell_ptr[i+1]-1
        push!(incid[cell_verts[k]], i)
    end
    # An interior point where only two cells meet (H3 class III edges bent at an
    # icosahedron edge) is not a vertex of the cell complex: drop it, renumber.
    keep = trues(nv)
    _pn(c, v) = (r = cell_ptr[c]:cell_ptr[c+1]-1; k = findfirst(j -> cell_verts[j] == v, r);
        L = length(r); (cell_verts[r[mod1(k - 1, L)]], cell_verts[r[mod1(k + 1, L)]]))
    for v in 1:nv
        length(incid[v]) == 2 || continue
        (a, b) = incid[v]; pa = _pn(a, v); pb = _pn(b, v)
        (pa[1] == pb[2] && pb[1] == pa[2]) && (keep[v] = false)
    end
    if !all(keep)
        newid = cumsum(keep)
        points = points[keep]
        cv = Int[]; cp = similar(cell_ptr); cp[1] = 1
        for i in 1:n
            for k in cell_ptr[i]:cell_ptr[i+1]-1
                keep[cell_verts[k]] && push!(cv, newid[cell_verts[k]])
            end
            cp[i+1] = length(cv) + 1
        end
        cell_verts, cell_ptr, incid = cv, cp, incid[keep]
        nv = length(points)
    end
    vert_ptr = Vector{Int}(undef, nv + 1)
    vert_cells = Int[]
    vert_ptr[1] = 1
    for v in 1:nv
        cs = unique(incid[v])
        az = [_azimuth(m, points[v], cell_centroid(grid, c)) for c in cs]
        append!(vert_cells, cs[sortperm(az)])
        vert_ptr[v+1] = length(vert_cells) + 1
    end
    return GrI.VertexTable{P}(points, cell_ptr, cell_verts, vert_ptr, vert_cells)
end

function GrI.clear_caches!()
    lock(_TREE_LOCK) do
        empty!(_TREE_CACHE)
    end
    lock(GrI._VERTEX_LOCK) do
        empty!(GrI._VERTEX_CACHE)
    end
    return nothing
end


# --- query regions ----------------------------------------------------------
# A lon/lat box as a covering cap: centre at the mean of its densified boundary,
# radius the largest boundary distance plus half the sample spacing (the far
# point of a box from an interior centre lies on its boundary).
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
_ll(x, y) = GO.UnitSphereFromGeographic()((Float64(x), Float64(y)))

end
