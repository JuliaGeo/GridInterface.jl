# ===========================================================================
# The face protocol and the verbs derived from it. The lattice conventions
# (corners, sides, the vertex a corner names) are fixed once, at
# `_lattice_corners` below.
# ===========================================================================

"""
    abstract type AbstractStructuredGrid <: AbstractGrid

Supertype of the structured family. A subtype implements [`nfaces`](@ref),
[`facesize`](@ref), [`facevertex`](@ref), [`faceposition`](@ref) and
[`facecell`](@ref); every other verb is derived here (and may be overridden).
Cells are `ConvexCells` by default (every lattice cell geodesically convex,
what ConservativeRegridding's clipper assumes); positions within a face must
be column-major (`i` fastest) for the tile chunks.
"""
abstract type AbstractStructuredGrid <: AbstractGrid end

"Whether `g` is a member of the structured family (an [`AbstractStructuredGrid`](@ref))."
isstructured(g) = g isa AbstractStructuredGrid

"""
    nfaces(grid) -> Int

Number of lattice faces (1 for a raster, 6 for a cubed sphere).
"""
function nfaces end

"""
    facesize(grid, f) -> (ni, nj)

Cells along each lattice axis of face `f`.
"""
function facesize end

"""
    facevertex(grid, f, a, b) -> point

Lattice vertex `(a, b)` of face `f`, `a in 1:ni+1`, `b in 1:nj+1`
(`UnitSphericalPoint`s on the sphere). Cell `(i, j)` has the corners
`(i, j), (i+1, j), (i+1, j+1), (i, j+1)`, counter-clockwise.
"""
function facevertex end

"""
    faceposition(grid, f, i, j) -> Int

Position of cell `(i, j)` of face `f`; a bijection onto `1:ncells(grid)`.
"""
function faceposition end

"""
    facecell(grid, p) -> (f, i, j)

The inverse of [`faceposition`](@ref).
"""
function facecell end

"""
    hasfaceaxes(grid, f) -> Bool

Whether the two coordinates of face `f` are separable (a raster, lon/lat):
the face answers [`faceaxes`](@ref), [`facecoords`](@ref),
[`faceperiod`](@ref) and point location is the axis floor. Default `false`.
"""
hasfaceaxes(grid, f) = false

"""
    faceaxes(grid, f) -> (xs, ys)

Ascending vertex coordinates along each axis of face `f`, of lengths `ni+1`
and `nj+1`.
"""
function faceaxes end

"""
    facecoords(grid, f, p) -> (x, y) | nothing

The chart coordinates of point `p` on face `f`, or `nothing` when `p` is not
on that face's chart.
"""
function facecoords end

"""
    faceperiod(grid, f) -> (px, py)

The period of each chart coordinate (`360` for longitude), or `nothing` for
a non-periodic one. Default `(nothing, nothing)`.
"""
faceperiod(grid, f) = (nothing, nothing)

"""
    facelocalindex(grid, f, x, y) -> (i, j) | nothing

The cell of face `f` whose axis intervals bracket chart coordinates `(x, y)`
(the floor on each axis, honouring the period; the last cell takes its
closing vertex), or `nothing` outside the axes.
"""
function facelocalindex(grid, f, x, y)
    xs, ys = faceaxes(grid, f)
    px, py = faceperiod(grid, f)
    i = _axisfloor(xs, x, px)
    i === nothing && return nothing
    j = _axisfloor(ys, y, py)
    j === nothing && return nothing
    return (i, j)
end
function _axisfloor(xs, v, period)
    n = length(xs)
    n < 2 && return nothing
    lo = first(xs)
    period === nothing || (v = lo + mod(v - lo, period))
    (lo <= v <= last(xs)) || return nothing
    return min(searchsortedlast(xs, v), n - 1)
end

"""
    facecolumnmajor(grid, f) -> Bool

Whether the positions of face `f` ascend with `i` fastest, then `j`
(`faceposition(grid, f, i, j) == faceposition(grid, f, 1, 1) + (i - 1) + (j - 1) * ni`),
so a lattice block's positions are strictly ascending in lattice order and
a tile chunk keeps its block in its type. Default `true` (the three concrete
grids); a face whose data order is something else (ring order, a
space-filling curve) answers `false` and its tiles are sorted subsets.
"""
facecolumnmajor(grid, f) = true

"""
    facecross(grid, f, i, j, side) -> (f2, i2, j2, side2) | nothing

The cell across lattice side `side` of cell `(i, j)` on face `f`, with the
side `side2` of that cell that is the same seam, when the side is a seam onto
another face or onto the same face (a periodic axis, a fold); `nothing` when
the side is a true boundary or a collapsed edge (a pole). Sides are numbered
as the corners are, `1:4` counter-clockwise, side `s` running from corner `s`
to corner `s+1`; the two cells run their shared side in opposite directions,
so corner `side` of `(i, j)` is corner `side2 + 1` of the landing cell and
corner `side + 1` is corner `side2` (a periodic wrap lands on the opposite
side: east onto west). Crossing back from `(f2, i2, j2, side2)` returns
`(f, i, j, side)`. Only called when `(i, j)` is on the border in that
direction; interior sides never reach it. Default `nothing`: a grid with no
seams declared is regional.
"""
facecross(grid, f, i, j, side) = nothing

"""
    facefan(grid, f, a, b) -> AbstractVector{Int}

Positions of the cells around lattice vertex `(a, b)` of face `f`,
counter-clockwise (an open fan starts just after its gap). Default: walk
around the vertex through [`facecross`](@ref) (at most 4 cells on a regular
vertex, 3 at a cube corner, stopping at a boundary for an open fan), a
`SmallVector{4,Int}`. A provider overrides it only where the walk cannot reach
every cell, such as a lon/lat pole where every cell of the polar row meets,
in which case the result may be a `Vector`.
"""
facefan(grid, f, a, b) = _fanwalk(grid, f, a, b)

# ---------------------------------------------------------------------------
# Lattice corners and the collapsed ring.
# ---------------------------------------------------------------------------

@noinline _throw_facevertex(g) = throw(ArgumentError(
    "$(typeof(g)) defines no `facevertex`: its cells are not a lattice of shared vertices"))

@inline function _facecell(g, k::Integer)
    1 <= k <= ncells(g) || GrI._throw_bounds(g, k)
    return facecell(g, Int(k))
end

# Cell (i, j) of face f has the lattice corners c1 = (i, j), c2 = (i+1, j),
# c3 = (i+1, j+1), c4 = (i, j+1), counter-clockwise, and the lattice sides
# s1 = (c1, c2) bottom, s2 = (c2, c3) right, s3 = (c3, c4) top, s4 = (c4, c1)
# left: side s runs from corner s to corner s+1 and faces the cell (i, j-1),
# (i+1, j), (i, j+1), (i-1, j) respectively. Two cells sharing a side run it
# in opposite directions (both are counter-clockwise), so when cell 2 lies
# across side s of cell 1 and cell 1 across side s2 of cell 2, corner s of
# cell 1 is corner s2+1 of cell 2 and corner s+1 of cell 1 is corner s2 of
# cell 2. A side whose two corners coincide (the top of a lon/lat pole cell)
# is degenerate: it is dropped from the ring and has no neighbour across it.
const _CORNER_OFFSET = ((0, 0), (1, 0), (1, 1), (0, 1))
@inline _corner_vertex(i, j, c) = (i + _CORNER_OFFSET[c][1], j + _CORNER_OFFSET[c][2])
@inline _lattice_corners(g, f, i, j) = (facevertex(g, f, i, j), facevertex(g, f, i + 1, j),
    facevertex(g, f, i + 1, j + 1), facevertex(g, f, i, j + 1))

"""
    Corners(c1, c2, c3, c4)

The counter-clockwise ring of a lattice cell with every degenerate side
dropped: corner `t` is the start of the `t`-th non-degenerate side
(`sides[t]`), so ring edge `t` is lattice side `sides[t]`. No heap storage.
"""
struct Corners{P} <: AbstractVector{P}
    pts::NTuple{4,P}
    sides::NTuple{4,Int}
    n::Int
end
function Corners(c::NTuple{4,P}) where {P}
    s = (0, 0, 0, 0); n = 0
    for k in 1:4
        GrI._near(c[k], c[mod1(k + 1, 4)]) && continue
        n += 1
        s = Base.setindex(s, k, n)
    end
    return Corners{P}(c, s, n)
end
Base.size(r::Corners) = (r.n,)
Base.IndexStyle(::Type{<:Corners}) = IndexLinear()
Base.@propagate_inbounds function Base.getindex(r::Corners, t::Int)
    @boundscheck 1 <= t <= r.n || throw(BoundsError(r, t))
    return r.pts[r.sides[t]]
end
# the ring corner equal to p, or 0
function _corner_of(r::Corners, p)
    for t in 1:r.n
        GrI._near(r[t], p) && return t
    end
    return 0
end

_corners(g, f, i, j) = Corners(_lattice_corners(g, f, i, j))
_corners(g, k::Integer) = ((f, i, j) = _facecell(g, k); _corners(g, f, i, j))

GrI.cell_boundary(g::AbstractStructuredGrid, k::Integer) = _corners(g, k)
GrI.cell_corners(g::AbstractStructuredGrid, k::Integer) = _corners(g, k)

# ---------------------------------------------------------------------------
# Identity, shape, extents.
# ---------------------------------------------------------------------------

function GrI.ncells(g::AbstractStructuredGrid)
    n = 0
    for f in 1:nfaces(g)
        ni, nj = facesize(g, f)
        n += ni * nj
    end
    return n
end
GrI.cellshape(::AbstractStructuredGrid) = ConvexCells()
GrI.extent_is_expensive(::AbstractStructuredGrid) = true
GrI.cell_extent(g::AbstractStructuredGrid, k::Integer) =
    ((f, i, j) = _facecell(g, k); blockextent(g, f, i:i, j:j))

# Lazy walk of the border vertices of the point block [a0..a1] × [b0..b1] of
# face f: west column, east column, then the interiors of the south and north
# rows (each corner once). CR's `CurvilinearGridPerimeterPoints` on the hooks.
struct Perimeter{G}
    g::G
    f::Int
    a0::Int
    a1::Int
    b0::Int
    b1::Int
end
Base.IteratorSize(::Type{<:Perimeter}) = Base.HasLength()
Base.length(p::Perimeter) = 2 * (p.a1 - p.a0 + 1) + 2 * (p.b1 - p.b0 + 1) - 4
@inline function Base.iterate(p::Perimeter, n::Int=1)
    n > length(p) && return nothing
    H = p.b1 - p.b0 + 1
    W = p.a1 - p.a0 + 1
    if n <= H
        a, b = p.a0, p.b0 + n - 1
    elseif n <= 2H
        a, b = p.a1, p.b0 + (n - H) - 1
    elseif n <= 2H + (W - 2)
        a, b = p.a0 + (n - 2H), p.b0
    else
        a, b = p.a0 + (n - 2H - (W - 2)), p.b1
    end
    return facevertex(p.g, p.f, a, b), n + 1
end

# The extent of the cell block ir × jr of face f: a box on the plane, a cap on
# the sphere, from the block's perimeter vertices (a lattice block's interior
# vertices lie inside its perimeter).
function _blockextent(::Planar, g, f, ir::UnitRange{Int}, jr::UnitRange{Int})
    xmin = ymin = Inf; xmax = ymax = -Inf
    for p in Perimeter(g, f, first(ir), last(ir) + 1, first(jr), last(jr) + 1)
        x, y = Float64(GI.x(p)), Float64(GI.y(p))
        xmin = min(xmin, x); xmax = max(xmax, x); ymin = min(ymin, y); ymax = max(ymax, y)
    end
    return Extents.Extent(X=(xmin, xmax), Y=(ymin, ymax))
end
function _blockextent(::Spherical, g, f, ir::UnitRange{Int}, jr::UnitRange{Int})
    a0, a1 = first(ir), last(ir) + 1
    b0, b1 = first(jr), last(jr) + 1
    return _spherical_cap(facevertex(g, f, a0, b0), facevertex(g, f, a1, b0),
        facevertex(g, f, a1, b1), facevertex(g, f, a0, b1), Perimeter(g, f, a0, a1, b0, b1))
end
_blockextent(m, g, f, ir, jr) = _blockextent_other(m, g, f, ir, jr)
_blockextent_other(m, g, f, ir, jr) = iscurved(m) ? _blockextent(Spherical(), g, f, ir, jr) :
    _blockextent(Planar(), g, f, ir, jr)

"""
    blockextent(grid, f, ir::UnitRange{Int}, jr::UnitRange{Int})

The extent of the cell block `ir × jr` of face `f`: an `Extents.Extent` box
on the plane, a `SphericalCap` on a curved manifold, containing every cell of
the block. The one override point for a provider whose block extents have a
closed form (a geographic raster's analytic cap): the quadtree's
`node_extent`, its leaf entries (`ir` and `jr` of length one),
[`cell_extent`](@ref GridInterface.cell_extent) and a tile's `grid_extent`
all read it. The default derives a block's extent from its perimeter vertices
and a single cell's from its four corners (the private `_blockextent` and
`_cellextent`, which a provider may also specialize per manifold).
"""
blockextent(g, f, ir::UnitRange{Int}, jr::UnitRange{Int}) =
    (length(ir) == 1 && length(jr) == 1) ? _cellextent(best_manifold(g), g, f, first(ir), first(jr)) :
    _blockextent(best_manifold(g), g, f, ir, jr)

# a single cell: its four corners, no perimeter
function _cellextent(m, g, f, i, j)
    c = _lattice_corners(g, f, i, j)
    iscurved(m) && return _spherical_cap(c[1], c[2], c[3], c[4], ())
    xmin = ymin = Inf; xmax = ymax = -Inf
    for p in c
        x, y = Float64(GI.x(p)), Float64(GI.y(p))
        xmin = min(xmin, x); xmax = max(xmax, x); ymin = min(ymin, y); ymax = max(ymax, y)
    end
    return Extents.Extent(X=(xmin, xmax), Y=(ymin, ymax))
end

# --- the bounding cap, as CR computes it (bit for bit) --------------------------

@inline function _whole_sphere_cap(::Type{T}) where {T<:AbstractFloat}
    return SphericalCap(UnitSphericalPoint(zero(T), zero(T), one(T)), nextfloat(T(pi)))
end
const _WHOLE_SPHERE_CAP32 = _whole_sphere_cap(Float32)
const _WHOLE_SPHERE_CAP64 = _whole_sphere_cap(Float64)
@inline _whole_sphere_cap(::Type{Float32}) = _WHOLE_SPHERE_CAP32
@inline _whole_sphere_cap(::Type{Float64}) = _WHOLE_SPHERE_CAP64

# StaticArrays' `norm` and `normalize` (`inv(norm(a)) * a`) spelled out, so the
# caps match CR's without a LinearAlgebra dependency.
@inline _norm3(p) = sqrt(abs2(p[1]) + abs2(p[2]) + abs2(p[3]))
@inline function _normalize3(p::UnitSphericalPoint)
    s = inv(_norm3(p))
    return UnitSphericalPoint(s * p[1], s * p[2], s * p[3])
end

# A cap covering four CCW corners (SW, SE, NE, NW) and every point of
# `perimeter`. A cap no wider than a hemisphere is geodesically convex, so
# covering the vertices covers the geodesic edges between them; wider
# candidates become the whole sphere.
@inline function _spherical_cap(p1, p2, p3, p4, perimeter)
    meanpoint = (p1 + p2 + p3 + p4) / 4
    T = eltype(meanpoint)
    mean_norm = _norm3(meanpoint)
    (!isfinite(mean_norm) || mean_norm <= eps(T)) && return _whole_sphere_cap(T)
    center = _normalize3(meanpoint)
    (isfinite(center[1]) && isfinite(center[2]) && isfinite(center[3])) || return _whole_sphere_cap(T)
    d = zero(T)
    for p in (p1, p2, p3, p4)
        distance = spherical_distance(center, p)
        !isfinite(distance) && return _whole_sphere_cap(T)
        d = max(d, distance)
    end
    for p in perimeter
        distance = spherical_distance(center, p)
        !isfinite(distance) && return _whole_sphere_cap(T)
        d = max(d, distance)
    end
    radius = nextfloat(T(d * T(1.0001)))
    (!isfinite(radius) || radius > T(pi) / T(2)) && return _whole_sphere_cap(T)
    return SphericalCap(center, radius)
end

"""
    circle_from_four_corners(corner_points, other_points) -> SphericalCap

Bounding cap of four cell corners given in (BL, TL, BR, TR) order plus an
iterable of further perimeter points; `UnitSphericalPoint`s or `(lon, lat)`
tuples. ConservativeRegridding's public helper, kept for its extensions.
"""
function circle_from_four_corners(corner_points, other_points)
    raw = UnitSphereFromGeographic().(corner_points)
    p1, p2, p3, p4 = raw[1], raw[3], raw[4], raw[2]
    return _spherical_cap(p1, p2, p3, p4, (UnitSphereFromGeographic()(p) for p in other_points))
end

# ---------------------------------------------------------------------------
# The tree family: the per-face quadtree (quadtree.jl) and tiles (tiles.jl).
# ---------------------------------------------------------------------------

GrI.leaf_capacity(::AbstractStructuredGrid, ::ForSweeps) = 16     # 4x4 leaves
GrI.leaf_capacity(::AbstractStructuredGrid, ::ForQueries) = 4     # 2x2: CR's default

# a cell count per leaf -> the quadtree's per-side leaf size
_leafsize(n::Integer) = (s = max(1, isqrt(Int(n))); (s, s))

GrI.treeify(g::AbstractStructuredGrid, n::Integer) = nfaces(g) == 1 ?
    TopDownQuadtreeCursor(g, 1; leafsize=_leafsize(n)) : FaceRootTree(g, _leafsize(n))

# ---------------------------------------------------------------------------
# Adjacency on the lattice: the cell across a side, within the face by index
# arithmetic, across a face border by `facecross`.
# ---------------------------------------------------------------------------

# (f2, i2, j2) across side s, or nothing
@inline function _acrosscell(g, f, i, j, s)
    x = _cross(g, f, i, j, s)
    x === nothing && return nothing
    return (x[1], x[2], x[3])
end
@inline function _facecross(g, f, i, j, s)
    x = facecross(g, f, i, j, s)
    x === nothing && return nothing
    return (Int(x[1]), Int(x[2]), Int(x[3]), Int(x[4]))
end

# the position across side s, or 0
@inline function _across(g, f, i, j, s)
    x = _acrosscell(g, f, i, j, s)
    x === nothing && return 0
    return faceposition(g, x[1], x[2], x[3])
end

# The crossing with the side it lands on: (f2, i2, j2, s2), or nothing. Within
# the face the landing side is the opposite one; across a seam `facecross`
# reports it.
@inline function _cross(g, f, i, j, s)
    ni, nj = facesize(g, f)
    s == 1 && return j > 1 ? (f, i, j - 1, 3) : _facecross(g, f, i, j, 1)
    s == 2 && return i < ni ? (f, i + 1, j, 4) : _facecross(g, f, i, j, 2)
    s == 3 && return j < nj ? (f, i, j + 1, 1) : _facecross(g, f, i, j, 3)
    return i > 1 ? (f, i - 1, j, 2) : _facecross(g, f, i, j, 4)
end
@noinline _throw_fan(g, f, a, b) = throw(ArgumentError(
    "$(typeof(g)): more than 4 cells meet at lattice vertex ($a, $b) of face $f: override `facefan`"))

# The default `facefan`: the cells around lattice vertex (a, b), walking from
# the first in-face cell at the vertex. Counter-clockwise about a corner the
# next cell lies across the side ending at it (side t-1 of corner t), and
# the vertex is corner s2 of the cell landed on; clockwise it lies across the
# side starting at it (side t), and the vertex is corner s2+1.
function _fanwalk(g, f, a, b)
    ni, nj = facesize(g, f)
    if a <= ni && b <= nj
        i, j, t = a, b, 1
    elseif a > 1 && b <= nj
        i, j, t = a - 1, b, 2
    elseif a > 1 && b > 1
        i, j, t = a - 1, b - 1, 3
    else
        i, j, t = a, b - 1, 4
    end
    k0 = faceposition(g, f, i, j)
    ccw = push(SmallVector{4,Int}(), k0)
    f1, i1, j1, t1 = f, i, j, t
    while true
        x = _cross(g, f1, i1, j1, mod1(t1 - 1, 4))
        x === nothing && break
        f1, i1, j1, s2 = x
        k = faceposition(g, f1, i1, j1)
        k == k0 && return ccw                      # closed
        length(ccw) < 4 || _throw_fan(g, f, a, b)
        ccw = push(ccw, k); t1 = s2
    end
    # open: the cells clockwise of k0 come first, nearest the gap first
    cw = SmallVector{4,Int}()
    f1, i1, j1, t1 = f, i, j, t
    while true
        x = _cross(g, f1, i1, j1, t1)
        x === nothing && break
        f1, i1, j1, s2 = x
        k = faceposition(g, f1, i1, j1)
        k == k0 && break
        length(cw) + length(ccw) < 4 || _throw_fan(g, f, a, b)
        cw = push(cw, k); t1 = mod1(s2 + 1, 4)
    end
    out = SmallVector{4,Int}()
    for q in length(cw):-1:1
        out = push(out, cw[q])
    end
    for q in 1:length(ccw)
        out = push(out, ccw[q])
    end
    return out
end

# ---------------------------------------------------------------------------
# Neighbours: counter-clockwise from the neighbour across side 1, so Edge
# neighbour t lies across ring edge t; Vertex neighbours add, after the
# neighbour across side t, the cells of corner t+1's fan strictly between it
# and the neighbour across side t+1.
# ---------------------------------------------------------------------------

GrI.neighbors_are_analytic(::AbstractStructuredGrid) = true
GrI.winding(::AbstractStructuredGrid, ::Connectivity) = CounterClockwise()
GrI.maxneighbors(::AbstractStructuredGrid, ::Edge) = 4
GrI.maxneighbors(::AbstractStructuredGrid, ::Vertex) = 8

# x appended unless already present; SmallVector stays on the stack. A Vector
# (a pole fan of thousands) is deduplicated once at the end instead, keeping
# first occurrences as the scan would.
@inline function _add(s::SmallVector{N,Int}, x::Int) where {N}
    for q in 1:length(s)
        s[q] == x && return s
    end
    return push(s, x)
end
@inline _add(v::Vector{Int}, x::Int) = (push!(v, x); v)
@inline _finish(s::SmallVector) = s
@inline _finish(v::Vector{Int}) = unique!(v)

GrI.neighbors(g::AbstractStructuredGrid, k::Integer; connectivity::Connectivity=Vertex()) =
    _structured_neighbors(g, Int(k), connectivity)
# the derivation on its own, for a provider that post-processes it
function _structured_neighbors(g, k::Int, connectivity::Connectivity)
    f, i, j = _facecell(g, k)
    across = (_across(g, f, i, j, 1), _across(g, f, i, j, 2), _across(g, f, i, j, 3), _across(g, f, i, j, 4))
    if connectivity isa Edge
        out = SmallVector{4,Int}()
        for t in 1:4
            q = across[t]
            (q == 0 || q == k) || (out = _add(out, q))
        end
        return out
    end
    # the fans at corners 2, 3, 4, 1 (corner t+1 after side t); each checked
    # on its own so the stack path sees a concrete tuple
    F2 = _cornerfan(g, f, i, j, 2); F3 = _cornerfan(g, f, i, j, 3)
    F4 = _cornerfan(g, f, i, j, 4); F1 = _cornerfan(g, f, i, j, 1)
    if F2 isa SmallVector{4,Int} && F3 isa SmallVector{4,Int} && F4 isa SmallVector{4,Int} && F1 isa SmallVector{4,Int}
        return _vertex_neighbors(SmallVector{8,Int}(), k, across, (F2, F3, F4, F1))
    end
    return _vertex_neighbors(Int[], k, across, (F2, F3, F4, F1))
end
@inline _cornerfan(g, f, i, j, c) = ((a, b) = _corner_vertex(i, j, c); facefan(g, f, a, b))

function _vertex_neighbors(out, k, across, fans)
    for t in 1:4
        A = across[t]; B = across[mod1(t + 1, 4)]
        (A == 0 || A == k) || (out = _add(out, A))
        F = fans[t]
        m = length(F)
        pos = findfirst(==(k), F)
        pos === nothing && _throw_fan_owner(k)
        # the fan from k, counter-clockwise: A first and B last when the sides are not open
        for q in 1:m-1
            c = F[mod1(pos + q, m)]
            (q == 1 && c == A) && continue
            (q == m - 1 && c == B) && continue
            c == k || (out = _add(out, c))
        end
    end
    return _finish(out)
end
@noinline _throw_fan_owner(k) = throw(ArgumentError("`facefan` at a corner of cell $k does not contain it"))

# ---------------------------------------------------------------------------
# Vertices: every lattice vertex is a vertex; the key is its corner on the
# lowest cell of its fan.
# ---------------------------------------------------------------------------

GrI.vertices_are_analytic(::AbstractStructuredGrid) = true

function GrI.vertex_key(g::AbstractStructuredGrid, k::Integer, c::Integer)
    k = Int(k); c = Int(c)
    f, i, j = _facecell(g, k)
    r = _corners(g, f, i, j)
    1 <= c <= r.n || throw(BoundsError(r, c))
    a, b = _corner_vertex(i, j, r.sides[c])
    owner = minimum(facefan(g, f, a, b))
    owner == k && return VertexKey(k, c)
    t = _corner_of(_corners(g, owner), r[c])
    t == 0 && error("vertex of cell $k is no corner of its owner $owner")
    return VertexKey(owner, t)
end

function GrI.vertex_cells(g::AbstractStructuredGrid, key::VertexKey)
    f, i, j = _facecell(g, key.cell)
    r = _corners(g, f, i, j)
    1 <= key.corner <= r.n || throw(BoundsError(r, key.corner))
    a, b = _corner_vertex(i, j, r.sides[key.corner])
    return facefan(g, f, a, b)
end

# ---------------------------------------------------------------------------
# Point location: the axis floor on faces with axes (repaired against the
# geodesic edges on a curved manifold), a greedy walk over the lattice
# otherwise, with the tree as the last resort.
# ---------------------------------------------------------------------------

function GrI.location_is_analytic(g::AbstractStructuredGrid)
    for f in 1:nfaces(g)
        hasfaceaxes(g, f) || return false
    end
    return true
end

function GrI.cellat(g::AbstractStructuredGrid, p)
    m = best_manifold(g)
    p = check_point(m, p)
    allaxes = true
    for f in 1:nfaces(g)
        if !hasfaceaxes(g, f)
            allaxes = false
            continue
        end
        xy = facecoords(g, f, p)
        xy === nothing && continue
        ij = facelocalindex(g, f, xy[1], xy[2])
        ij === nothing && continue
        return _repair(m, g, f, ij[1], ij[2], p)
    end
    allaxes && return nothing
    return _walk_cellat(m, g, p)
end

# On the plane the chart cell is the cell. On the sphere a cell's edges are
# geodesics, not the chart's parallels: a point just across a parallel can
# lie under the geodesic bulge of the cell on the other side.
_repair(::Planar, g, f, i, j, p) = faceposition(g, f, i, j)
function _repair(m, g, f, i, j, p)
    k = faceposition(g, f, i, j)
    cell_contains(g, k, p) && return k
    for s in 1:4
        q = _across(g, f, i, j, s)
        q == 0 && continue
        cell_contains(g, q, p) && return q
    end
    for q in neighbors(g, k)
        cell_contains(g, q, p) && return q
    end
    return k
end

# squared chord / Euclidean distance from the mean of a cell's corners to p
@inline function _corner_dist(m, g, k, p)
    f, i, j = facecell(g, k)
    c = _lattice_corners(g, f, i, j)
    if iscurved(m)
        s = c[1] + c[2] + c[3] + c[4]
        return -(s[1] * p[1] + s[2] * p[2] + s[3] * p[3]) / _norm3(s)
    end
    x = (Float64(GI.x(c[1])) + Float64(GI.x(c[2])) + Float64(GI.x(c[3])) + Float64(GI.x(c[4]))) / 4
    y = (Float64(GI.y(c[1])) + Float64(GI.y(c[2])) + Float64(GI.y(c[3])) + Float64(GI.y(c[4]))) / 4
    return (x - Float64(GI.x(p)))^2 + (y - Float64(GI.y(p)))^2
end

function _walk_cellat(m, g, p)
    ni, nj = facesize(g, 1)
    k = faceposition(g, 1, (ni + 1) ÷ 2, (nj + 1) ÷ 2)
    d = _corner_dist(m, g, k, p)
    for _ in 1:ncells(g)
        cell_contains(g, k, p) && return k
        f, i, j = facecell(g, k)
        best = 0; bd = d
        for s in 1:4
            q = _across(g, f, i, j, s)
            q == 0 && continue
            dq = _corner_dist(m, g, q, p)
            dq < bd && (bd = dq; best = q)
        end
        best == 0 && break
        k = best; d = bd
    end
    return GrI.fallback_cellat(m, g, p)
end
