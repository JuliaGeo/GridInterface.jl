# ===========================================================================
# Vertex identity and the dual grid.
#
# A vertex is a point where three or more cells meet, or a corner on the outer
# boundary of a regional grid. A corner where only two cells meet (H3's class
# III bend points; a densification point) is not a vertex, and a repeated
# corner (a lon/lat pole cell) is one vertex. So `cell_vertices` is a
# subsequence of `cell_corners`.
#
# Identity has two layers:
#  - keys: `VertexKey(owner, k)`, corner `k` (of `cell_corners`) of the
#    canonical owner cell (the lowest position meeting at the vertex). A
#    provider with analytic vertices (`vertices_are_analytic`) names them
#    locally in O(1), with no table.
#  - the dense numbering `1:nvertices(grid)`: the provider's own
#    (`vertex_numbering_is_analytic`), derived from the keys (an O(n) count of
#    owned corners), or read from the fallback vertex table.
#
# Neither precomputation is cached here: `vertex_table` and `vertex_numbering`
# build per call, and the caller holds the result through `GridAndVertices`
# (which the dense verbs read) or `dual` (which stores it in the DualGrid).
# A provider with its own numbering needs neither, and nothing is built for it.
# ===========================================================================

"""
    abstract type Conformity

`Conforming()`: cells meet edge to edge (no hanging vertices) and tessellate,
so vertex incidence is well formed. `NonConforming()`: the dual family is
refused. Answered by [`conformity`](@ref).
"""
abstract type Conformity end
"Cells meet edge to edge and tessellate: vertices and the dual exist."
struct Conforming <: Conformity end
"Hanging vertices or overlaps: the dual family is refused."
struct NonConforming <: Conformity end

"""
    conformity(grid) -> Conformity

Default: conforming iff `SingleRing` and `Tessellating`. A multi-level set is
non-conforming (hanging vertices).
"""
conformity(grid) = _conformity(cellshape(grid), tessellation(grid))
_conformity(::SingleRing, ::Tessellating) = Conforming()
_conformity(::CellShape, ::Tessellation) = NonConforming()

@noinline _refuse_dual(grid) = throw(ArgumentError(
    "vertex identity and the dual are refused for $(typeof(grid)): it declares " *
    "$(cellshape(grid)), $(tessellation(grid)), $(conformity(grid))"))

"""
    VertexKey(cell, corner)

Typed vertex key: corner `corner` (an index into [`cell_corners`](@ref)) of
cell position `cell`. Canonical when `cell` is the lowest position meeting at
the vertex. Ordered by `(cell, corner)`, which is the order of the derived
dense numbering.
"""
struct VertexKey
    cell::Int
    corner::Int
end
Base.isless(a::VertexKey, b::VertexKey) = isless((a.cell, a.corner), (b.cell, b.corner))

"""
    vertices_are_analytic(grid) -> Bool

Whether the grid names its vertices itself through keys: it defines
[`vertex_key`](@ref)`(grid, i, k)` and [`vertex_cells`](@ref)`(grid, key)` in
O(1) with no table, and the dense numbering is derived. Default `false`: the
geometric fallback merges corners into a [`VertexTable`](@ref).
"""
vertices_are_analytic(grid) = false

"""
    vertex_numbering_is_analytic(grid) -> Bool

Whether the grid answers the dense vertex verbs itself in O(1), by index
arithmetic or its own connectivity, so no numbering or table is ever built
for it. With analytic keys it defines `nvertices`, `vertex_key(grid, v)`,
and `vertex_index(grid, key)`; without keys it defines `nvertices`,
`cell_vertices`, `vertex_cells(grid, v)`, `vertex_point(grid, v)`, and
`vertex_isopen(grid, v)`, and the keyed verbs are derived from those.
Default `false`: the numbering is derived from the keys, or the geometric
table is built.
"""
vertex_numbering_is_analytic(grid) = false

"""
    VertexTable

Merged vertex incidence of a conforming single-ring grid, built by the
geometric fallback (extension) and held by [`GridAndVertices`](@ref) or the
[`DualGrid`](@ref).
"""
struct VertexTable{P}
    points::Vector{P}           # vertex v -> its point
    cell_ptr::Vector{Int}       # cell i -> cell_verts[cell_ptr[i]:cell_ptr[i+1]-1]
    cell_verts::Vector{Int}
    cell_corner::Vector{Int}    # parallel to cell_verts: the corner index in cell_corners(grid, i)
    vert_ptr::Vector{Int}       # vertex v -> vert_cells[vert_ptr[v]:vert_ptr[v+1]-1], CCW (fan order when open)
    vert_cells::Vector{Int}
    isopen::BitVector
end

"""
    vertex_table(grid) -> VertexTable

The geometric vertex incidence (fallback), built per call in O(ncells); a
[`GridAndVertices`](@ref) answers its stored table. Grids with analytic vertex
identity never reach this.
"""
function vertex_table(grid)
    conformity(grid) isa Conforming || _refuse_dual(grid)
    P = _pointtype(grid)
    return fallback_vertex_table(best_manifold(grid), grid)::VertexTable{P}
end
_pointtype(grid) = eltype(Base.promote_op(cell_corners, typeof(grid), Int))

# --- the derived dense numbering for analytic vertices -------------------------

struct VertexNumbering
    offsets::Vector{Int}        # vertices owned by cells 1:i-1
end
"The dense numbering of a grid with analytic vertices, derived per call in O(ncells)."
function vertex_numbering(grid)
    n = ncells(grid)
    offs = Vector{Int}(undef, n + 1)
    offs[1] = 0
    for i in 1:n
        offs[i+1] = offs[i] + _nowned(grid, i)
    end
    return VertexNumbering(offs)
end
"The precomputation the dense vertex verbs of `grid` need: its numbering, its table, or nothing."
_vertex_precomputation(grid) = vertex_numbering_is_analytic(grid) ? nothing :
    vertices_are_analytic(grid) ? vertex_numbering(grid) : vertex_table(grid)
function _nowned(grid, i)
    c = 0
    for k in 1:length(cell_corners(grid, i))
        key = vertex_key(grid, i, k)
        key !== nothing && key.cell == i && key.corner == k && (c += 1)
    end
    return c
end

# --- keys derived from a provider's own dense connectivity ------------------------
# (`vertex_numbering_is_analytic` without `vertices_are_analytic`)

_derives_keys(grid) = vertex_numbering_is_analytic(grid) && !vertices_are_analytic(grid)
@noinline _missing_dense(grid, f) = throw(ArgumentError(
    "$(typeof(grid)) declares vertex_numbering_is_analytic but defines no `$f(grid, ...)` dense method"))

# the canonical key of dense vertex v: its corner on the lowest cell of its fan
function _dense_key(grid, v)
    c = minimum(vertex_cells(grid, v))
    vs = cell_vertices(grid, c)
    j = findfirst(==(v), vs)
    j === nothing && error("vertex $v not on its owner cell $c")
    cs = cell_corners(grid, c)
    length(cs) == length(vs) && return VertexKey(c, j)
    p = vertex_point(grid, v)
    k = findfirst(q -> _near(q, p), cs)
    k === nothing && error("vertex $v is no corner of its owner cell $c")
    return VertexKey(c, k)
end
# the dense vertex at corner key.corner of cell key.cell, or 0 when that corner is not a vertex
function _dense_index(grid, key::VertexKey)
    vs = cell_vertices(grid, key.cell)
    cs = cell_corners(grid, key.cell)
    length(cs) == length(vs) && return Int(vs[key.corner])
    p = cs[key.corner]
    for v in vs
        _near(vertex_point(grid, v), p) && return Int(v)
    end
    return 0
end

# --- the dense API ----------------------------------------------------------------

"""
    nvertices(grid) -> Int

Dense count of the grid's vertices (conforming single-ring grids). Analytic
vertices: derived from the keys; otherwise read from the geometric vertex
table. Either is O(ncells) per call on a plain grid: hold it through
[`GridAndVertices`](@ref) where the dense verbs are called repeatedly.
"""
function nvertices(grid)
    vertex_numbering_is_analytic(grid) && _missing_dense(grid, :nvertices)
    return vertices_are_analytic(grid) ? vertex_numbering(grid).offsets[end] : length(vertex_table(grid).points)
end

"""
    cell_vertices(grid, i) -> AbstractVector{Int}

Dense vertex numbers of cell `i`, counter-clockwise: a subsequence of
[`cell_corners`](@ref) (corners where only two cells meet are dropped,
repeated corners appear once).
"""
function cell_vertices(grid, i::Integer)
    vertices_are_analytic(grid) && return MappedVec(key -> vertex_index(grid, key), cell_vertex_keys(grid, i))
    vertex_numbering_is_analytic(grid) && _missing_dense(grid, :cell_vertices)
    t = vertex_table(grid)
    return view(t.cell_verts, t.cell_ptr[i]:t.cell_ptr[i+1]-1)
end

"""
    vertex_cells(grid, v) -> AbstractVector{Int}
    vertex_cells(grid, key::VertexKey)

Positions of the cells meeting at vertex `v` (dense number or key),
counter-clockwise around it; for an open fan ([`vertex_isopen`](@ref)), in fan
order from the boundary edge where the fan starts.
"""
function vertex_cells(grid, v::Integer)
    vertices_are_analytic(grid) && return vertex_cells(grid, vertex_key(grid, v))
    vertex_numbering_is_analytic(grid) && _missing_dense(grid, :vertex_cells)
    t = vertex_table(grid)
    return view(t.vert_cells, t.vert_ptr[v]:t.vert_ptr[v+1]-1)
end
vertex_cells(grid, key::VertexKey) = vertices_are_analytic(grid) ? _missing_keyed(grid, :vertex_cells) :
    vertex_cells(grid, vertex_index(grid, key))
@noinline _missing_keyed(grid, f) = throw(ArgumentError(
    "$(typeof(grid)) declares vertices_are_analytic but defines no `$f(grid, ...)` keyed method"))

"""
    vertex_point(grid, v)
    vertex_point(grid, key::VertexKey)

The point of a vertex.
"""
function vertex_point(grid, v::Integer)
    vertices_are_analytic(grid) && return vertex_point(grid, vertex_key(grid, v))
    vertex_numbering_is_analytic(grid) && _missing_dense(grid, :vertex_point)
    return vertex_table(grid).points[v]
end
vertex_point(grid, key::VertexKey) = cell_corners(grid, key.cell)[key.corner]

"""
    vertex_isopen(grid, v) -> Bool
    vertex_isopen(grid, key::VertexKey)

Whether the fan of cells around a vertex is open (the vertex lies on the outer
boundary of a regional grid). The dual cell of an open vertex is the median
dual: its ring holds points that are no primal cell's centroid, so
[`dual_ring_cells`](@ref) (not `vertex_cells`) pairs ring points with cells.
The boundary band (between the outer centroids and the grid edge) is consumer
policy.
"""
function vertex_isopen(grid, v::Integer)
    vertices_are_analytic(grid) && return vertex_isopen(grid, vertex_key(grid, v))
    vertex_numbering_is_analytic(grid) && _missing_dense(grid, :vertex_isopen)
    return vertex_table(grid).isopen[v]
end
function vertex_isopen(grid, key::VertexKey)
    _derives_keys(grid) && return vertex_isopen(grid, vertex_index(grid, key))
    isglobal(grid) && return false
    cs = vertex_cells(grid, key)
    m = length(cs)
    m < 3 && return true
    p = vertex_point(grid, key)
    for k in 1:m
        _prev_next_corner(grid, cs[k], p)[1] == _prev_next_corner(grid, cs[mod1(k + 1, m)], p)[2] || return true
    end
    return false
end
# keys of the vertices before and after point p around cell c
function _prev_next_corner(grid, c, p)
    ks = cell_vertex_keys(grid, c)
    n = length(ks)
    for j in 1:n
        _near(vertex_point(grid, ks[j]), p) && return (ks[mod1(j - 1, n)], ks[mod1(j + 1, n)])
    end
    throw(ArgumentError("vertex not on cell $c"))
end

# --- the keyed API ------------------------------------------------------------------

"""
    vertex_key(grid, i, k) -> Union{VertexKey,Nothing}
    vertex_key(grid, v::Integer) -> VertexKey

The canonical key of corner `k` of cell `i` (`nothing` when that corner is not
a vertex), or of dense vertex `v`. A grid with [`vertices_are_analytic`](@ref)
defines the three-argument form (and the two-argument one when its numbering
is analytic too); a grid with its own dense connectivity has both derived from
it; otherwise they are read from the fallback table.
"""
function vertex_key(grid, i::Integer, k::Integer)
    vertices_are_analytic(grid) && _missing_keyed(grid, :vertex_key)
    if vertex_numbering_is_analytic(grid)
        v = _dense_index(grid, VertexKey(i, k))
        return v == 0 ? nothing : _dense_key(grid, v)
    end
    t = vertex_table(grid)
    for j in t.cell_ptr[i]:t.cell_ptr[i+1]-1
        t.cell_corner[j] == k && return _table_key(t, t.cell_verts[j])
    end
    return nothing
end
function vertex_key(grid, v::Integer)
    if vertices_are_analytic(grid)
        vertex_numbering_is_analytic(grid) && _missing_dense(grid, :vertex_key)
        offs = vertex_numbering(grid).offsets
        1 <= v <= offs[end] || throw(BoundsError(grid, v))
        c = searchsortedlast(offs, v - 1)
        r = v - offs[c]
        for k in 1:length(cell_corners(grid, c))
            key = vertex_key(grid, c, k)
            key !== nothing && key.cell == c && key.corner == k && (r -= 1) == 0 && return key
        end
        error("inconsistent vertex numbering")
    end
    vertex_numbering_is_analytic(grid) && return _dense_key(grid, v)
    return _table_key(vertex_table(grid), v)
end
function _table_key(t::VertexTable, v)
    c = t.vert_cells[t.vert_ptr[v]]
    for j in t.vert_ptr[v]:t.vert_ptr[v+1]-1
        c = min(c, t.vert_cells[j])
    end
    for j in t.cell_ptr[c]:t.cell_ptr[c+1]-1
        t.cell_verts[j] == v && return VertexKey(c, t.cell_corner[j])
    end
    error("vertex $v not on its owner")
end

"""
    vertex_index(grid, key::VertexKey) -> Int

The dense number of a vertex key (any corner key of the vertex is accepted).
"""
function vertex_index(grid, key::VertexKey)
    if vertices_are_analytic(grid)
        vertex_numbering_is_analytic(grid) && _missing_dense(grid, :vertex_index)
        key = vertex_key(grid, key.cell, key.corner)
        key === nothing && throw(ArgumentError("corner is not a vertex"))
        c = key.cell
        r = vertex_numbering(grid).offsets[c]
        for k in 1:key.corner
            kk = vertex_key(grid, c, k)
            kk !== nothing && kk.cell == c && kk.corner == k && (r += 1)
        end
        return r
    end
    if vertex_numbering_is_analytic(grid)
        v = _dense_index(grid, key)
        v == 0 || return v
    else
        t = vertex_table(grid)
        for j in t.cell_ptr[key.cell]:t.cell_ptr[key.cell+1]-1
            t.cell_corner[j] == key.corner && return t.cell_verts[j]
        end
    end
    throw(ArgumentError("corner $(key.corner) of cell $(key.cell) is not a vertex"))
end

"""
    cell_vertex_keys(grid, i) -> AbstractVector{VertexKey}

The canonical keys of cell `i`'s vertices, in [`cell_vertices`](@ref) order.
Generic for analytic vertices: one key per corner (lazy).
"""
function cell_vertex_keys(grid, i::Integer)
    if vertices_are_analytic(grid)
        return FuncVec(k -> vertex_key(grid, i, k)::VertexKey, length(cell_corners(grid, i)))
    end
    vertex_numbering_is_analytic(grid) && return MappedVec(v -> _dense_key(grid, v), cell_vertices(grid, i))
    t = vertex_table(grid)
    return MappedVec(v -> _table_key(t, v), view(t.cell_verts, t.cell_ptr[i]:t.cell_ptr[i+1]-1))
end

# ---------------------------------------------------------------------------
# The dual grid.
# ---------------------------------------------------------------------------

"""
    DualGrid(primal)

One cell per primal vertex, numbered by the dense vertex numbering. A closed
fan's ring is the centroids of `vertex_cells(primal, v)` in order (lazy) and its
centroid is the primal vertex. An open fan (a boundary vertex of a regional
grid) is closed the median-dual way: the vertex, the midpoint of the boundary
edge where the fan starts, the centroids in order, the midpoint of the boundary
edge where it ends ([`dual_ring_cells`](@ref) gives the correspondence).
Declared `SimpleCells`. The dual of a global conforming grid tiles the sphere.
"""
struct DualGrid{G,T} <: AbstractGrid
    primal::G
    # the primal's VertexTable or VertexNumbering (the dual owns its
    # precomputation), or nothing when the primal numbers its vertices itself
    table::T
end

"""
    dual(grid) -> DualGrid

The dual grid, holding the primal's vertex precomputation (built once here;
nothing when the primal's numbering is analytic). Refused up front for `GeneralCells`, overlapping, or non-conforming grids.
"""
function dual(grid)
    conformity(grid) isa Conforming || _refuse_dual(grid)
    return DualGrid(grid, _vertex_precomputation(grid))
end

# The primal's dense vertex verbs, answered from the stored precomputation, or
# by the primal itself when it numbers its own vertices.
_vg(d::DualGrid) = GridAndVertices(d.primal, d.table)
_vg(d::DualGrid{G,Nothing}) where {G} = d.primal
_nv(d::DualGrid) = nvertices(_vg(d))
_vpoint(d::DualGrid, v) = vertex_point(_vg(d), v)
_vcells(d::DualGrid, v) = vertex_cells(_vg(d), v)
_cverts(d::DualGrid, c) = cell_vertices(_vg(d), c)
_isopen(d::DualGrid, v) = vertex_isopen(_vg(d), v)

Base.parent(d::DualGrid) = d.primal
ncells(d::DualGrid) = _nv(d)
GOC.best_manifold(d::DualGrid) = best_manifold(d.primal)
cellshape(::DualGrid) = SimpleCells()
isglobal(d::DualGrid) = isglobal(d.primal)
cell_centroid(d::DualGrid, v::Integer) = _vpoint(d, v)
location_is_analytic(d::DualGrid) = location_is_analytic(d.primal) && vertices_are_analytic(d.primal)

# prev/next vertex of v around cell c (cell rings are CCW); repeats stepped over
function _ring_neighbours(d, c, v)
    vs = _cverts(d, c)
    k = findfirst(==(v), vs)
    n = length(vs)
    p = k; while p > k - n && vs[mod1(p - 1, n)] == v; p -= 1; end
    q = k; while q < k + n && vs[mod1(q + 1, n)] == v; q += 1; end
    return vs[mod1(p - 1, n)], vs[mod1(q + 1, n)]
end

# An open fan in walking order: (cells, vertex across the start edge, vertex across the end edge)
function _open_fan(d::DualGrid, v)
    cs = collect(_vcells(d, v))
    pn = [_ring_neighbours(d, c, v) for c in cs]
    start = findfirst(k -> !any(q -> q[1] == pn[k][2], pn), eachindex(cs))
    start === nothing && (start = 1)
    order = Int[start]
    while length(order) < length(cs)
        prev_edge = pn[order[end]][1]
        k = findfirst(q -> q[2] == prev_edge, pn)
        (k === nothing || k in order) && break
        push!(order, k)
    end
    return cs[order], pn[start][2], pn[order[end]][1]
end

function cell_boundary(d::DualGrid, v::Integer)
    1 <= v <= ncells(d) || _throw_bounds(d, v)
    g = d.primal
    P = typeof(_vpoint(d, v))
    if !_isopen(d, v)        # closed fan: the centroids in order, lazily
        return _centroid_ring(g, P, _vcells(d, v))
    end
    cells, a, b = _open_fan(d, v)
    vp = _vpoint(d, v)
    m = best_manifold(g)
    pts = P[vp, _midpoint(m, vp, _vpoint(d, a))]
    for c in cells
        push!(pts, convert(P, cell_centroid(g, c)))
    end
    push!(pts, _midpoint(m, vp, _vpoint(d, b)))
    return pts
end
_centroid_ring(g, ::Type{P}, cs) where {P} = MappedVec(c -> convert(P, cell_centroid(g, c)), cs)

"""
    dual_ring_cells(d::DualGrid, v) -> AbstractVector{Int}

For each point of `cell_boundary(d, v)`, the primal cell whose centroid it is,
or [`OUTSIDE`](@ref) (`0`) for the vertex itself and the edge midpoints of an
open fan. Closed fans: exactly `vertex_cells(primal, v)`.
"""
function dual_ring_cells(d::DualGrid, v::Integer)
    _isopen(d, v) || return collect(Int, _vcells(d, v))
    cells, _, _ = _open_fan(d, v)
    return vcat(OUTSIDE, OUTSIDE, cells, OUTSIDE)
end

"""
    cellat(d::DualGrid, p)

Dual point location through the primal host: locate `p` in the primal, then
test the dual cells of the host's vertices (and, rarely, of its neighbours'
vertices); the lowest containing vertex wins. Analytic when the primal's
location and vertices are.
"""
function cellat(d::DualGrid, p)
    m = best_manifold(d)
    check_point(m, p)
    h = cellat(d.primal, p)
    h === nothing && return nothing
    best = _lowest_containing(d, _cverts(d, h), p, 0)
    best == 0 || return best
    for c in neighbors(d.primal, h)
        best = _lowest_containing(d, _cverts(d, c), p, best)
    end
    return best == 0 ? nothing : best
end
function _lowest_containing(d, vs, p, best)
    for v in vs
        (best == 0 || v < best) && cell_contains(d, v, p) && (best = v)
    end
    return best
end

"""
    vertexat(grid, p) -> Union{VertexKey,Nothing}

The key of the dual cell containing `p`, by the primal-host method, without
building the dense numbering or the dual grid: what a point-based regridder on
an analytic-vertex grid calls per destination point. Closed fans only (a point
in an open fan's band returns `nothing`).
"""
function vertexat(grid, p)
    m = best_manifold(grid)
    check_point(m, p)
    h = cellat(grid, p)
    h === nothing && return nothing
    P = _pointtype(grid)
    best = nothing
    for key in cell_vertex_keys(grid, h)
        (best === nothing || key < best) || continue
        cs = vertex_cells(grid, key)
        length(cs) >= 3 || continue
        ring = MappedVec(c -> convert(P, cell_centroid(grid, c)), cs)
        fallback_ring_contains(m, ring, p) && (best = key)
    end
    return best
end
fallback_ring_contains(::Any, ring, p) = throw(NoGeometryBackend(:vertexat))

_midpoint(::Planar, a, b) = _like(a, (_coords(a)[1] + _coords(b)[1]) / 2, (_coords(a)[2] + _coords(b)[2]) / 2)
function _midpoint(::Spherical, a, b)
    x = _coords(a); y = _coords(b)
    s = (x[1] + y[1], x[2] + y[2], x[3] + y[3]); r = sqrt(s[1]^2 + s[2]^2 + s[3]^2)
    return _like(a, s[1] / r, s[2] / r, s[3] / r)
end
_like(a::Tuple, x...) = x
_like(a, x...) = typeof(a)(x...)
