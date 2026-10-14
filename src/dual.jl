# ===========================================================================
# Vertex identity and the dual grid.
# ===========================================================================

"""
    abstract type Conformity

`Conforming()`: cells meet edge to edge (no hanging vertices) and tessellate,
so vertex incidence is well formed. `NonConforming()`: the dual family is
refused. Default: conforming iff `SingleRing` and `Tessellating`. A multi-level
set is non-conforming (hanging vertices). Name chosen here (PRD: "conformity
trait", name to be settled).
"""
abstract type Conformity end
struct Conforming <: Conformity end
struct NonConforming <: Conformity end

conformity(grid) = _conformity(cellshape(grid), tessellation(grid))
_conformity(::SingleRing, ::Tessellating) = Conforming()
_conformity(::CellShape, ::Tessellation) = NonConforming()

@noinline _refuse_dual(grid) = throw(ArgumentError(
    "vertex identity and the dual are refused for $(typeof(grid)): it declares " *
    "$(cellshape(grid)), $(tessellation(grid)), $(conformity(grid))"))

"""
    VertexTable

Merged vertex incidence of a single-ring grid, built once by the geometric
fallback (extension) and cached per grid object.
"""
struct VertexTable{P}
    points::Vector{P}           # vertex v -> its point
    cell_ptr::Vector{Int}       # cell i -> cell_verts[cell_ptr[i]:cell_ptr[i+1]-1]
    cell_verts::Vector{Int}
    vert_ptr::Vector{Int}       # vertex v -> vert_cells[vert_ptr[v]:vert_ptr[v+1]-1], CCW
    vert_cells::Vector{Int}
end

const _VERTEX_CACHE = IdDict{Any,Any}()
const _VERTEX_LOCK = ReentrantLock()

"""
    vertex_table(grid) -> VertexTable

The cached geometric vertex incidence (fallback). Grids with analytic vertex
identity never reach this.
"""
function vertex_table(grid)
    conformity(grid) isa Conforming || _refuse_dual(grid)
    t = lock(_VERTEX_LOCK) do
        get!(() -> fallback_vertex_table(best_manifold(grid), grid), _VERTEX_CACHE, grid)
    end
    return t::VertexTable{_pointtype(grid)}
end
_pointtype(grid) = eltype(Base.promote_op(cell_corners, typeof(grid), Int))

"""
    nvertices(grid) -> Int

Dense count of the grid's vertices. Optional family on conforming single-ring
grids; fallback (extension) merges coincident boundary points once and caches.
"""
nvertices(grid) = length(vertex_table(grid).points)

"""
    cell_vertices(grid, i) -> AbstractVector{Int}

Vertex numbers of cell `i`, in the same counter-clockwise order as
[`cell_corners`](@ref) (which is `cell_boundary` unless the grid densifies).
"""
function cell_vertices(grid, i::Integer)
    t = vertex_table(grid)
    return view(t.cell_verts, t.cell_ptr[i]:t.cell_ptr[i+1]-1)
end

"""
    vertex_cells(grid, v) -> AbstractVector{Int}

Positions of the cells meeting at vertex `v`, counter-clockwise around it.
"""
function vertex_cells(grid, v::Integer)
    t = vertex_table(grid)
    return view(t.vert_cells, t.vert_ptr[v]:t.vert_ptr[v+1]-1)
end

"""
    vertex_point(grid, v)

The point of vertex `v`. (Not named in the PRD; the dual's centroid needs it.)
"""
vertex_point(grid, v::Integer) = vertex_table(grid).points[v]

"""
    DualGrid(primal)

One cell per primal vertex. Its boundary is the centroids of
`vertex_cells(primal, v)` in order, its centroid is the primal vertex, and it
inherits the primal manifold.

At a vertex on the outer boundary of an open grid the fan of cells does not
close, and the "partial polygon those centroids span" (PRD) is degenerate at
corners (one centroid) and wound arbitrarily along edges (two centroids). This
implementation closes the fan the median-dual way instead: the vertex, the
midpoint of the boundary edge where the fan starts, the centroids in order, the
midpoint of the boundary edge where it ends. Declared `SimpleCells`.
"""
struct DualGrid{G,T} <: AbstractGrid
    primal::G
    table::T       # the primal's fallback VertexTable, or `nothing` when the primal answers analytically
end

"""
    dual(grid) -> DualGrid

The dual grid. Refused up front for `GeneralCells`, overlapping, or
non-conforming grids.
"""
function dual(grid)
    conformity(grid) isa Conforming || _refuse_dual(grid)
    return DualGrid(grid, uses_vertex_fallback(grid) ? vertex_table(grid) : nothing)
end

"Whether `nvertices` on this grid is answered by the geometric fallback."
uses_vertex_fallback(grid) = which(nvertices, Tuple{typeof(grid)}).sig.parameters[2] === Any

_nv(d::DualGrid{G,Nothing}) where {G} = nvertices(d.primal)
_nv(d::DualGrid) = length(d.table.points)
_vpoint(d::DualGrid{G,Nothing}, v) where {G} = vertex_point(d.primal, v)
_vpoint(d::DualGrid, v) = @inbounds d.table.points[v]
_vcells(d::DualGrid{G,Nothing}, v) where {G} = vertex_cells(d.primal, v)
_vcells(d::DualGrid, v) = (t = d.table; view(t.vert_cells, t.vert_ptr[v]:t.vert_ptr[v+1]-1))
_cverts(d::DualGrid{G,Nothing}, c) where {G} = cell_vertices(d.primal, c)
_cverts(d::DualGrid, c) = (t = d.table; view(t.cell_verts, t.cell_ptr[c]:t.cell_ptr[c+1]-1))

Base.parent(d::DualGrid) = d.primal
ncells(d::DualGrid) = _nv(d)
GOC.best_manifold(d::DualGrid) = best_manifold(d.primal)
cellshape(::DualGrid) = SimpleCells()
cell_centroid(d::DualGrid, v::Integer) = _vpoint(d, v)

# prev/next vertex of v around cell c (cell rings are CCW)
function _ring_neighbours(d, c, v)
    vs = _cverts(d, c)
    k = findfirst(==(v), vs)
    n = length(vs)
    # a degenerate ring may repeat v (lon/lat polar cells): step past the repeats
    p = k; while p > k - n && vs[mod1(p - 1, n)] == v; p -= 1; end
    q = k; while q < k + n && vs[mod1(q + 1, n)] == v; q += 1; end
    return vs[mod1(p - 1, n)], vs[mod1(q + 1, n)]
end

"""
    vertex_isopen(grid, v) -> Bool

Whether the fan of cells around vertex `v` is open (`v` lies on the outer
boundary of a regional grid). The dual cell of an open vertex is the median
dual: its ring carries points (the vertex, edge midpoints) that are not source
centroids, so `vertex_cells(grid, v)` and `cell_boundary(dual(grid), v)` no
longer correspond one to one. Fallback: from the vertex table.
"""
function vertex_isopen(grid, v::Integer)
    cs = vertex_cells(grid, v)
    m = length(cs)
    m < 3 && return true
    # closed iff the shared-edge law closes the cycle: prev(c_k) == next(c_{k+1})
    for k in 1:m
        _prevnext(grid, cs[k], v)[1] == _prevnext(grid, cs[k == m ? 1 : k + 1], v)[2] || return true
    end
    return false
end
function _prevnext(grid, c, v)
    vs = cell_vertices(grid, c)
    k = findfirst(==(v), vs)
    n = length(vs)
    return vs[mod1(k - 1, n)], vs[mod1(k + 1, n)]
end

function cell_boundary(d::DualGrid, v::Integer)
    1 <= v <= ncells(d) || _throw_bounds(d, v)
    g = d.primal
    P = typeof(_vpoint(d, v))
    if !vertex_isopen(g, v)      # closed fan: the centroids in order, nothing else
        return P[convert(P, cell_centroid(g, c)) for c in _vcells(d, v)]
    end
    cs = collect(_vcells(d, v))
    pn = [_ring_neighbours(d, c, v) for c in cs]
    # the fan is open iff some cell's "next" edge is nobody's "prev" edge
    start = findfirst(k -> !any(q -> q[1] == pn[k][2], pn), eachindex(cs))
    if start === nothing
        return P[convert(P, cell_centroid(g, c)) for c in cs]
    end
    # walk the open fan counter-clockwise from its starting boundary edge
    order = Int[start]
    while length(order) < length(cs)
        prev_edge = pn[order[end]][1]
        k = findfirst(q -> q[2] == prev_edge, pn)
        (k === nothing || k in order) && break
        push!(order, k)
    end
    vp = _vpoint(d, v)
    pts = P[vp, _midpoint(best_manifold(g), vp, _vpoint(d, pn[start][2]))]
    for k in order
        push!(pts, convert(P, cell_centroid(g, cs[k])))
    end
    push!(pts, _midpoint(best_manifold(g), vp, _vpoint(d, pn[order[end]][1])))
    return pts
end

_midpoint(::Planar, a, b) = typeof(a)(((a[1] + b[1]) / 2, (a[2] + b[2]) / 2))
_midpoint(::Spherical, a, b) = (s = (a[1] + b[1], a[2] + b[2], a[3] + b[3]); r = sqrt(s[1]^2 + s[2]^2 + s[3]^2);
    _like(a, s[1] / r, s[2] / r, s[3] / r))
_like(a::Tuple, x...) = x
_like(a, x...) = typeof(a)(x...)
