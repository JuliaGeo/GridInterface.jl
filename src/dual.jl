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
    lock(_VERTEX_LOCK) do
        get!(() -> fallback_vertex_table(best_manifold(grid), grid), _VERTEX_CACHE, grid)
    end
end

"""
    nvertices(grid) -> Int

Dense count of the grid's vertices. Optional family on conforming single-ring
grids; fallback (extension) merges coincident boundary points once and caches.
"""
nvertices(grid) = length(vertex_table(grid).points)

"""
    cell_vertices(grid, i) -> AbstractVector{Int}

Vertex numbers of cell `i`, in the same counter-clockwise order as
[`cell_boundary`](@ref).
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
inherits the primal manifold. At a vertex with fewer than three cells (the
outer boundary of an open grid) the vertex itself is inserted so the partial
polygon has a chance of being non-degenerate; it is declared `SimpleCells`.
"""
struct DualGrid{G} <: AbstractGrid
    primal::G
end

"""
    dual(grid) -> DualGrid

The dual grid. Refused up front for `GeneralCells`, overlapping, or
non-conforming grids.
"""
function dual(grid)
    conformity(grid) isa Conforming || _refuse_dual(grid)
    return DualGrid(grid)
end

Base.parent(d::DualGrid) = d.primal
ncells(d::DualGrid) = nvertices(d.primal)
GOC.best_manifold(d::DualGrid) = best_manifold(d.primal)
cellshape(::DualGrid) = SimpleCells()
function cell_boundary(d::DualGrid, v::Integer)
    1 <= v <= ncells(d) || _throw_bounds(d, v)
    cs = vertex_cells(d.primal, v)
    pts = [cell_centroid(d.primal, c) for c in cs]
    length(pts) < 3 && pushfirst!(pts, vertex_point(d.primal, v))
    return pts
end
cell_centroid(d::DualGrid, v::Integer) = vertex_point(d.primal, v)
