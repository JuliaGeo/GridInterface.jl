# ===========================================================================
# Neighbours: every polygonal grid has them; a trait says whether they are
# analytic or derived. DGG's contract, lifted, plus the disc form.
# ===========================================================================

"""
    abstract type Connectivity

Which cells count as adjacent: [`Vertex`](@ref) or [`Edge`](@ref).
"""
abstract type Connectivity end
"Cells sharing at least a vertex (Moore). The default."
struct Vertex <: Connectivity end
"Cells sharing a whole edge (von Neumann)."
struct Edge <: Connectivity end

"""
    abstract type Winding

The order [`neighbors`](@ref) returns, answered by [`winding`](@ref).
"""
abstract type Winding end
"Neighbours come counter-clockwise around the cell, seen from outside the sphere / from above the plane."
struct CounterClockwise <: Winding end
"Neighbours come clockwise."
struct Clockwise <: Winding end
"Neighbours come in a deterministic, documented order that is not a rotation."
struct CustomOrder <: Winding end
"No order is promised."
struct Unordered <: Winding end

"""
    neighbors_are_analytic(grid) -> Bool

Cheapness trait: `true` when [`neighbors`](@ref) is answered from the grid's own
structure. A consumer sweeping many times over a grid where this is `false`
should precompute an adjacency once.
"""
neighbors_are_analytic(grid) = false

"""
    maxneighbors(grid, connectivity = Vertex()) -> Union{Int,Nothing}
    maxneighbors(grid, k, connectivity) -> Union{Int,Nothing}

Capacity hint for static stencil buffers: an upper bound on
`length(neighbors(grid, i; connectivity))`, or (three-argument form) on the
disc `length(neighbors(grid, i, k; connectivity))`, which is the sum of
[`maxring`](@ref) over `1:k`. `nothing` means unbounded/undeclared.
"""
maxneighbors(grid, ::Connectivity) = nothing
maxneighbors(grid) = maxneighbors(grid, Vertex())
function maxneighbors(grid, k::Integer, conn::Connectivity=Vertex())
    k >= 0 || throw(ArgumentError("disc radius must be non-negative"))
    total = 0
    for j in 1:k
        m = maxring(grid, j, conn)
        m === nothing && return nothing
        total += m
    end
    return total
end

"""
    maxring(grid, k, connectivity = Vertex()) -> Union{Int,Nothing}

Upper bound on `length(ring(grid, i, k; connectivity))`. Generic: `1` for
`k == 0`, `maxneighbors(grid, connectivity)` for `k == 1`, `nothing` beyond
(a provider with a ring law, e.g. `6k` on a hexagonal grid, declares it).
"""
function maxring(grid, k::Integer, conn::Connectivity=Vertex())
    k == 0 && return 1
    k == 1 && return maxneighbors(grid, conn)
    return nothing
end

"""
    winding(grid, connectivity = Vertex()) -> Winding

The order [`neighbors`](@ref) returns. Default: `CounterClockwise()` for a
tessellating single-ring grid answered by the geometric fallback (which sorts
by azimuth, starting at the edge from corner 1 to corner 2), `Unordered()`
otherwise.

On conforming single-ring grids, `CounterClockwise` Edge winding promises
increasing boundary-edge indices, where edge `k` runs from
`cell_corners(grid, i)[k]` to the next corner (wrapping at the end). Edges
without a neighbour, including neighbours clipped out of a [`Subset`](@ref),
are omitted. When every edge has a neighbour, neighbour `k` lies across edge
`k`. The compact vector has no placeholders for missing edges.
"""
winding(grid, conn::Connectivity) = _default_winding(cellshape(grid), grid)
winding(grid) = winding(grid, Vertex())
# (on an overlapping grid "touching" includes overlaps, which have no ring order)
_default_winding(::SingleRing, grid) =
    uses_neighbor_fallback(grid) && tessellation(grid) isa Tessellating ? CounterClockwise() : Unordered()
_default_winding(::CellShape, grid) = Unordered()

"Whether `neighbors` on this grid type resolves to the generic (geometric fallback) method."
uses_neighbor_fallback(grid) = _uses_nbr_fb(typeof(grid))
function _uses_nbr_fb(G)
    m = which(neighbors, Tuple{G,Int})
    return Base.unwrap_unionall(m.sig).parameters[2] === Any
end

"""
    neighbors(grid, i; connectivity = Vertex()) -> AbstractVector{Int}
    neighbors(grid, i, k::Integer; connectivity = Vertex())        # the disc
    neighbors(grid, i, sel::Index; connectivity = Vertex())        # in an index space

Positions of the cells sharing an edge (`Edge()`) or at least a vertex
(`Vertex()`) with cell `i`, excluding `i`. The disc form returns every cell
within `k` steps, rings `1:k` concatenated, nearest first. Analytic grids
override them; any other polygonal grid gets the geometric fallback
(extension): on `SingleRing` grids, cells sharing boundary points (ordered
CCW); on `GeneralCells` grids, cells **touching** cell `i` (sharing boundary
or overlapping in area), `Unordered`, and `Edge`/`Vertex` coincide.

On a [`Subset`](@ref) the plain and disc forms are **clipped** (local positions
of members only; distance measured in the parent); the [`Index`](@ref) forms are
**unclipped** (every parent slot, `Index(Local())` with [`OUTSIDE`](@ref) for
non-members, `Index(Global())` in complete-grid positions).
"""
neighbors(grid, i::Integer; connectivity::Connectivity=Vertex()) =
    _neighbors(cellshape(grid), grid, Int(i), connectivity)
_neighbors(::Polygonal, grid, i, conn) = fallback_neighbors(best_manifold(grid), grid, i, conn)

neighbors(grid, i::Integer, k::Integer; connectivity::Connectivity=Vertex()) =
    _disc(grid, Int(i), Int(k), connectivity)

# Breadth-first over `neighbors`; discs are small, so membership is a linear scan.
function _disc(grid, i, k, conn)
    k >= 0 || throw(ArgumentError("disc radius must be non-negative"))
    out = Int[]
    k == 0 && return out
    lo = 0                     # out[lo+1:end] is the current frontier
    for j in neighbors(grid, i; connectivity=conn)
        j == i || j in out || push!(out, j)
    end
    for _ in 2:k
        hi = length(out)
        for c in lo+1:hi, j in neighbors(grid, out[c]; connectivity=conn)
            (j == i || j in out) || push!(out, j)
        end
        lo = hi
    end
    return out
end

# Index forms. Complete grids: every space is the grid's own (typed ids map
# lazily). Provider subsets (a DGG PartialGrid): through the parent, unclipped.
neighbors(grid, i::Integer, ::Index{Local}; connectivity::Connectivity=Vertex()) =
    _local_nbrs(gridrelation(grid), grid, Int(i), connectivity)
_local_nbrs(::GridRelation, grid, i, conn) = neighbors(grid, i; connectivity=conn)
function _local_nbrs(::ParentGrid, grid, i, conn)
    P = parent(grid)
    nb = neighbors(P, localindex(P, cellindex(grid, i))::Int; connectivity=conn)
    return MappedVec(j -> something(localindex(grid, cellindex(P, j)), OUTSIDE), nb)
end
neighbors(grid, i::Integer, ::Index{Global}; connectivity::Connectivity=Vertex()) =
    _global_nbrs(gridrelation(grid), grid, Int(i), connectivity)
_global_nbrs(::GridRelation, grid, i, conn) = neighbors(grid, i; connectivity=conn)
function _global_nbrs(::ParentGrid, grid, i, conn)
    P = parent(grid)
    nb = neighbors(P, localindex(P, cellindex(grid, i))::Int; connectivity=conn)
    gridrelation(P) isa NoParent && return nb
    return MappedVec(j -> toindex(P, j, Index(Global())), nb)
end
neighbors(grid, i::Integer, ::Index{Type{T}}; connectivity::Connectivity=Vertex()) where {T} =
    _typed_nbrs(gridrelation(grid), grid, Int(i), T, connectivity)
_typed_nbrs(::GridRelation, grid, i, T, conn) =
    MappedVec(j -> cellindex(grid, j, T), neighbors(grid, i; connectivity=conn))
function _typed_nbrs(::ParentGrid, grid, i, T, conn)
    P = parent(grid)
    nb = neighbors(P, localindex(P, cellindex(grid, i))::Int; connectivity=conn)
    return MappedVec(j -> cellindex(P, j, T), nb)
end

"""
    ring(grid, i, k; connectivity = Vertex()) -> AbstractVector{Int}

Positions at adjacency distance exactly `k` from `i` (`ring(grid, i, 0) == [i]`).
Generic method: the last ring of the disc, in discovery order (so ring 1 has the
order of `neighbors`). On a [`Subset`](@ref), distance is measured in the
parent and only member cells are returned, just as for the disc form of
[`neighbors`](@ref); paths may pass through non-members.
"""
function ring(grid, i::Integer, k::Integer; connectivity::Connectivity=Vertex())
    k >= 0 || throw(ArgumentError("ring distance must be non-negative"))
    k == 0 && return Int[Int(i)]
    k == 1 && return collect(Int, neighbors(grid, i; connectivity))
    inner = length(neighbors(grid, i, k - 1; connectivity))
    d = neighbors(grid, i, k; connectivity)
    return d[inner+1:end]
end

# --- subsets: clipped forms (local positions) and unclipped Index forms ----------

# parent's neighbours of member i, mapped to local positions (OUTSIDE for non-members)
_local_slots(s::Subset, v) = MappedVec(j -> something(_local_of_parent(s, j), OUTSIDE), v)

neighbors(s::Subset, i::Integer; connectivity::Connectivity=Vertex()) =
    Unsentinel(_local_slots(s, neighbors(s.parent, _ppos(s, i); connectivity)))
neighbors(s::Subset, i::Integer, k::Integer; connectivity::Connectivity=Vertex()) =
    Unsentinel(_local_slots(s, neighbors(s.parent, _ppos(s, i), k; connectivity)))
ring(s::Subset, i::Integer, k::Integer; connectivity::Connectivity=Vertex()) =
    k == 0 ? Int[Int(i)] : collect(Int, Unsentinel(_local_slots(s, ring(s.parent, _ppos(s, i), k; connectivity))))
neighbors(s::Subset, i::Integer, ::Index{Local}; connectivity::Connectivity=Vertex()) =
    _local_slots(s, neighbors(s.parent, _ppos(s, i); connectivity))
neighbors(s::Subset, i::Integer, sel::Index{Global}; connectivity::Connectivity=Vertex()) =
    neighbors(s.parent, _ppos(s, i), sel; connectivity)
neighbors(s::Subset, i::Integer, sel::Index{Type{T}}; connectivity::Connectivity=Vertex()) where {T} =
    neighbors(s.parent, _ppos(s, i), sel; connectivity)
maxneighbors(s::Subset, c::Connectivity) = maxneighbors(s.parent, c)
maxring(s::Subset, k::Integer, c::Connectivity=Vertex()) = maxring(s.parent, k, c)
winding(s::Subset, c::Connectivity) = winding(s.parent, c)

# ---------------------------------------------------------------------------
# Region verbs. The primitive is on the parent and a set of its positions, so
# the provider that owns the structure overrides it (DGG's O(border) subtree
# engines for rooted chunks); the Subset forms route through it.
# ---------------------------------------------------------------------------

@inline function _member(positions, p::Int)
    positions isa AbstractUnitRange && return p in positions
    k = searchsortedfirst(positions, p)
    return k <= length(positions) && @inbounds(positions[k]) == p
end

"""
    halo(parent, positions, k = 1; connectivity = Vertex()) -> Vector{Int}
    halo(s::Subset, k = 1; connectivity = Vertex())

Ascending **parent** positions of the cells outside the set `positions`
(ascending parent positions) within `k` steps of a member. The
three-argument form on the parent is the primitive a provider overrides; the
subset form calls it with `parentpositions(s)`. A complete grid has no outside.
"""
function halo(P, positions::AbstractVector{<:Integer}, k::Integer=1; connectivity::Connectivity=Vertex())
    out = Int[]
    frontier = positions
    for _ in 1:k
        new = Int[]
        for p in frontier, q in neighbors(P, p; connectivity)
            _member(positions, q) || push!(new, q)
        end
        unique!(sort!(new))
        new = _merge(:setdiff, new, out)
        isempty(new) && break
        out = _merge(:union, out, new)
        frontier = new
    end
    return out
end
halo(s::Subset, k::Integer=1; connectivity::Connectivity=Vertex()) = halo(s.parent, s.cells, k; connectivity)
halo(grid, k::Integer=1; connectivity::Connectivity=Vertex()) = _halo(gridrelation(grid), grid, k, connectivity)
_halo(::ParentGrid, grid, k, conn) = halo(parent(grid), parentpositions(grid), k; connectivity=conn)
_halo(::GridRelation, grid, k, conn) = Int[]

"""
    border(parent, positions; connectivity = Vertex()) -> Vector{Int}
    border(s::Subset; connectivity = Vertex())

Ascending **own** positions (indices into `positions`) of members with a
neighbour outside the set.
"""
function border(P, positions::AbstractVector{<:Integer}; connectivity::Connectivity=Vertex())
    out = Int[]
    for (k, p) in enumerate(positions)
        any(q -> !_member(positions, q), neighbors(P, p; connectivity)) && push!(out, k)
    end
    return out
end
border(s::Subset; connectivity::Connectivity=Vertex()) = border(s.parent, s.cells; connectivity)
border(grid; connectivity::Connectivity=Vertex()) = _border(gridrelation(grid), grid, connectivity)
_border(::ParentGrid, grid, conn) = border(parent(grid), parentpositions(grid); connectivity=conn)
_border(::GridRelation, grid, conn) = Int[]

"""
    interior(parent, positions; connectivity = Vertex()) -> AbstractVector{Int}
    interior(s::Subset; connectivity = Vertex())

Ascending own positions of members whose every neighbour is a member: the
complement of [`border`](@ref) in `1:length(positions)`, lazy (storage is the
border's). A complete grid's interior is `1:ncells(grid)` (a range).
"""
function interior(P, positions::AbstractVector{<:Integer}; connectivity::Connectivity=Vertex())
    return Complement(length(positions), border(P, positions; connectivity))
end
interior(s::Subset; connectivity::Connectivity=Vertex()) = interior(s.parent, s.cells; connectivity)
interior(grid; connectivity::Connectivity=Vertex()) = _interior(gridrelation(grid), grid, connectivity)
_interior(::ParentGrid, grid, conn) = interior(parent(grid), parentpositions(grid); connectivity=conn)
_interior(::GridRelation, grid, conn) = Base.OneTo(ncells(grid))
