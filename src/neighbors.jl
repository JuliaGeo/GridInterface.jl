# ===========================================================================
# Neighbours: every single-ring grid has them; a trait says whether they are
# analytic or derived. DGG's contract, lifted as is.
# ===========================================================================

abstract type Connectivity end
"Cells sharing at least a vertex (Moore). The default."
struct Vertex <: Connectivity end
"Cells sharing a whole edge (von Neumann)."
struct Edge <: Connectivity end

abstract type Winding end
struct CounterClockwise <: Winding end
struct Clockwise <: Winding end
struct CustomOrder <: Winding end
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

Capacity hint for static stencil buffers; `nothing` means unbounded/undeclared.
"""
maxneighbors(grid, ::Connectivity) = nothing
maxneighbors(grid) = maxneighbors(grid, Vertex())

"""
    winding(grid, connectivity = Vertex()) -> Winding

The order [`neighbors`](@ref) returns: `CounterClockwise()`, `Clockwise()`,
`CustomOrder()` (deterministic, documented) or `Unordered()` (default).
"""
winding(grid, ::Connectivity) = Unordered()
winding(grid) = winding(grid, Vertex())

@noinline _refuse_neighbors(grid) = throw(ArgumentError(
    "neighbours exist only for SingleRing grids (or grids that define `neighbors`); " *
    "$(typeof(grid)) declares $(cellshape(grid))"))

"""
    neighbors(grid, i; connectivity = Vertex()) -> AbstractVector{Int}

Positions of the cells sharing an edge (`Edge()`) or at least a vertex
(`Vertex()`) with cell `i`, excluding `i`. Analytic grids override it; any
other `SingleRing` grid gets the geometric fallback (extension).
"""
neighbors(grid, i::Integer; connectivity::Connectivity=Vertex()) =
    _neighbors(cellshape(grid), grid, i, connectivity)
_neighbors(::SingleRing, grid, i, conn) = fallback_neighbors(best_manifold(grid), grid, i, conn)
_neighbors(::CellShape, grid, i, conn) = _refuse_neighbors(grid)

# Subsets clip the parent's answer to membership (distance measured in the parent).
function neighbors(s::Subset, i::Integer; connectivity::Connectivity=Vertex())
    out = Int[]
    for j in neighbors(s.parent, _ppos(s, i); connectivity)
        k = _local_of_parent(s, j)
        k === nothing || push!(out, k)
    end
    return out
end
maxneighbors(s::Subset, c::Connectivity) = maxneighbors(s.parent, c)
winding(s::Subset, c::Connectivity) = winding(s.parent, c)

"""
    neighbors(grid, i, sel::Index; connectivity)

Neighbours of position `i` in the index space `sel` selects
(`Index(Local())`, `Index(Global())`, `Index(T)`).
"""
neighbors(grid, i::Integer, sel::Index; connectivity::Connectivity=Vertex()) =
    [toindex(grid, j, sel) for j in neighbors(grid, i; connectivity)]

"""
    ring(grid, i, k; connectivity = Vertex()) -> Vector{Int}

Positions at adjacency distance exactly `k` from `i` (`ring(grid, i, 0) == [i]`).
Discs are rings concatenated. Generic method: breadth-first over
[`neighbors`](@ref), keeping each ring in discovery order (so ring 1 has the
order of `neighbors`).
"""
function ring(grid, i::Integer, k::Integer; connectivity::Connectivity=Vertex())
    k >= 0 || throw(ArgumentError("ring distance must be non-negative"))
    k == 0 && return Int[Int(i)]
    seen = Set{Int}((Int(i),))
    frontier = Int[Int(i)]
    for _ in 1:k
        next = Int[]
        for c in frontier, j in neighbors(grid, c; connectivity)
            j in seen && continue
            push!(seen, j)
            push!(next, j)
        end
        frontier = next
    end
    return frontier
end

# ---------------------------------------------------------------------------
# Region verbs over subsets. Generic over `neighbors` and membership; a grid
# with structure overrides them.
# ---------------------------------------------------------------------------

"""
    halo(s, k = 1; connectivity = Vertex()) -> Vector{Int}

Ascending **parent** positions of the cells outside subset `s` within `k`
steps of a member. A complete grid has no outside: empty.
"""
halo(grid, k::Integer=1; connectivity::Connectivity=Vertex()) = Int[]
function halo(s::Subset, k::Integer=1; connectivity::Connectivity=Vertex())
    P = s.parent
    out = Set{Int}()
    frontier = collect(s.cells)
    seen = Set{Int}(s.cells)
    for _ in 1:k
        next = Int[]
        for p in frontier, q in neighbors(P, p; connectivity)
            q in seen && continue
            push!(seen, q); push!(next, q); push!(out, q)
        end
        frontier = next
    end
    return sort!(collect(out))
end

"""
    border(s; connectivity = Vertex()) -> Vector{Int}

Ascending own positions of members with a neighbour outside `s`.
"""
border(grid; connectivity::Connectivity=Vertex()) = Int[]
function border(s::Subset; connectivity::Connectivity=Vertex())
    out = Int[]
    for (k, p) in enumerate(s.cells)
        any(q -> _local_of_parent(s, q) === nothing, neighbors(s.parent, p; connectivity)) && push!(out, k)
    end
    return out
end

"""
    interior(s; connectivity = Vertex()) -> Vector{Int}

Ascending own positions of members whose every neighbour is a member.
"""
interior(grid; connectivity::Connectivity=Vertex()) = collect(1:ncells(grid))
function interior(s::Subset; connectivity::Connectivity=Vertex())
    b = border(s; connectivity)
    return setdiff(1:ncells(s), b)
end
