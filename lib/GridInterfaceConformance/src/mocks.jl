"""
    GridInterfaceConformance.Mocks

The suite's own test subjects, and a provider's catalogue of what each family
catches: a conforming hierarchical [`CubeSystem`](@ref) and single-fault
wrappers that each take a correct grid or system and break exactly one law
family.

| mock | family it fails | the one lie |
|---|---|---|
| `BrokenNeighbors(grid)` | neighbors | one cell forgets its last neighbour |
| `BrokenTree(grid)` | tree | a leaf index is emitted twice |
| `BrokenIntent(grid)` | tree | `ForSweeps` is ignored: one-cell leaves |
| `LeakyWrapper(grid)` | tree | a hand-written pass-through forgets to forward `maxneighbors` |
| `BrokenAreas(grid)` | geometry | areas off by a factor |
| `BrokenExtent(grid)` | geometry | `cell_extent` does not cover the cell (the tree stays honest) |
| `BrokenLocation(grid)` | location | `cellat` answers the next position (wrap an analytic grid) |
| `FalseAnalytic(grid)` | location | declares `location_is_analytic` with no `cellat` method |
| `BrokenDual(grid)` | dual | vertex fans reversed (wrap a grid with analytic vertices) |
| `BrokenCovering(sys)` | hierarchy | the first root's children miss one cell |
| `BrokenDescendantRange(sys)` | hierarchy | one position too many at the deepest level |
| `stateful_pair(small, large)` | performance | a verb's bytes grow after a call on the bigger grid |

A wrapper is an `AbstractWrappedGrid` (everything else forwards); `FalseAnalytic`,
`LeakyWrapper` and the level grids of a wrapped system are thin explicit types.
"""
module Mocks

import GridInterface as GrI
using GridInterface: ncells, cellindex, localindex, cell_extent, cell_area, cellat, neighbors, treeify, subtree,
    cell_centroid, best_manifold, cellshape, cell_boundary, Vertex, Connectivity, VertexKey, AbstractWrappedGrid
const Structured = GrI.Structured
import GeometryOpsCore as GOC
using GeometryOps.UnitSpherical: UnitSphericalPoint, SphericalCap, spherical_distance
import Extents

export CubeSystem, CubeCell, CubeGrid
export BrokenNeighbors, BrokenTree, BrokenIntent, LeakyWrapper, BrokenAreas, BrokenExtent, BrokenLocation,
    FalseAnalytic, BrokenDual, BrokenCovering, BrokenDescendantRange, Stateful, stateful_pair

# ===========================================================================
# CubeSystem: the gnomonic cube, six root faces refined 4-way. Ids are
# `CubeCell(level, idx)` with `idx` 0-based, face-major and Morton-ordered
# within a face, so the hierarchy is arithmetic (children of `idx` are
# `4idx .+ (0:3)`, the parent `idx ÷ 4`), subtrees are contiguous
# (`has_sorted_subtrees`) and `descendant_range` is closed-form. The level grid
# is a member of the structured family (six faces with seams), which answers
# its neighbours, vertices, point location and tree. DGG's test mock, in
# GridInterface's terms.
# ===========================================================================

struct CubeCell <: GrI.AbstractCellIndex
    level::Int32
    idx::Int64
end
GrI.level(c::CubeCell) = Int(c.level)
Base.isless(a::CubeCell, b::CubeCell) = isless((a.level, a.idx), (b.level, b.idx))

struct CubeSystem <: GrI.AbstractHierarchicalGridSystem
    maxlevel::Int
end
struct CubeGrid <: Structured.AbstractStructuredGrid
    sys::CubeSystem
    level::Int
end

const _FACES = (((1, 0, 0), (0, 1, 0), (0, 0, 1)), ((-1, 0, 0), (0, 0, 1), (0, 1, 0)),
                ((0, 1, 0), (0, 0, 1), (1, 0, 0)), ((0, -1, 0), (1, 0, 0), (0, 0, 1)),
                ((0, 0, 1), (1, 0, 0), (0, 1, 0)), ((0, 0, -1), (0, 1, 0), (1, 0, 0)))
@inline function _cubept(f, u, v)
    N, U, V = _FACES[f]
    x = N[1] + u * U[1] + v * V[1]; y = N[2] + u * U[2] + v * V[2]; z = N[3] + u * U[3] + v * V[3]
    r = sqrt(x^2 + y^2 + z^2)
    return UnitSphericalPoint(x / r, y / r, z / r)
end
_samept(a, b) = all(abs.(a .- b) .< 1e-9)
# the face and side across each face's side, and whether it runs the other way
const _CUBE_SEAMS = let
    corner(f, c) = _cubept(f, ((-1, 1, 1, -1))[c], ((-1, -1, 1, 1))[c])
    ends(f, s) = (corner(f, s), corner(f, mod1(s + 1, 4)))
    ntuple(f -> ntuple(s -> begin
            p, q = ends(f, s)
            for f2 in 1:6, s2 in 1:4
                (f2 == f && s2 == s) && continue
                p2, q2 = ends(f2, s2)
                _samept(p, q2) && _samept(q, p2) && return (f2, s2, true)
                _samept(p, p2) && _samept(q, q2) && return (f2, s2, false)
            end
            error("cube edge ($f, $s) has no partner")
        end, 4), 6)
end
_sidepos(n, i, j, s) = s == 1 ? i : s == 2 ? j : s == 3 ? n + 1 - i : n + 1 - j
_sidecell(n, s, p) = s == 1 ? (p, 1) : s == 2 ? (n, p) : s == 3 ? (n + 1 - p, n) : (1, n + 1 - p)

function _demorton(code::Int64)
    ix = 0; iy = 0
    for b in 0:30
        ix |= ((code >> (2b)) & 1) << b
        iy |= ((code >> (2b + 1)) & 1) << b
    end
    return ix, iy
end
function _morton(ix::Int, iy::Int)
    code = Int64(0)
    for b in 0:30
        code |= Int64((ix >> b) & 1) << (2b)
        code |= Int64((iy >> b) & 1) << (2b + 1)
    end
    return code
end

_n(g::CubeGrid) = 1 << g.level
_facecells(g::CubeGrid) = 1 << (2 * g.level)

# --- the system
GrI.levels(s::CubeSystem) = 0:s.maxlevel
GrI.rootcells(::CubeSystem) = ntuple(f -> CubeCell(0, f - 1), 6)
Base.parent(::CubeSystem, c::CubeCell) =
    c.level == 0 ? throw(ArgumentError("$c is a root and has no parent")) : CubeCell(c.level - 1, c.idx >> 2)
function GrI.children(s::CubeSystem, c::CubeCell)
    c.level < s.maxlevel || throw(ArgumentError("$c is at maxlevel $(s.maxlevel) and has no children"))
    b = c.idx << 2
    return ntuple(k -> CubeCell(c.level + 1, b + k - 1), 4)
end
function GrI.levelgrid(s::CubeSystem, l::Integer)
    0 <= l <= s.maxlevel || throw(ArgumentError("level $l outside 0:$(s.maxlevel)"))
    return CubeGrid(s, Int(l))
end
GrI.has_sorted_subtrees(::CubeSystem) = true
GrI.has_congruent_refinement(::CubeSystem) = true
function GrI.descendant_range(::CubeSystem, c::CubeCell, l::Integer)
    d = 2 * (Int(l) - Int(c.level))
    d >= 0 || throw(ArgumentError("level $l is above $c"))
    return Int(c.idx << d) + 1:Int((c.idx + 1) << d)
end
# the (u, v) box of a cell on its face
function _box(c::CubeCell)
    n = 1 << Int(c.level)
    f, within = divrem(Int(c.idx), n * n)
    ix, iy = _demorton(Int64(within))
    w = 2.0 / n
    return f + 1, -1.0 + w * ix, -1.0 + w * iy, w
end
# The cap from the cell's gnomonic centre to its farthest corner: convex (radius
# below 90°), so it holds the geodesic quad and every descendant inside it.
function GrI.node_extent(::CubeSystem, c::CubeCell)
    f, u0, v0, w = _box(c)
    centre = _cubept(f, u0 + w / 2, v0 + w / 2)
    r = maximum(spherical_distance(centre, _cubept(f, u, v)) for u in (u0, u0 + w), v in (v0, v0 + w))
    return SphericalCap(centre, r + 1e-12)
end

# --- the level grid: identity and the face hooks
GrI.system(g::CubeGrid) = g.sys
GrI.level(g::CubeGrid) = g.level
GOC.best_manifold(::CubeGrid) = GOC.Spherical(; radius=1.0)
GrI.isglobal(::CubeGrid) = true
function GrI.cellindex(g::CubeGrid, i::Integer)
    1 <= i <= ncells(g) || throw(BoundsError(g, i))
    return CubeCell(g.level, i - 1)
end
GrI.localindex(g::CubeGrid, c::CubeCell) =
    (Int(c.level) == g.level && 0 <= c.idx < ncells(g)) ? Int(c.idx) + 1 : nothing
Structured.nfaces(::CubeGrid) = 6
Structured.facesize(g::CubeGrid, f) = (_n(g), _n(g))
Structured.facevertex(g::CubeGrid, f, a, b) = _cubept(f, -1 + 2(a - 1) / _n(g), -1 + 2(b - 1) / _n(g))
Structured.faceposition(g::CubeGrid, f, i, j) = (f - 1) * _facecells(g) + Int(_morton(i - 1, j - 1)) + 1
function Structured.facecell(g::CubeGrid, k)
    f, within = divrem(k - 1, _facecells(g))
    ix, iy = _demorton(Int64(within))
    return (f + 1, ix + 1, iy + 1)
end
function Structured.facecross(g::CubeGrid, f, i, j, s)
    f2, s2, rev = _CUBE_SEAMS[f][s]
    p = _sidepos(_n(g), i, j, s)
    rev && (p = _n(g) + 1 - p)
    i2, j2 = _sidecell(_n(g), s2, p)
    return (f2, i2, j2, s2)
end

# ===========================================================================
# Single-fault wrappers over a correct grid. Each is an AbstractWrappedGrid
# (every other verb forwards) with one lie.
# ===========================================================================

"""
    BrokenNeighbors(grid, cell = 1)

`cell` forgets its last neighbour: the neighbour still lists it, so symmetry
fails (neighbors family).
"""
struct BrokenNeighbors{G} <: AbstractWrappedGrid
    grid::G
    cell::Int
end
BrokenNeighbors(grid) = BrokenNeighbors(grid, 1)
Base.parent(b::BrokenNeighbors) = b.grid
function GrI.neighbors(b::BrokenNeighbors, i::Integer; connectivity::Connectivity=Vertex())
    ns = collect(Int, neighbors(b.grid, i; connectivity))
    Int(i) == b.cell && !isempty(ns) && pop!(ns)
    return ns
end
# the Index forms tell the same lie (a wrapper forwards them to the parent)
GrI.neighbors(b::BrokenNeighbors, i::Integer, ::GrI.Index{GrI.Local}; connectivity::Connectivity=Vertex()) =
    neighbors(b, i; connectivity)
GrI.neighbors(b::BrokenNeighbors, i::Integer, ::GrI.Index{GrI.Global}; connectivity::Connectivity=Vertex()) =
    neighbors(b, i; connectivity)
GrI.neighbors(b::BrokenNeighbors, i::Integer, ::GrI.Index{Type{T}}; connectivity::Connectivity=Vertex()) where {T} =
    [cellindex(b, j, T) for j in neighbors(b, i; connectivity)]

"""
    BrokenTree(grid)

`treeify` emits position 1 twice (tree family). Wrap a plain grid: the lie is
a subtree over `vcat(1:n, 1)`.
"""
struct BrokenTree{G} <: AbstractWrappedGrid
    grid::G
end
Base.parent(b::BrokenTree) = b.grid
GrI.treeify(b::BrokenTree, ::Integer) = subtree(b.grid, vcat(1:ncells(b.grid), 1))
GrI.treeify(b::BrokenTree, h::GrI.LeafHint) = treeify(b, GrI.leaf_capacity(b, h))

"""
    BrokenIntent(grid)

`ForSweeps` is ignored: its tree has one-cell leaves, not the tree of
`leaf_capacity(grid, ForSweeps())` (tree family). Wrap a grid whose tree
honours capacities (a level grid of a system).
"""
struct BrokenIntent{G} <: AbstractWrappedGrid
    grid::G
end
Base.parent(b::BrokenIntent) = b.grid
GrI.treeify(b::BrokenIntent, ::GrI.ForSweeps) = treeify(b.grid, 1)

"""
    BrokenAreas(grid, factor = 1.1)

Areas off by `factor` (geometry family).
"""
struct BrokenAreas{G} <: AbstractWrappedGrid
    grid::G
    factor::Float64
end
BrokenAreas(grid) = BrokenAreas(grid, 1.1)
Base.parent(b::BrokenAreas) = b.grid
GrI.cell_area(b::BrokenAreas, i::Integer) = b.factor * cell_area(b.grid, i)

"""
    BrokenExtent(grid)

`cell_extent` covers only the middle of the cell (a box shrunk to its inner
half, a cap to half its radius) while the tree stays the parent's honest one,
so only "cell_extent contains the cell" fails (geometry family).
"""
struct BrokenExtent{G} <: AbstractWrappedGrid
    grid::G
end
Base.parent(b::BrokenExtent) = b.grid
GrI.cell_extent(b::BrokenExtent, i::Integer) = _shrink(cell_extent(b.grid, i))
_shrink(e::SphericalCap) = SphericalCap(e.point, e.radius / 2)
function _shrink(e::Extents.Extent)
    x0, x1 = e.X; y0, y1 = e.Y
    return Extents.Extent(X=(x0 + (x1 - x0) / 4, x1 - (x1 - x0) / 4), Y=(y0 + (y1 - y0) / 4, y1 - (y1 - y0) / 4))
end
GrI.treeify(b::BrokenExtent, n::Integer) = treeify(b.grid, n)
GrI.treeify(b::BrokenExtent, h::GrI.LeafHint) = treeify(b.grid, h)

"""
    BrokenLocation(grid)

`cellat` answers the position after the right one (location family). Wrap a
grid with analytic location, so that `GridAndTree` forwards the lie rather
than answering from its tree.
"""
struct BrokenLocation{G} <: AbstractWrappedGrid
    grid::G
end
Base.parent(b::BrokenLocation) = b.grid
function GrI.cellat(b::BrokenLocation, p)
    j = cellat(b.grid, p)
    j === nothing && return nothing
    return j == ncells(b.grid) ? j - 1 : j + 1
end

"""
    BrokenDual(grid)

Vertex fans in reverse (clockwise) order, so the dual's rings wind the wrong
way (dual family). Wrap a grid with analytic vertices.
"""
struct BrokenDual{G} <: AbstractWrappedGrid
    grid::G
end
Base.parent(b::BrokenDual) = b.grid
GrI.vertex_cells(b::BrokenDual, key::VertexKey) = reverse!(collect(Int, GrI.vertex_cells(b.grid, key)))
GrI.vertex_cells(b::BrokenDual, v::Integer) = reverse!(collect(Int, GrI.vertex_cells(b.grid, v)))
# a wrapper forwards `dual` to its parent; the dual must be built over the lie
GrI.dual(b::BrokenDual) = invoke(GrI.dual, Tuple{Any}, b)

"""
    FalseAnalytic(grid)

The required four verbs of `grid` plus `location_is_analytic = true` and no
`cellat` method: the generic fallback answers in place of the declared
analytic location (location family). A thin type, not a wrapper.
"""
struct FalseAnalytic{G}
    grid::G
end
GrI.isgrid(::FalseAnalytic) = true
GrI.ncells(f::FalseAnalytic) = ncells(f.grid)
GOC.best_manifold(f::FalseAnalytic) = best_manifold(f.grid)
GrI.cellshape(f::FalseAnalytic) = cellshape(f.grid)
GrI.cell_boundary(f::FalseAnalytic, i::Integer) = cell_boundary(f.grid, i)
GrI.location_is_analytic(::FalseAnalytic) = true

"""
    LeakyWrapper(grid)

A hand-written pass-through (not an `AbstractWrappedGrid`) that forwards every
verb but `maxneighbors`, which the generic answers with `nothing`: a verb the
parent provides is demoted to a fallback (tree family, wrappers group). Wrap a
grid that declares `maxneighbors`.
"""
struct LeakyWrapper{G}
    grid::G
end
Base.parent(w::LeakyWrapper) = w.grid
GrI.isgrid(::LeakyWrapper) = true
for (f, Ts) in GrI.WRAPPED_VERBS
    f === :maxneighbors && continue
    args = GrI._fwd_args(Ts)
    @eval GrI.$f(w::LeakyWrapper, $(GrI._fwd_sig(Ts)...); kw...) = GrI.$f(w.grid, $(args...); kw...)
end
GrI.reindex(::Type{T}, w::LeakyWrapper, c) where {T} = GrI.reindex(T, w.grid, c)
GrI.cell_boundary!(buf::AbstractVector, w::LeakyWrapper, i::Integer) = GrI.cell_boundary!(buf, w.grid, i)

"""
    Stateful(grid, scratch)
    stateful_pair(small, large) -> (Stateful(small, s), Stateful(large, s))

`cell_extent` allocates a buffer sized to the grid the shared scratch last
saw, so a call on the bigger grid costs far more than one on the smaller: the
laziness law fails (performance family), and so does the allocation-free
extent sweep.
"""
mutable struct Scratch
    seen::Int
end
struct Stateful{G} <: AbstractWrappedGrid
    grid::G
    scratch::Scratch
end
Stateful(grid) = Stateful(grid, Scratch(0))
stateful_pair(small, large) = (s = Scratch(0); (Stateful(small, s), Stateful(large, s)))
Base.parent(s::Stateful) = s.grid
function GrI.cell_extent(s::Stateful, i::Integer)
    s.scratch.seen = ncells(s.grid)
    buf = Vector{UInt8}(undef, 64 * s.scratch.seen)
    buf[1] = 0x01
    return cell_extent(s.grid, i)
end

# ===========================================================================
# Single-fault wrappers over a correct system: every system verb forwards, the
# level grids are thin pass-throughs that answer the wrapped system.
# ===========================================================================

abstract type WrappedSystem <: GrI.AbstractHierarchicalGridSystem end
for f in (:levels, :rootcells, :root_index, :has_sorted_subtrees, :has_congruent_refinement)
    @eval GrI.$f(w::WrappedSystem) = GrI.$f(w.sys)
end
GrI.children(w::WrappedSystem, c) = GrI.children(w.sys, c)
Base.parent(w::WrappedSystem, c) = Base.parent(w.sys, c)
GrI.node_extent(w::WrappedSystem, c) = GrI.node_extent(w.sys, c)
GrI.descendant_range(w::WrappedSystem, c, l::Integer) = GrI.descendant_range(w.sys, c, l)
GrI.levelgrid(w::WrappedSystem, l::Integer) = SysGrid(GrI.levelgrid(w.sys, l), w)

struct SysGrid{G,S}
    grid::G
    sys::S
end
GrI.isgrid(::SysGrid) = true
GrI.system(g::SysGrid) = g.sys
for (f, Ts) in GrI.WRAPPED_VERBS
    f === :system && continue
    args = GrI._fwd_args(Ts)
    @eval GrI.$f(g::SysGrid, $(GrI._fwd_sig(Ts)...); kw...) = GrI.$f(g.grid, $(args...); kw...)
end
GrI.reindex(::Type{T}, g::SysGrid, c) where {T} = GrI.reindex(T, g.grid, c)
GrI.cell_boundary!(buf::AbstractVector, g::SysGrid, i::Integer) = GrI.cell_boundary!(buf, g.grid, i)

"""
    BrokenCovering(sys)

The first root's children miss their last cell: the children of level 0 no
longer partition level 1 (hierarchy family).
"""
struct BrokenCovering{S} <: WrappedSystem
    sys::S
end
function GrI.children(w::BrokenCovering, c)
    ch = GrI.children(w.sys, c)
    c == first(GrI.rootcells(w.sys)) || return ch
    return collect(ch)[1:end-1]
end

"""
    BrokenDescendantRange(sys)

`descendant_range` of a non-root cell to the deepest level is one position
too wide (hierarchy family: the children-recursion oracle and the sibling
partition disagree with it).
"""
struct BrokenDescendantRange{S} <: WrappedSystem
    sys::S
end
function GrI.descendant_range(w::BrokenDescendantRange, c, l::Integer)
    r = GrI.descendant_range(w.sys, c, l)
    GrI.level(c) >= 1 && l == last(GrI.levels(w.sys)) || return r
    return first(r):last(r)+1
end

end # module Mocks
