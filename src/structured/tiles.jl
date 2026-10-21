# ===========================================================================
# Tiles: the storage partition of a structured grid is a lazy vector of
# (ti × tj) blocks per face, each a `Subset` whose positions stay in their
# type, so its extent, its subtree in the grid and its own tree are O(1).
# ===========================================================================

"""
    TileIndices(grid, face, ir, jr)

The positions of the `ir × jr` block of face `face`, lazily, `i` fastest:
strictly ascending when the face's positions are column-major (the three
concrete grids), so `Subset(grid, TileIndices(...))` is a valid subset that
keeps its block in its type.
"""
struct TileIndices{G} <: AbstractVector{Int}
    grid::G
    face::Int
    ir::UnitRange{Int}
    jr::UnitRange{Int}
end
Base.size(t::TileIndices) = (length(t.ir) * length(t.jr),)
Base.IndexStyle(::Type{<:TileIndices}) = IndexLinear()
@inline function Base.getindex(t::TileIndices, k::Int)
    @boundscheck checkbounds(t, k)
    a, b = divrem(k - 1, length(t.ir))
    return faceposition(t.grid, t.face, first(t.ir) + b, first(t.jr) + a)
end
# the block's corners directly, without the divrem (the Subset constructor's range check)
Base.first(t::TileIndices) = faceposition(t.grid, t.face, first(t.ir), first(t.jr))
Base.last(t::TileIndices) = faceposition(t.grid, t.face, last(t.ir), last(t.jr))
GrI._strictly_ascending(::TileIndices) = true
function Base.in(p::Integer, t::TileIndices)
    1 <= p <= ncells(t.grid) || return false
    f, i, j = facecell(t.grid, Int(p))
    return f == t.face && i in t.ir && j in t.jr
end

# the same block numbered in its own column-major order (a tile's own tree)
struct LocalTileIndices <: AbstractVector{Int}
    ir::UnitRange{Int}
    jr::UnitRange{Int}
    i0::Int
    j0::Int
    ni::Int
end
Base.size(t::LocalTileIndices) = (length(t.ir) * length(t.jr),)
Base.IndexStyle(::Type{LocalTileIndices}) = IndexLinear()
@inline function Base.getindex(t::LocalTileIndices, k::Int)
    @boundscheck checkbounds(t, k)
    a, b = divrem(k - 1, length(t.ir))
    return (first(t.ir) + b - t.i0 + 1) + (first(t.jr) + a - t.j0) * t.ni
end

const TileSubset{G} = Subset{G,TileIndices{G}}

"""
    TileChunks(grid, (ti, tj) = (64, 64))

Lazy vector of the `ti × tj` tiles of every face (edge tiles smaller),
face-major and column-major over the tiles of a face; what `chunks(grid)`
returns for a structured grid. `chunkat` is O(1).
"""
struct TileChunks{G} <: AbstractVector{Subset{G}}
    grid::G
    ti::Int
    tj::Int
    offs::Vector{Int}        # cumulative tile counts per face
end
function TileChunks(g, tile::NTuple{2,Int}=(64, 64))
    (tile[1] > 0 && tile[2] > 0) || throw(ArgumentError("tile sides must be positive; got $tile"))
    offs = Int[]; s = 0
    for f in 1:nfaces(g)
        ni, nj = facesize(g, f)
        s += cld(ni, tile[1]) * cld(nj, tile[2])
        push!(offs, s)
    end
    return TileChunks{typeof(g)}(g, tile[1], tile[2], offs)
end
Base.size(c::TileChunks) = (isempty(c.offs) ? 0 : last(c.offs),)
Base.IndexStyle(::Type{<:TileChunks}) = IndexLinear()
function Base.getindex(c::TileChunks, k::Int)
    @boundscheck checkbounds(c, k)
    f = length(c.offs) == 1 ? 1 : searchsortedfirst(c.offs, k)
    k0 = k - (f == 1 ? 0 : c.offs[f-1])
    ni, nj = facesize(c.grid, f)
    a, b = divrem(k0 - 1, cld(ni, c.ti))
    ir = (b * c.ti + 1):min((b + 1) * c.ti, ni)
    jr = (a * c.tj + 1):min((a + 1) * c.tj, nj)
    t = TileIndices(c.grid, f, ir, jr)
    facecolumnmajor(c.grid, f) || return Subset(c.grid, sort!(collect(t)))     # positions not in lattice order
    return Subset(c.grid, t)
end
function GrI.chunkat(c::TileChunks, p::Integer)
    f, i, j = _facecell(c.grid, p)
    ni, _ = facesize(c.grid, f)
    return (f == 1 ? 0 : c.offs[f-1]) + (i - 1) ÷ c.ti + 1 + ((j - 1) ÷ c.tj) * cld(ni, c.ti)
end

GrI.chunks(g::AbstractStructuredGrid) = TileChunks(g)
GrI.chunks(g::AbstractStructuredGrid, tile::NTuple{2,Int}) = TileChunks(g, tile)

# a tile's extent, its subtree in the grid and its own tree are the block's
GrI.grid_extent(s::TileSubset) = (t = s.cells; blockextent(t.grid, t.face, t.ir, t.jr))
GrI.subtree(g::AbstractStructuredGrid, t::TileIndices; leaf_capacity=Balanced()) =
    TopDownQuadtreeCursor(g, t.face, t.ir, t.jr; leafsize=_leafsize(GrI._capacity(g, leaf_capacity)))
GrI.subtree(g::AbstractStructuredGrid, s::TileSubset; leaf_capacity=Balanced()) = subtree(g, s.cells; leaf_capacity)
function GrI.treeify(s::TileSubset, n::Integer)
    t = s.cells
    return TopDownQuadtreeCursor(t.grid, t.face, t.ir, t.jr; leafsize=_leafsize(n), local_positions=true)
end
