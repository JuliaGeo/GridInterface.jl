# GridInterface.jl prototype environment (revision 2)

`prototype/env` develops, by path, GridInterface (this package), its
`lib/SpatialTreeInterfaceShim` and `lib/GridInterfaceConformance`,
ConservativeRegridding (CR), GlobalRegridding and DiscreteGlobalGrids (DGG),
and adds GeometryOps, GeometryOpsCore, Extents, GeoInterface, Healpix,
BenchmarkTools and Test. Neither CR nor DGG is modified: everything that makes
them GridInterface grids lives in `prototype/adapters/`.

| file | what it does |
|---|---|
| `setup_env.jl` | (re)creates `env/` (run once) |
| `adapters/DGGAdapter.jl` | DGG systems, level grids, partial grids and multi-order sets as GI systems/grids; DGG's own trees; aligned HEALPix neighbours; analytic HEALPix vertices |
| `adapters/CRAdapter.jl` | (a) CR's curvilinear grids (with raster-tile chunks) and `Healpix.HealpixMap` as GI grids; (b) `CRTree`, the bridge that lets `CR.Regridder` consume any `GridAndTree` |
| `conformance.jl` | conformance on every adapted grid (about 3 min) |
| `bench.jl` | round-one performance laws and regrid parity (about 6 min) |
| `../test/revision_bench.jl` | revision benchmarks: analytic vertices, root index, split_weight, lazy chunks, halo, queries, four-method regrid |

## Loading

```julia
# once:  julia /home/anshul/geo/gridinterface/GridInterface.jl/prototype/setup_env.jl
# then:  julia --project=/home/anshul/geo/gridinterface/GridInterface.jl/prototype/env
using GridInterface, GridInterfaceConformance
import GridInterface as GrI
include("/home/anshul/geo/gridinterface/GridInterface.jl/prototype/adapters/DGGAdapter.jl"); using .DGGAdapter
include("/home/anshul/geo/gridinterface/GridInterface.jl/prototype/adapters/CRAdapter.jl"); using .CRAdapter
import Healpix; adapt_healpix!(Healpix)   # optional: HealpixMap through CR's Healpix extension
import Extents, GeometryOps as GO, GeometryOpsCore as GOC
```

The geometry-backed fallbacks (centroid, area, extent, R-tree, point location,
geometric neighbours, vertex merge, query predicates) live in
`ext/GridInterfaceGeometryOpsExt.jl` and switch on when GeometryOps is loaded.

## What changed since round one (in one screen)

- **Required set is four methods**: `ncells`, `best_manifold`, `cellshape`, and `cell_boundary` (single ring) or `cell_geometry`. Such a grid passes `check_grid` and regrids through CR (`test/grids.jl`: `MinLonLat`).
- **Spherical types are enforced**: on a `Spherical` manifold points are `UnitSphericalPoint` and extents `SphericalCap`; the core refuses anything else with an `ArgumentError` naming `GO.UnitSphereFromGeographic()`, and conformance fails it first.
- **Engines take a grid paired with a tree** (`GridAndTree`). `split_weight` forwards through every wrapper and the CR bridge. A provider's `treeify` wins (DGG grids use DGG's cursors); the generic `hierarchy_tree` is the default for other level grids.
- **Neighbours**: disc `neighbors(g, i, k)`, `maxneighbors(g, k, conn)`, `maxring`; on subsets the plain form is clipped (local), `Index(Local())` is unclipped with `OUTSIDE` (0), `Index(Global())` gives complete-grid positions. `CounterClockwise` now promises alignment: Edge neighbour k lies across edge (corner k, corner k+1).
- **Collections**: local / parent (`parentpositions`) / complete (`globalindex`) spaces; `union`, `intersect`, `setdiff`, `grow`; `cellsat`; `query(g, Intersects|Within|Contains|CentroidIn(region))`.
- **Dual**: vertices are points where ≥ 3 cells meet (or boundary corners); `cell_vertices` ⊆ `cell_corners`; `vertex_isopen`, `dual_ring_cells`; typed `VertexKey`s with derived dense numbering; HEALPix vertices are analytic.
- **Chunks**: `chunks(g)` is the storage partition (lazy), `chunks(g, level)` on hierarchies, `chunkat`, cheap covering extents; chunks of subsets and multi-level sets; CR rasters chunk into `TileIndices` subsets.

## One example per family

### Identity and the three index spaces

```julia
hp = DGG.HEALPixSystem()
g  = levelgrid(hp, 4)                    # a DGG level grid is a GI grid
ncells(g), level(g), system(g) === hp    # (3072, 4, true)
c  = cellindex(g, 17); localindex(g, c)  # typed id at position 17; back to 17
pg = DGG.subtree(hp, rootcells(hp)[2], 4) # a partial grid: ParentGrid relation
s  = subset(pg, 10:20)                   # a subset of a subset
parentpositions(s)[1], globalindex(s, cellindex(s, 1))   # (10, position in the complete level grid)
```

### Geometry and location

```julia
cellshape(g), edgekind(g), tessellation(g)   # (SimpleCells(), CurvedEdges(), Tessellating())
cell_boundary(g, 17)                     # densified ring (HEALPix: 32 points)
cell_corners(g, 17)                      # the 4 true corners (default: the ring with repeats collapsed)
buf = UnitSphericalPoint{Float64}[]; cell_boundary!(buf, g, 17)   # reuse a buffer
cellsize(g), isglobal(g), cell_geometry_is_expensive(g)
p = cell_centroid(g, 17); cellat(g, p)   # 17, analytic, 0 bytes
ll = lonlat_grid(36, 18)                 # CR CellBasedGrid on the sphere
cellat(ll, cell_centroid(ll, 100))       # 100, through the cached fallback R-tree
cellsat(ll, cell_corners(ll, 100)[3])    # every cell containing a point (4 at a corner)
```

### Trees: provider first, GridAndTree for engines

```julia
treeify(g)                               # DGG's HierarchicalGridCursor (provider wins)
gt = GridAndTree(g)                      # pair with treeify(g); or GridAndTree(g, GrI.hierarchy_tree(g))
GrI.split_weight(gt) == ncells(g)        # forwarded
leaf_capacity(g), root_index(DGG.CopernicusDEMSystem(90)) isa RootIndex
```

### Queries

```julia
box = Extents.Extent(X=(-10.0, 30.0), Y=(35.0, 60.0))   # read on the sphere as DGG does:
query(g, Intersects(box))                # boundary sampled every ≤ 2°, joined by great circles
query(g, CentroidIn(box)), query(g, Within(box))
query(g, box)                            # == query(g, Intersects(box))
query_candidates(g, box)                 # the tree's candidates (cap test only)
```

### Subsets, algebra and chunks

```julia
a = subset(g, 100:200); b = subset(g, 150:300)
parentpositions(union(a, b)), parentpositions(intersect(a, b)), parentpositions(setdiff(a, b))
grow(a, 1)                               # a plus its halo
ch = chunks(levelgrid(hp, 8))            # lazy: 192 chunks of 4096 cells, 0 bytes
chunkat(ch, 300_000), grid_extent(ch[5]) # O(1); a cap from the chunk's root
chunks(levelgrid(hp, 8), 2)              # choose the level
chunks(subset(levelgrid(hp, 8), 1000:90_000), 3)   # chunks of a subset: windows
tc = chunks(lonlat_grid(1440, 720))      # CR raster: TileChunks of TileIndices subsets
subtree(parent(tc[3]), tc[3]), treeify(tc[3])      # O(1) quadtree cursors
```

### Neighbours

```julia
neighbors(g, 17)                         # Vertex, CCW, starting across edge (corner 1, corner 2)
neighbors(g, 17; connectivity=Edge())    # Edge neighbour k across edge k
neighbors(g, 17, 2)                      # the disc: rings 1:2
maxneighbors(g, 2, Vertex()), maxring(g, 1), winding(g)
s = subset(g, 100:200)
neighbors(s, 1)                          # clipped: local positions of members
neighbors(s, 1, Index(Local()))          # unclipped: OUTSIDE (0) for non-members
neighbors(s, 1, Index(Global()))         # unclipped: complete-grid positions
halo(s), border(s), interior(s)          # parent positions around; own positions inside (lazy)
```

### Dual and vertices

```julia
nvertices(g), vertices_are_analytic(g)   # (3074, true): HEALPix names vertices itself
k = vertex_key(g, 17, 2)                 # VertexKey(owner cell, corner)
vertex_cells(g, k), vertex_point(g, k), vertex_index(g, k)
cell_vertices(g, 17), cell_vertex_keys(g, 17), vertex_isopen(g, 5)
d = dual(g); cellat(d, vertex_point(g, 5))   # dual location through the primal host
h3 = levelgrid(DGG.H3System(), 1)        # odd level: two-cell bend points are not vertices
nvertices(h3), dual_ring_cells(dual(h3), 1)
```

### Hierarchy and multi-level sets

```julia
levels(hp), rootcells(hp)[1:3], children(hp, rootcells(hp)[1])
descendant_range(hp, rootcells(hp)[1], 2), ancestor(hp, c, 1)
set = DGG.query(hp, DGG.MultiOrderCoverage(box); level=8)
mls = dgg_multilevel(set)                # O(1) wrap of DGG's MultiOrderCellSet
chunks(mls, 4)                           # a member coarser than level 4 is its own chunk
treeify(subset(mls, 1:50))               # window mode on a subset of a set
em = expand_to_reference(mls)            # behaviour 1: expand, regrid, aggregate (aggregation_pairs)
rows, cols = aggregation_pairs(mls)      # behaviour 2 (direct) is Regridder(dst, mls)
```

### Regridding through CR

```julia
dst = levelgrid(DGG.IGeo7System(), 3)
R1 = CR.Regridder(dst, levelgrid(hp, 3); threaded=false)                 # DGG's own path
R2 = regridder(dst, GridAndTree(levelgrid(hp, 3), GrI.hierarchy_tree(levelgrid(hp, 3))); threaded=false)
R1.intersections == R2.intersections     # true: the generic cursor through the bridge
R3 = regridder(dst, mls; threaded=false) # a multi-level set as a source (direct form)
```

Mind the radius: the manifold's radius is the grid's area unit. DGG grids use
`Spherical(radius=1)` and CR's `Spherical()` defaults to Earth's radius; a
two-argument consumer never reconciles them, so mix them through the
three-argument `Regridder(manifold, dst, src)`.

### Conformance

```julia
rep = check_grid(g; name="HEALPix level 4"); print_report(rep); family_status(rep)
rep2 = check_system(DGG.H3System(); name="H3", grid_level=1)  # hierarchy + level grid + multi-level set
```

### A four-method grid

```julia
struct LL <: GrI.AbstractGrid; n::Int; end
GrI.ncells(g::LL) = 2g.n^2
GOC.best_manifold(::LL) = GOC.Spherical(; radius=1.0)
GrI.cellshape(::LL) = GrI.SimpleCells()
function GrI.cell_boundary(g::LL, i::Integer)
    a = (i - 1) % 2g.n; b = (i - 1) ÷ 2g.n; d = 180 / g.n; f = GO.UnitSphereFromGeographic()
    x0, y0 = -180 + a * d, -90 + b * d
    return [f((x0, y0)), f((x0 + d, y0)), f((x0 + d, y0 + d)), f((x0, y0 + d))]
end
family_status(check_grid(LL(6)))         # passes every family
regridder(levelgrid(hp, 2), LL(6); threaded=false)
```
