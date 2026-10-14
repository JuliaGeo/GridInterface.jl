# GridInterface.jl prototype environment

`prototype/env` develops, by path, GridInterface (this package), its
`lib/SpatialTreeInterfaceShim` and `lib/GridInterfaceConformance`,
ConservativeRegridding (CR), GlobalRegridding and DiscreteGlobalGrids (DGG),
and adds GeometryOps, GeometryOpsCore, Extents, GeoInterface, Healpix,
BenchmarkTools and Test. Neither CR nor DGG is modified: everything that makes
them GridInterface grids lives in `prototype/adapters/`.

| file | what it does |
|---|---|
| `setup_env.jl` | (re)creates `env/` (run once) |
| `adapters/DGGAdapter.jl` | DGG systems, level grids, partial grids and multi-order sets as GI systems/grids |
| `adapters/CRAdapter.jl` | (a) CR's curvilinear grids and `Healpix.HealpixMap` as GI grids; (b) `CRTree`, a Trees-compatible bridge so `CR.Regridder` accepts any GI grid |
| `conformance.jl` | conformance on every adapted grid (about 3 min) |
| `bench.jl` | performance laws and the regrid parity check (about 6 min) |

## Loading

```julia
# once:  julia /home/anshul/geo/gridinterface/GridInterface.jl/prototype/setup_env.jl
# then:  julia --project=/home/anshul/geo/gridinterface/GridInterface.jl/prototype/env
using GridInterface, GridInterfaceConformance
import GridInterface as GrI
include("/home/anshul/geo/gridinterface/GridInterface.jl/prototype/adapters/DGGAdapter.jl"); using .DGGAdapter
include("/home/anshul/geo/gridinterface/GridInterface.jl/prototype/adapters/CRAdapter.jl"); using .CRAdapter
import Healpix; adapt_healpix!(Healpix)   # optional: HealpixMap through CR's Healpix extension
import Extents, GeometryOps as GO
```

The geometry-backed fallbacks (centroid, area, extent, R-tree, point location,
geometric neighbours and vertex merge) live in `ext/GridInterfaceGeometryOpsExt.jl`
and switch on when GeometryOps is loaded; `GrI.has_geometry_backend()` says
whether they are active.

## One example per family

### Identity

```julia
hp = DGG.HEALPixSystem()
g  = levelgrid(hp, 4)                    # a DGG level grid is a GI grid
ncells(g), level(g)                      # (3072, 4)
c  = cellindex(g, 17)                    # the typed id at position 17
localindex(g, c)                         # 17
cellindextypes(g)                        # the id scheme(s) the grid declares
pg = DGG.subtree(hp, rootcells(hp)[2], 4) # a partial grid: ParentGrid relation
gridrelation(pg), globalindex(pg, cellindex(pg, 1))  # (ParentGrid(), position in the complete level grid)
```

### Geometry and location

```julia
cellshape(g), tessellation(g)            # (SimpleCells(), Tessellating())
cell_boundary(g, 17)                     # densified ring (HEALPix: 32 points)
cell_corners(g, 17)                      # the 4 true corners
cell_centroid(g, 17), cell_area(g, 17), cell_extent(g, 17)
p = cell_centroid(g, 17)
location_is_analytic(g), cellat(g, p)    # (true, 17), allocation free
ll = lonlat_grid(36, 18)                 # CR CellBasedGrid on the sphere
cellat(ll, cell_centroid(ll, 100))       # 100, through the cached fallback R-tree
```

### Tree

```julia
t = treeify(g)                           # HierarchyCursor: O(1), 0 bytes
typeof(t), GrI.STI.isspatialtree(t)
grid_extent(g)                           # root node extent
q = query(g, Extents.Extent(X=(0.0, 10.0), Y=(40.0, 50.0)))  # sorted positions; the lon/lat box becomes a cap (query_region)
st = subtree(g, q)                       # tree over just those cells
gt = GridAndTree(ll, treeify(ll))        # pair a grid with a stored tree
```

### Subsets and chunks

```julia
s = subset(g, q)                         # a Subset: positions 1:length(q)
parent(s) === g, parentpositions(s) == q, gridrelation(s)
cell_centroid(s, 1) == cell_centroid(g, q[1])
ch = chunks(levelgrid(hp, 5))            # hierarchy_chunks at level 2 for sorted systems
length(ch), ncells(first(ch)), parentpositions(first(ch))   # each chunk is a contiguous range
```

### Neighbours

```julia
neighbors(g, 17)                         # Vertex connectivity, provider answer
neighbors(g, 17; connectivity=Edge())
maxneighbors(g, Vertex()), winding(g, Vertex()), neighbors_are_analytic(g)
ring(g, 17, 2)                           # all cells within 2 steps
halo(s, 1), border(s), interior(s)       # parent positions around / inside a subset
```

### Dual

```julia
h3 = levelgrid(DGG.H3System(), 1)
conformity(h3)                           # Conforming()
d = dual(h3)                             # one cell per primal vertex
ncells(d), nvertices(h3), cell_vertices(h3, 1), vertex_cells(h3, 1)
cell_boundary(d, 1)                      # centroids of the cells around vertex 1
```

### Hierarchy: window cursor, selection cursor, multi-level set

```julia
levels(hp), rootcells(hp)[1:3], children(hp, rootcells(hp)[1])
has_sorted_subtrees(hp), descendant_range(hp, rootcells(hp)[1], 2)   # 1:16
node_extent(hp, rootcells(hp)[1]), ancestor(hp, c, 1)

a5 = levelgrid(DGG.A5System(), 3)        # unsorted system: selection mode
typeof(treeify(a5)) <: SelectionCursor

set = DGG.query(hp, DGG.MultiOrderCoverage(Extents.Extent(X=(-10.0, 30.0), Y=(35.0, 60.0))); level=8)
mls = dgg_multilevel(set)                # O(1) wrap of DGG's MultiOrderCellSet
ml2 = MultiLevelSet(hp, set.cells)       # or build one from any antichain of ids
GrI.check_antichain(mls), gridrelation(mls)          # (true, ParentSystem())
cell_area(mls, 1), cellat(mls, cell_centroid(mls, 10))  # forwards to levelgrid(sys, level(c))
typeof(treeify(mls)) <: HierarchyCursor              # window mode, O(1)
rows, cols = GrI.aggregation_pairs(mls)              # reference-level cell -> member
```

### Regridding through CR with any GI grid

```julia
dst = levelgrid(DGG.IGeo7System(), 3)
R1 = CR.Regridder(dst, levelgrid(hp, 3); threaded=false)                    # DGG's own path
R2 = CR.Regridder(cr_tree(dst), cr_tree(levelgrid(hp, 3)); threaded=false)  # GI generic cursor
R1.intersections == R2.intersections                                        # true
R3 = CR.Regridder(cr_tree(dst), cr_tree(mls); threaded=false)               # a multi-level set as a source
```

Mind the radius: DGG grids use `Spherical(radius=1)` and CR's `Spherical()`
defaults to Earth's radius, so mixing the two needs the three-argument
`Regridder(manifold, dst, src)`.

### Conformance

```julia
rep = check_grid(g; name="HEALPix level 4")
print_report(rep)                        # per-family pass/fail/skip + metrics + fallbacks answering
family_status(rep)
rep2 = check_system(DGG.H3System(); name="H3", grid_level=2)   # hierarchy + level grid + multi-level set
test_grid(levelgrid(DGG.IGeo7System(), 2))                     # the same laws as @test inside a testset
```
