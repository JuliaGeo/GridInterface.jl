```@meta
CurrentModule = GridInterface
```

# API

```@index
Modules = [GridInterface]
```

```@docs
GridInterface
```

## Re-exports

`SmallVector` (SmallCollections), `UnitSphericalPoint` and `SphericalCap` (UnitSphericalGeometryCore) are re-exported so that providers need no separate dependency: small answers (neighbour rings, vertex fans, corners) are `SmallVector`s, and spherical points and extents are these two types.

## Internals

The GeometryOps extension fills in these hooks; the lazy vectors are what verbs return where they would otherwise materialise O(ncells) positions.

```@docs
NoGeometryBackend
iscurved
check_point
check_extent
check_ring
MappedVec
FuncVec
Unsentinel
RangeCat
Complement
ClosedRing
CollapsedRing
ring_polygon
```
