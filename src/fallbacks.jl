# ---------------------------------------------------------------------------
# Hooks the GeometryOps extension fills in.
#
# The extension never adds methods to a public verb: every public verb's generic
# method calls one of these hooks, and the extension adds methods to the hook
# that dispatch on the manifold. A provider method on the public verb is always
# strictly more specific than the core's generic one, so provider methods and
# extension methods can never be ambiguous with each other.
# ---------------------------------------------------------------------------

"""
    NoGeometryBackend

Thrown when a fallback that needs real geometry is called without GeometryOps
loaded.
"""
struct NoGeometryBackend <: Exception
    verb::Symbol
end
Base.showerror(io::IO, e::NoGeometryBackend) = print(io,
    "GridInterface: the generic `", e.verb, "` needs real geometry. Either the grid ",
    "defines `", e.verb, "` itself, or `import GeometryOps` to activate ",
    "GridInterfaceGeometryOpsExt.")

# Each hook takes the manifold first so the extension's methods are strictly
# more specific than these stubs.
fallback_centroid(::Any, grid, i) = throw(NoGeometryBackend(:cell_centroid))
fallback_area(::Any, grid, i) = throw(NoGeometryBackend(:cell_area))
fallback_extent(::Any, grid, i) = throw(NoGeometryBackend(:cell_extent))
fallback_treeify(::Any, grid) = throw(NoGeometryBackend(:treeify))
fallback_subtree(::Any, grid, cells) = throw(NoGeometryBackend(:subtree))
fallback_cellat(::Any, grid, p) = throw(NoGeometryBackend(:cellat))
fallback_cellsat(::Any, grid, p) = throw(NoGeometryBackend(:cellsat))
fallback_contains(::Any, grid, i, p) = throw(NoGeometryBackend(:cell_contains))
fallback_neighbors(::Any, grid, i, conn) = throw(NoGeometryBackend(:neighbors))
fallback_vertex_table(::Any, grid) = throw(NoGeometryBackend(:cell_vertices))
fallback_query(::Any, grid, pred) = throw(NoGeometryBackend(:query))
fallback_point_extent(::Any, p) = throw(NoGeometryBackend(:cellat))
fallback_region_extent(::Any, region) = throw(NoGeometryBackend(:query))

"""
    has_geometry_backend() -> Bool

Whether GridInterfaceGeometryOpsExt is loaded.
"""
has_geometry_backend() = _HAS_BACKEND[]
const _HAS_BACKEND = Ref(false)

"""
    clear_caches!()

Drop the per-grid trees and vertex tables the geometric fallbacks cache.
"""
function clear_caches! end

# ---------------------------------------------------------------------------
# Spherical type contract. On a spherical or ellipsoidal manifold every point a
# grid hands out or a consumer hands in is a `UnitSphericalPoint`, and every
# cell or node extent is a `SphericalCap`. A plain 3-tuple is silently read as
# (lon, lat, z) by GeometryOps and refused by ConservativeRegridding's
# clipper, so the core refuses it up front, here, at the fallback boundary.
# ---------------------------------------------------------------------------

"""
    UnitSphericalPoint

The point type of every spherical grid: a unit 3-vector (GeometryOps'
`UnitSpherical.UnitSphericalPoint`). Boundaries, corners, centroids and query
points on a `Spherical` manifold must be these; convert lon/lat with
`GeometryOps.UnitSphereFromGeographic()`.
"""
const UnitSphericalPoint = UnitSpherical.UnitSphericalPoint
"""
    SphericalCap

The extent type of every spherical grid (GeometryOps' `UnitSpherical.SphericalCap`):
a centre `UnitSphericalPoint` and an angular radius. `cell_extent`, `node_extent`
and chunk extents on a `Spherical` manifold must be caps; `Extents.union` and
`Extents.intersects` are defined on them.
"""
const SphericalCap = UnitSpherical.SphericalCap

"Whether manifold `m` is curved (spherical or ellipsoidal): the unit-sphere types are required."
iscurved(::Spherical) = true
iscurved(::GOC.Geodesic) = true
iscurved(::Any) = false

@noinline _throw_point_type(m, p) = throw(ArgumentError(
    "GridInterface: on $(nameof(typeof(m))) manifolds points must be " *
    "`UnitSphericalPoint`s (GeometryOps.UnitSpherical), got $(typeof(p)). " *
    "Convert lon/lat with `GeometryOps.UnitSphereFromGeographic()((lon, lat))`."))
@noinline _throw_extent_type(m, e) = throw(ArgumentError(
    "GridInterface: on $(nameof(typeof(m))) manifolds extents must be " *
    "`SphericalCap`s (GeometryOps.UnitSpherical), got $(typeof(e))."))

"""
    check_point(manifold, p) -> p

`p` unchanged when its type is allowed on `manifold`; an `ArgumentError`
naming the conversion otherwise. Free on the plane; a type check on the sphere.
"""
@inline check_point(m, p) = iscurved(m) ? _check_usp(m, p) : p
@inline _check_usp(m, p::UnitSphericalPoint) = p
_check_usp(m, p) = _throw_point_type(m, p)

"""
    check_extent(manifold, e) -> e

`e` unchanged when it is a `SphericalCap` on a curved manifold (or anything on
the plane); an `ArgumentError` otherwise.
"""
@inline check_extent(m, e) = iscurved(m) ? _check_cap(m, e) : e
@inline _check_cap(m, e::SphericalCap) = e
_check_cap(m, e) = _throw_extent_type(m, e)

"""
    check_ring(manifold, ring) -> ring

Point-type check for a whole ring: compile-time when the ring's element type
is concrete, one check of the first point otherwise.
"""
@inline function check_ring(m, ring)
    iscurved(m) || return ring
    eltype(ring) <: UnitSphericalPoint && return ring
    isempty(ring) || _check_usp(m, first(ring))
    return ring
end
