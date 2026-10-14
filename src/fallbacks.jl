# ---------------------------------------------------------------------------
# Hooks the GeometryOps extension fills in.
#
# The extension never adds methods to a public verb: every public verb's generic
# method calls one of these hooks, and the extension adds methods to the hook
# that dispatch on the manifold. A provider method on the public verb is always
# strictly more specific than the core's generic one, so provider methods and
# extension methods can never be ambiguous with each other. This is the answer
# to the PRD's deferred "how the extension registers fallback methods".
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
fallback_contains(::Any, grid, i, p) = throw(NoGeometryBackend(:cellat))
fallback_neighbors(::Any, grid, i, conn) = throw(NoGeometryBackend(:neighbors))
fallback_vertex_table(::Any, grid) = throw(NoGeometryBackend(:cell_vertices))
fallback_query(::Any, grid, tree, ext) = throw(NoGeometryBackend(:query))
fallback_point_extent(::Any, p) = throw(NoGeometryBackend(:cellat))

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
