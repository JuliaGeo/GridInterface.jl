# --- sampling ----------------------------------------------------------------
sample_positions(n, k) = n <= k ? collect(1:n) : unique(round.(Int, range(1, n; length=k)))

# --- points on either manifold ------------------------------------------------
_v(p) = (Float64(GI.x(p)), Float64(GI.y(p)), Float64(GI.z(p)))
_xy(p) = (Float64(GI.x(p)), Float64(GI.y(p)))
_dot(a, b) = a[1] * b[1] + a[2] * b[2] + a[3] * b[3]
_cross(a, b) = (a[2] * b[3] - a[3] * b[2], a[3] * b[1] - a[1] * b[3], a[1] * b[2] - a[2] * b[1])
_norm(a) = sqrt(_dot(a, a))
_usp(a) = (n = _norm(a); UnitSphericalPoint(a[1] / n, a[2] / n, a[3] / n))

# point a fraction t of the way from p to q (geodesic on the sphere)
_toward(::Planar, p, q, t) = (a = _xy(p); b = _xy(q); (a[1] + t * (b[1] - a[1]), a[2] + t * (b[2] - a[2])))
_toward(::Spherical, p, q, t) = (a = _v(p); b = _v(q); _usp(((1 - t) * a[1] + t * b[1], (1 - t) * a[2] + t * b[2], (1 - t) * a[3] + t * b[3])))

_point_extent(m, p) = GrI.fallback_point_extent(m, p)

# Does extent `e` contain point `p`? (closed; caps by distance, boxes by bounds)
function extent_contains_point(::Spherical, e::SphericalCap, p; slack=1e-9)
    e.radius >= π && return true
    return spherical_distance(e.point, _usp(_v(p))) <= e.radius + slack
end
function extent_contains_point(::Spherical, e::Extents.Extent, p; slack=1e-9)
    v = _v(p)
    return e.X[1] - slack <= v[1] <= e.X[2] + slack && e.Y[1] - slack <= v[2] <= e.Y[2] + slack &&
           e.Z[1] - slack <= v[3] <= e.Z[2] + slack
end
function extent_contains_point(::Planar, e::Extents.Extent, p; slack=1e-9)
    x, y = _xy(p)
    s = slack * max(1.0, abs(x), abs(y))
    return e.X[1] - s <= x <= e.X[2] + s && e.Y[1] - s <= y <= e.Y[2] + s
end

function extent_contains_extent(::Spherical, big::SphericalCap, small::SphericalCap; slack=1e-9)
    big.radius >= π && return true
    small.radius >= π && return false
    return spherical_distance(big.point, small.point) + small.radius <= big.radius + slack
end
extent_contains_extent(::Planar, big::Extents.Extent, small::Extents.Extent; slack=1e-9) =
    big.X[1] - slack <= small.X[1] && small.X[2] <= big.X[2] + slack &&
    big.Y[1] - slack <= small.Y[1] && small.Y[2] <= big.Y[2] + slack
extent_contains_extent(::Spherical, big::Extents.Extent, small::Extents.Extent; slack=1e-9) =
    all(k -> getproperty(big, k)[1] - slack <= getproperty(small, k)[1] &&
             getproperty(small, k)[2] <= getproperty(big, k)[2] + slack, (:X, :Y, :Z))
extent_contains_extent(::Spherical, big::Extents.Extent, small::SphericalCap; slack=1e-9) =
    extent_contains_extent(Spherical(), big, convert(Extents.Extent, small); slack)
extent_contains_extent(::Spherical, big::SphericalCap, small::Extents.Extent; slack=1e-9) = true  # not decidable cheaply; not checked

all_points(geom) = collect(GI.getpoint(geom))

# --- orientation -----------------------------------------------------------------
# signed turning of the ring about point c: +2π for CCW (seen from outside)
function _turning(::Spherical, ring, c)
    p = _v(c); n = length(ring); total = 0.0
    for k in 1:n
        a = _v(ring[k]); b = _v(ring[mod1(k + 1, n)])
        total += atan(_dot(p, _cross(a, b)), _dot(a, b) - _dot(a, p) * _dot(b, p))
    end
    return total
end
function _signed_area(::Planar, ring)
    n = length(ring); s = 0.0
    for k in 1:n
        a = _xy(ring[k]); b = _xy(ring[mod1(k + 1, n)])
        s += a[1] * b[2] - b[1] * a[2]
    end
    return s / 2
end
ring_is_ccw(m::Planar, ring, c) = _signed_area(m, ring) > 0
ring_mean(ring) = (s = reduce((a, b) -> a .+ b, (_v(p) for p in ring)); s ./ sqrt(_dot(s, s)))
ring_is_ccw(m::Spherical, ring, c) = _turning(m, ring, c) > π
# A point to wind about: the centroid, unless it lies on the ring (an open dual
# fan has its vertex on the ring), then the ring's own mean.
function ccw_centre(m, ring, c)
    m isa Spherical || return c
    n = length(ring)
    for k in 1:n
        a = _v(ring[k]); b = _v(ring[mod1(k + 1, n)]); p = _v(c)
        nab = _cross(a, b); nn = _norm(nab)
        nn == 0 && continue
        abs(_dot(nab, p)) <= 1e-12 * nn && _dot(_cross(a, p), nab) >= 0 && _dot(_cross(p, b), nab) >= 0 &&
            return _usp(ring_mean(ring))
    end
    return c
end
_close(::Planar, p, q; tol=1e-9) = (a = _xy(p); b = _xy(q); abs(a[1] - b[1]) + abs(a[2] - b[2]) <= tol * max(1.0, abs(a[1]) + abs(a[2])))
_close(::Spherical, p, q; tol=1e-9) = _norm(_v(p) .- _v(q)) <= tol
# is `sub` (in order, cyclically) a subsequence of `ring` under approximate equality?
function is_cyclic_subsequence(m, sub, ring)
    isempty(sub) && return true
    n = length(ring)
    for s in 1:n
        _close(m, ring[s], sub[1]) || continue
        k = 1
        for j in 0:n-1
            k <= length(sub) && _close(m, ring[mod1(s + j, n)], sub[k]) && (k += 1)
        end
        k > length(sub) && return true
    end
    return false
end

function ring_is_convex(::Planar, ring; tol=1e-12)
    n = length(ring)
    for k in 1:n
        a = _xy(ring[k]); b = _xy(ring[mod1(k + 1, n)]); c = _xy(ring[mod1(k + 2, n)])
        cr = (b[1] - a[1]) * (c[2] - b[2]) - (b[2] - a[2]) * (c[1] - b[1])
        cr < -tol * max(1.0, abs(a[1]) + abs(a[2])) && return false
    end
    return true
end
function ring_is_convex(::Spherical, ring; tol=1e-12)
    n = length(ring)
    s = (0.0, 0.0, 0.0)
    for k in 1:n
        a = _v(ring[k]); b = _v(ring[mod1(k + 1, n)]); c = _v(ring[mod1(k + 2, n)])
        _dot(_cross(a, b), c) < -tol && return false
        s = (s[1] + a[1], s[2] + a[2], s[3] + a[3])
    end
    return all(k -> _dot(_v(ring[k]), s) > 0, 1:n)    # inside an open hemisphere
end

# azimuth for winding checks
_azimuth(::Planar, c, q) = (a = _xy(c); b = _xy(q); atan(b[2] - a[2], b[1] - a[1]))
function _azimuth(::Spherical, c, q)
    n = _v(c)
    ref = abs(n[3]) < 0.9 ? (0.0, 0.0, 1.0) : (1.0, 0.0, 0.0)
    d = _dot(ref, n)
    e1 = (ref[1] - d * n[1], ref[2] - d * n[2], ref[3] - d * n[3]); e1 = e1 ./ _norm(e1)
    e2 = _cross(n, e1)
    v = _v(q)
    return atan(_dot(v, e2), _dot(v, e1))
end
# CCW order: increasing azimuth wrapping exactly once
function is_ccw_sequence(az)
    length(az) <= 2 && return true
    wraps = count(k -> az[mod1(k + 1, length(az))] < az[k], eachindex(az))
    return wraps == 1
end

# --- fallback detection ------------------------------------------------------------
# A verb is "answered by the generic" when the method selected for this grid type
# has `Any` as its first argument type.
function _is_generic(f, argtypes::Tuple)
    m = try
        which(f, Tuple{argtypes...})
    catch
        return true
    end
    sig = Base.unwrap_unionall(m.sig)
    return sig.parameters[2] === Any
end

"""
    fallback_verbs(grid) -> Vector{Symbol}

The verbs that resolve to GridInterface's generic method for this grid type.
"""
function fallback_verbs(grid)
    G = typeof(grid)
    out = Symbol[]
    for (f, args) in ((GrI.cell_centroid, (G, Int)), (GrI.cell_area, (G, Int)),
                      (GrI.cell_extent, (G, Int)), (GrI.cell_geometry, (G, Int)),
                      (GrI.cellat, (G, Any)), (GrI.treeify, (G,)),
                      (GrI.cell_contains, (G, Int, Any)), (GrI.nvertices, (G,)),
                      (GrI.cell_vertices, (G, Int)), (GrI.vertex_cells, (G, Int)))
        _is_generic(f, args) && push!(out, nameof(f))
    end
    # analytic vertices: the generic dense verbs derive from the provider's keys
    GrI.vertices_are_analytic(grid) && filter!(f -> !(f in (:nvertices, :cell_vertices, :vertex_cells)), out)
    # On a grid whose ids are its positions the generic identity is the intended
    # answer, not a fallback: list cellindex/localindex only for typed ids.
    c1 = ncells(grid) > 0 ? cellindex(grid, 1) : 1
    if !(c1 isa Integer)
        _is_generic(GrI.cellindex, (G, Int)) && push!(out, :cellindex)
        _is_generic(GrI.localindex, (G, typeof(c1))) && push!(out, :localindex)
    end
    # keyword-argument functions: check the underlying positional method
    _is_generic(GrI.neighbors, (G, Int)) && push!(out, :neighbors)
    return out
end

# geometric containing set of point p: all cells whose interior test accepts p
function containing_cells(grid, tree, m, p)
    q = _point_extent(m, p)
    out = Int[]
    STI.depth_first_search(e -> Extents.intersects(e, q), tree) do i
        GrI.cell_contains(grid, i, p) && push!(out, i)
        nothing
    end
    return sort!(unique!(out))
end

# all leaves of a tree, with their extents
function collect_leaves(tree)
    out = Int[]
    STI.depth_first_search(_ -> true, tree) do i
        push!(out, i)
        nothing
    end
    return out
end
