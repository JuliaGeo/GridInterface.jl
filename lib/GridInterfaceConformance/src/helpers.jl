# --- seeded sampling -------------------------------------------------------------
"""
    DEFAULT_SEED

The seed every `check_*` / `test_*` call samples with unless given an `rng`.
The report header prints the seed a run used, so a failure is reproduced by
running the same call again; passing `rng = MersenneTwister(s)` explores a
different sample (an RNG advances as it is drawn from, so reproducing a run
means a fresh generator on the same seed, not the one it consumed).
"""
const DEFAULT_SEED = 20260813
_default_rng(seed::Integer=DEFAULT_SEED) = Random.MersenneTwister(seed)

"""
    sample_positions(rng, n, k) -> Vector{Int}

Up to `k` distinct positions of `1:n`, sorted, always holding `1` and `n` (the
first and last cells are where edge and pole arithmetic breaks); every
position when `n <= k`. Drawn with replacement and deduplicated, so a level
grid with more cells than fit in memory never allocates a permutation.
"""
function sample_positions(rng::Random.AbstractRNG, n::Integer, k::Integer)
    n <= 0 && return Int[]
    n <= k && return collect(1:Int(n))
    k <= 2 && return [1, Int(n)]
    return sort!(unique!(vcat([1, Int(n)], rand(rng, 1:Int(n), Int(k) - 2))))
end

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
"""
    has_nonfallback_method(f, args...) -> Bool

Whether dispatch for `f(args...)` reaches a method written for the grid or
system among `args` rather than an interface-wide generic. "A fallback
answering in place of a declared analytic trait" is a conformance failure, and
this is how the suite tells the two apart.

# Method selection

The subjects are the arguments that claim the contract (`isgrid` or
`issystem`). The matched method counts as specialized when its parameter at a
subject's position is **strictly narrower** than the interface-wide subject
type: GridInterface's generics take the grid as `Any`, so anything a provider
wrote (its own type, a family supertype such as `AbstractStructuredGrid` or
DGG's `AbstractGrid`) is narrower. Module ownership is irrelevant: a module may
define both generic fallbacks and specialized methods.

A transparent wrapper ([`GridInterface.AbstractWrappedGrid`](@ref)) forwards
every verb it does not hold to its parent, so a method whose subject parameter
is the wrapper supertype is a forward, not an answer: the question is asked of
the parent instead. A wrapper's own override (`cellat` on a `GridAndTree`) is
narrower than the supertype and counts as specialized.

This function reports only whether a specialized method exists; the laws
separately validate its result.
"""
function has_nonfallback_method(f, args...)
    applicable(f, args...) || return false
    m = try
        which(f, Base.typesof(args...))
    catch err
        # Ambiguous dispatch does not identify a specialized implementation.
        (err isa ArgumentError || err isa MethodError) || rethrow()
        return false
    end
    positions = _subject_positions(args)
    isempty(positions) && return false
    has_specific_subject(m, positions) && return true
    for k in positions
        a = args[k]
        a isa GrI.AbstractWrappedGrid && return has_nonfallback_method(f, Base.setindex(args, parent(a), k)...)
    end
    return false
end

_subject_positions(args) = [k for k in eachindex(args) if GrI.isgrid(args[k]) || GrI.issystem(args[k])]

"""
    has_specific_subject(m::Method, positions) -> Bool

Whether `m` takes, at one of the argument `positions` (1-based, after the
function), a type **strictly narrower** than the interface-wide subject types
(`Any`, and `AbstractWrappedGrid` for the forwards): the specificity test behind
[`has_nonfallback_method`](@ref).

Each parameter of the signature is rewrapped in the method's own `where`
clauses before the comparison, so a parametric subject type (`Subset{G,V,R}`,
`GridAndTree{G,T}`) is compared as the closed type `Subset{G,V,R} where {G,V,R}`
rather than as an open body whose free type variables make subtyping
undecidable.

"Strictly narrower" is spelled `!(base <: T)` rather than `T !== base` on
purpose: a method written `f(g::G) where {G}` has a parameter that is *equal*
to `Any` without being `===` to it, and that method is the interface-wide
generic, not an implementation.
"""
function has_specific_subject(m::Method, positions)
    sig = m.sig
    body = Base.unwrap_unionall(sig)
    body isa DataType || return false
    params = body.parameters
    # Parameter 1 is `typeof(f)`; subject k is parameter k + 1.
    for k in positions
        k + 1 <= length(params) || continue
        T = Base.rewrap_unionall(params[k+1], sig)
        T isa Type || continue
        _specific_subject(T) && return true
    end
    return false
end
_specific_subject(T) = !(Any <: T) && !(GrI.AbstractWrappedGrid <: T)

"""
    fallback_verbs(grid) -> Vector{Symbol}

The verbs that resolve to GridInterface's generic method for this grid
(through any transparent wrappers), by [`has_nonfallback_method`](@ref).
"""
function fallback_verbs(grid)
    out = Symbol[]
    n = ncells(grid)
    n > 0 || return out
    # a grid whose centroid throws (a failed type contract) lists the point
    # verbs as fallbacks rather than aborting the report
    p = try cell_centroid(grid, 1) catch err; err isa InterruptException && rethrow(); nothing end
    for (f, args) in ((GrI.cell_centroid, (1,)), (GrI.cell_area, (1,)), (GrI.cell_extent, (1,)),
                      (GrI.cell_geometry, (1,)), (GrI.cellat, (p,)), (GrI.treeify, (3,)),
                      (GrI.cell_contains, (1, p)), (GrI.nvertices, ()),
                      (GrI.cell_vertices, (1,)), (GrI.vertex_cells, (1,)))
        has_nonfallback_method(f, grid, args...) || push!(out, nameof(f))
    end
    # analytic vertices: the generic dense verbs derive from the provider's keys
    GrI.vertices_are_analytic(grid) && filter!(f -> !(f in (:nvertices, :cell_vertices, :vertex_cells)), out)
    # On a grid whose ids are its positions the generic identity is the intended
    # answer, not a fallback: list cellindex/localindex only for typed ids.
    c1 = cellindex(grid, 1)
    if !(c1 isa Integer)
        has_nonfallback_method(GrI.cellindex, grid, 1) || push!(out, :cellindex)
        has_nonfallback_method(GrI.localindex, grid, c1) || push!(out, :localindex)
    end
    # keyword-argument functions: check the underlying positional method
    has_nonfallback_method(GrI.neighbors, grid, 1) || push!(out, :neighbors)
    return out
end

# A hand-written pass-through (not an AbstractWrappedGrid) over a grid with the
# same positions: `parent` is a grid of the same size. Subsets, duals and
# multi-level sets remap positions and are not wrappers.
function handwritten_wrapper_parent(grid)
    (grid isa GrI.AbstractWrappedGrid || grid isa GrI.Subset || grid isa GrI.DualGrid || grid isa GrI.MultiLevelSet) && return nothing
    hasmethod(Base.parent, Tuple{typeof(grid)}) || return nothing
    p = parent(grid)
    return GrI.isgrid(p) && ncells(p) == ncells(grid) ? p : nothing
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

# the leaf indices under `node`, in traversal order (STI's default `subtree_indices`)
function leaf_walk(node)
    out = Int[]
    STI.depth_first_search(_ -> true, node) do i
        push!(out, i)
        nothing
    end
    return out
end

# the first leaf of a tree (descending first children), or nothing
function first_leaf(node)
    STI.isleaf(node) && return node
    for c in STI.getchild(node)
        return first_leaf(c)
    end
    return nothing
end

# the root plus the nodes on `k` random descents (a few internal nodes and leaves)
function sample_nodes(rng, tree, k)
    out = Any[tree]
    for _ in 1:k
        node = tree
        while !STI.isleaf(node)
            cs = collect(STI.getchild(node))
            isempty(cs) && break
            node = cs[rand(rng, 1:length(cs))]
            push!(out, node)
        end
    end
    return out
end
