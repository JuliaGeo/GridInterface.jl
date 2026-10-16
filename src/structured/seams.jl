# ===========================================================================
# Seams of a multi-face grid found from its vertices: the table of which
# face and side lies across every side of every face, built once from the
# faces' corner vertices, and the index arithmetic of crossing it. A
# provider whose faces share whole sides (a cubed sphere, the HEALPix base
# faces) stores the table and forwards `facecross` to `seamcross`.
# ===========================================================================

"""
    FaceSeams(grid; tol = 1e-12)

The seam table of a structured grid whose faces share whole sides: for every
face `f` and side `s`, the face and side `(f2, s2)` across it, or `(0, 0)` at
a boundary. Built once from the two end vertices of every side (the faces'
corner vertices through [`facevertex`](@ref)): side `s` of `f` lies against
side `s2` of `f2` when the ends coincide to `tol` (relative) and the sides
hold the same number of cells. Faces must be consistently oriented, so the
two faces run a shared side in opposite directions ([`facecross`](@ref)'s
contract); a side matched in the same direction throws. A side whose two
ends coincide (a collapsed row) is left a boundary: declare its fan with
[`facefan`](@ref). O(nfaces²) vertex reads.

    seamcross(seams, grid, f, i, j, s) -> (f2, i2, j2, s2) | nothing

The crossing of side `s` of cell `(i, j)` through the table: cell at
position `p` along side `s` lands on the cell at position `n + 1 - p` along
side `s2`, the cells of a side counted counter-clockwise about the face.
"""
struct FaceSeams
    across::Matrix{Tuple{Int,Int}}        # 4 × nfaces: (f2, s2), or (0, 0)
end

@inline _facecornervertex(g, f, ni, nj, c) =
    c == 1 ? facevertex(g, f, 1, 1) : c == 2 ? facevertex(g, f, ni + 1, 1) :
    c == 3 ? facevertex(g, f, ni + 1, nj + 1) : facevertex(g, f, 1, nj + 1)
@inline _sidelen(ni, nj, s) = isodd(s) ? ni : nj
# position p (1..len) of cell (i, j) along side s, counter-clockwise, and back
@inline _sidepos(ni, nj, i, j, s) = s == 1 ? i : s == 2 ? j : s == 3 ? ni + 1 - i : nj + 1 - j
@inline _sidecell(ni, nj, s, p) = s == 1 ? (p, 1) : s == 2 ? (ni, p) : s == 3 ? (ni + 1 - p, nj) : (1, nj + 1 - p)

@inline function _near_tol(a, b, tol)
    x = GrI._coords(a); y = GrI._coords(b)
    for k in eachindex(x)
        abs(x[k] - y[k]) <= tol * (1 + abs(x[k])) || return false
    end
    return true
end

function FaceSeams(g::AbstractStructuredGrid; tol::Real=1e-12)
    nf = nfaces(g)
    near(a, b) = _near_tol(a, b, tol)
    sizes = [facesize(g, f) for f in 1:nf]
    corners = [ntuple(c -> _facecornervertex(g, f, sizes[f][1], sizes[f][2], c), 4) for f in 1:nf]
    across = fill((0, 0), 4, nf)
    for f in 1:nf, s in 1:4
        p, q = corners[f][s], corners[f][mod1(s + 1, 4)]
        near(p, q) && continue                       # collapsed side: a fan, not a seam
        len = _sidelen(sizes[f]..., s)
        found = (0, 0)
        for f2 in 1:nf, s2 in 1:4
            (f2 == f && s2 == s) && continue
            p2, q2 = corners[f2][s2], corners[f2][mod1(s2 + 1, 4)]
            if near(p, q2) && near(q, p2)
                _sidelen(sizes[f2]..., s2) == len ||
                    throw(ArgumentError("$(typeof(g)): side $s of face $f and side $s2 of face $f2 share their ends but not their cell count"))
                found = (f2, s2)
                break
            elseif near(p, p2) && near(q, q2)
                throw(ArgumentError("$(typeof(g)): side $s of face $f runs the same way as side $s2 of face $f2: faces are not consistently oriented"))
            end
        end
        across[s, f] = found
    end
    # a seam crosses back
    for f in 1:nf, s in 1:4
        f2, s2 = across[s, f]
        f2 == 0 && continue
        across[s2, f2] == (f, s) ||
            throw(ArgumentError("$(typeof(g)): side $s of face $f lands on side $s2 of face $f2, which does not land back"))
    end
    return FaceSeams(across)
end

@inline function seamcross(t::FaceSeams, g, f, i, j, s)
    f2, s2 = @inbounds t.across[s, f]
    f2 == 0 && return nothing
    ni, nj = facesize(g, f)
    n2i, n2j = facesize(g, f2)
    p = _sidepos(ni, nj, i, j, s)
    i2, j2 = _sidecell(n2i, n2j, s2, _sidelen(n2i, n2j, s2) + 1 - p)
    return (f2, i2, j2, s2)
end
Base.show(io::IO, t::FaceSeams) = print(io, "FaceSeams(", size(t.across, 2), " faces, ",
    count(x -> x[1] != 0, t.across), " seam sides)")
