# One-off parity run (not a test): for every pair ConservativeRegridding's own
# tests regrid between (test/usecases/{constant_field,fullclenshaw,oceananigans,
# climacore}.jl and test/extensions/ringgrids_*.jl), the regridder on the native
# objects, `Regridder(dst, src)`, against the one on CR's GridInterface providers,
# `Regridder(cr_tree(grid_of(dst)), cr_tree(grid_of(src)))`: same sparse pattern
# and weights (bit-identical where the native path is deterministic, else within
# the tests' tolerances), including the tripolar fold rows (the provider regrids
# the unique cells; the native rows/columns outside them must be zero) and the
# ClimaCore SE <-> FV operators with the FV side on the provider grid.
#   julia --project=prototype/usecases/climate_grids/env prototype/usecases/conservative/cr_parity.jl [filter]
#   -> cr_parity_results.tsv next to this file
using GridInterface
import GridInterface as GrI
import ConservativeRegridding as CR
using ConservativeRegriddingTestHelpers: cr_tree, grid_of
import GeometryOps as GO, GeometryOpsCore as GOC
using SparseArrays, Printf
using Oceananigans
using Oceananigans.Grids: RightFaceFolded, RightCenterFolded
using ClimaCore: CommonSpaces, Fields, Spaces, Meshes, Topologies, Domains, ClimaComms
import Healpix
using RingGrids
using SpeedyWeather

const FILTER = length(ARGS) >= 1 ? ARGS[1] : ""
const HERE = @__DIR__
const CCExt = Base.get_extension(CR, :ConservativeRegriddingClimaCoreExt)

# --- the grids of CR's tests ------------------------------------------------------
ll(nx, ny) = LatitudeLongitudeGrid(size=(nx, ny, 1), longitude=(0, 360), latitude=(-90, 90), z=(0, 1))
rot(nx, ny) = RotatedLatitudeLongitudeGrid(size=(nx, ny, 1), longitude=(0, 360), latitude=(-90, 90), z=(0, 1), north_pole=(70, 55))
tri(nx, ny, T) = TripolarGrid(size=(nx, ny, 1), fold_topology=T)
rect(n, a, b) = RectilinearGrid(size=(n, n), x=(a, b), y=(a, b), topology=(Periodic, Periodic, Flat))
function cs_space(; h_elem=16, n_quad_points=4, use_sfc=false, radius=GO.Spherical().radius)
    context = ClimaComms.context()
    h_mesh = Meshes.EquiangularCubedSphere(Domains.SphereDomain{Float64}(radius), h_elem)
    h_topology = use_sfc ? Topologies.Topology2D(context, h_mesh, Topologies.spacefillingcurve(h_mesh)) :
        Topologies.Topology2D(context, h_mesh)
    return CommonSpaces.CubedSphereSpace(; radius=h_mesh.domain.radius, n_quad_points, h_elem, h_mesh, h_topology)
end
latlon_cc = LatitudeLongitudeGrid(size=(360, 180, 1), longitude=(0, 360), latitude=(-90, 90), z=(0, 1), radius=GO.Spherical().radius)

# (origin, name, dst thunk, src thunk, kind): kind :fv (both sides finite-volume
# cells), :se_src / :se_dst (ClimaCore spectral-element operators; the SE side
# stays native, the FV side is the provider)
PAIRS = Any[]
add!(origin, name, d, s; kind=:fv, manifold=nothing) = push!(PAIRS, (origin, name, d, s, kind, manifold))
# constant_field.jl: every spherical pair of its registry (SE<->SE dropped there too) and the planar pair
cf = [("LonLat(90x45)", () -> ll(90, 45)), ("LonLat(180x90)", () -> ll(180, 90)),
      ("TripolarF(120x60)", () -> tri(120, 60, RightFaceFolded)), ("TripolarC(120x60)", () -> tri(120, 60, RightCenterFolded)),
      ("Healpix(Nested,16)", () -> Healpix.HealpixMap{Float64,Healpix.NestedOrder}(16)),
      ("Healpix(Ring,16)", () -> Healpix.HealpixMap{Float64,Healpix.RingOrder}(16)),
      ("ClimaCore(regular)", () -> cs_space(; h_elem=16)), ("ClimaCore(Gilbert)", () -> cs_space(; h_elem=16, use_sfc=true))]
isse(k) = startswith(cf[k][1], "ClimaCore")
for i in 1:length(cf), j in i+1:length(cf)
    isse(i) && isse(j) && continue
    kind = isse(j) ? :se_src : isse(i) ? :se_dst : :fv
    add!("constant_field", cf[i][1] * " <- " * cf[j][1], cf[i][2], cf[j][2]; kind, manifold=GO.Spherical())   # the test's explicit manifold
end
add!("constant_field", "Rect(50) <- Rect(30) (planar)", () -> rect(50, 0, 1), () -> rect(30, 0, 1))
# fullclenshaw.jl
add!("fullclenshaw", "FullClenshaw(16) <- FullClenshaw(12)", () -> FullClenshawGrid(16), () -> FullClenshawGrid(12))
add!("fullclenshaw", "FullClenshaw(18) <- FullClenshaw(12)", () -> FullClenshawGrid(18), () -> FullClenshawGrid(12))
add!("fullclenshaw", "FullClenshaw(24) <- FullClenshaw(16)", () -> FullClenshawGrid(24), () -> FullClenshawGrid(16))
add!("fullclenshaw", "FullGaussian(16) <- FullClenshaw(12)", () -> FullGaussianGrid(16), () -> FullClenshawGrid(12))
add!("fullclenshaw", "FullClenshaw(12) <- FullGaussian(16)", () -> FullClenshawGrid(12), () -> FullGaussianGrid(16))
# oceananigans.jl
add!("oceananigans", "LonLat(90x45) <- LonLat(360x180)", () -> ll(90, 45), () -> ll(360, 180))
add!("oceananigans", "LonLat(360x180) <- TripolarF(360x180)", () -> ll(360, 180), () -> tri(360, 180, RightFaceFolded))
add!("oceananigans", "LonLat(360x180) <- TripolarU(360x180)", () -> ll(360, 180), () -> tri(360, 180, RightCenterFolded))
add!("oceananigans", "LonLat(90x45) <- Rotated(90x40)", () -> ll(90, 45), () -> rot(90, 40))
add!("oceananigans", "Rotated(90x40) <- LonLat(90x45)", () -> rot(90, 40), () -> ll(90, 45))
add!("oceananigans", "Rotated(90x40) <- TripolarF(120x60)", () -> rot(90, 40), () -> tri(120, 60, RightFaceFolded))
add!("oceananigans", "Rotated(90x40) <- TripolarU(120x60)", () -> rot(90, 40), () -> tri(120, 60, RightCenterFolded))
add!("oceananigans", "Rect(100,0:2) <- Rect(200,0:1) (planar)", () -> rect(100, 0, 2), () -> rect(200, 0, 1))
# climacore.jl
add!("climacore", "LonLat(360x180) <- SE regular h16", () -> latlon_cc, () -> cs_space(; h_elem=16); kind=:se_src)
add!("climacore", "LonLat(360x180) <- SE Gilbert h16", () -> latlon_cc, () -> cs_space(; h_elem=16, use_sfc=true); kind=:se_src)
add!("climacore", "LonLat(360x180) <- SE regular h8", () -> latlon_cc, () -> cs_space(; h_elem=8); kind=:se_src)
add!("climacore", "SE regular h16 <- LonLat(360x180)", () -> cs_space(; h_elem=16), () -> latlon_cc; kind=:se_dst)
add!("climacore", "SE Gilbert h16 <- LonLat(360x180)", () -> cs_space(; h_elem=16, use_sfc=true), () -> latlon_cc; kind=:se_dst)
const R_TRI = TripolarGrid(size=(360, 180, 1)).radius   # Oceananigans' default Earth radius
add!("climacore", "SE h32 <- TripolarU(360x180)", () -> CommonSpaces.CubedSphereSpace(; radius=R_TRI, n_quad_points=4, h_elem=32), () -> tri(360, 180, RightCenterFolded); kind=:se_dst)
add!("climacore", "SE h32 <- TripolarF(360x180)", () -> CommonSpaces.CubedSphereSpace(; radius=R_TRI, n_quad_points=4, h_elem=32), () -> tri(360, 180, RightFaceFolded); kind=:se_dst)
# ringgrids extensions
add!("ringgrids_healpix", "HEALPix(16) <- HEALPix(8)", () -> HEALPixGrid(16), () -> HEALPixGrid(8))
add!("ringgrids_healpix", "HEALPix(8) <- HEALPix(16)", () -> HEALPixGrid(8), () -> HEALPixGrid(16))
add!("ringgrids_healpix", "HEALPix(10) <- HEALPix(6)", () -> HEALPixGrid(10), () -> HEALPixGrid(6))
add!("ringgrids_healpix", "FullClenshaw(16) <- HEALPix(8)", () -> FullClenshawGrid(16), () -> HEALPixGrid(8))
add!("ringgrids_octahealpix", "OctaHEALPix(16) <- OctaHEALPix(8)", () -> OctaHEALPixGrid(16), () -> OctaHEALPixGrid(8))
add!("ringgrids_octahealpix", "OctaHEALPix(8) <- OctaHEALPix(16)", () -> OctaHEALPixGrid(8), () -> OctaHEALPixGrid(16))
add!("ringgrids_octahealpix", "FullClenshaw(16) <- OctaHEALPix(8)", () -> FullClenshawGrid(16), () -> OctaHEALPixGrid(8))
add!("ringgrids_octahealpix", "OctaHEALPix(10) <- OctaHEALPix(6)", () -> OctaHEALPixGrid(10), () -> OctaHEALPixGrid(6))
add!("ringgrids_reduced", "FullGaussian(9) <- OctahedralGaussian(12)", () -> FullGaussianGrid(9), () -> OctahedralGaussianGrid(12))
add!("ringgrids_reduced", "OctaminimalGaussian(10) <- OctahedralClenshaw(8)", () -> OctaminimalGaussianGrid(10), () -> OctahedralClenshawGrid(8))
add!("ringgrids_reduced", "OctaHEALPix(8) <- OctahedralGaussian(10)", () -> OctaHEALPixGrid(8), () -> OctahedralGaussianGrid(10))
add!("ringgrids_reduced", "OctahedralGaussian(8) <- HEALPix(8)", () -> OctahedralGaussianGrid(8), () -> HEALPixGrid(8))

# --- the comparison -----------------------------------------------------------------
_pp(g) = g isa GrI.Subset ? collect(GrI.parentpositions(g)) : nothing
function compare(Rn, Rp, rows, cols)
    A = Rn.intersections; da = Rn.dst_areas; sa = Rn.src_areas
    zr = rows === nothing || nnz(A[setdiff(1:size(A, 1), rows), :]) == 0
    zc = cols === nothing || nnz(A[:, setdiff(1:size(A, 2), cols)]) == 0
    rows === nothing || (A = A[rows, :]; da = da[rows])
    cols === nothing || (A = A[:, cols]; sa = sa[cols])
    B = Rp.intersections
    size(A) == size(B) || return (; ident=false, pattern=false, rel=NaN, rel_area=NaN, outside_zero=(zr, zc), size=(size(A), size(B)))
    ident = A == B && Rp.dst_areas == da && Rp.src_areas == sa
    pattern = findnz(A)[1:2] == findnz(B)[1:2]
    rel = maximum(abs, nonzeros(A - B); init=0.0) / maximum(abs, nonzeros(A); init=1.0)
    ra = max(maximum(abs.(Rp.dst_areas .- da) ./ max.(abs.(da), eps()); init=0.0),
             maximum(abs.(Rp.src_areas .- sa) ./ max.(abs.(sa), eps()); init=0.0))
    return (; ident, pattern, rel, rel_area=ra, outside_zero=(zr, zc), size=(size(A), size(B)))
end

out = open(joinpath(HERE, isempty(FILTER) ? "cr_parity_results.tsv" : "cr_parity_results_$(FILTER).tsv"), "w")
println(out, join(["origin", "pair", "kind", "ndst", "nsrc", "verdict", "pattern", "rel_weights", "rel_areas", "outside_zero", "nnz", "t_native_s", "t_provider_s", "note"], '\t'))
ntot = 0; nident = 0; nclose = 0; nbad = 0
for (origin, name, dmk, smk, kind, manifold) in PAIRS
    (isempty(FILTER) || occursin(FILTER, origin * " " * name)) || continue
    dst = dmk(); src = smk()
    D = kind === :se_dst ? dst : grid_of(dst)
    S = kind === :se_src ? src : grid_of(src)
    dT = kind === :se_dst ? dst : cr_tree(D)
    sT = kind === :se_src ? src : cr_tree(S)
    note = String[]
    local Rn, Rp, tn, tp
    try
        tn = @elapsed Rn = manifold === nothing ? CR.Regridder(dst, src; threaded=false) : CR.Regridder(manifold, dst, src; threaded=false)
    catch e
        println(rpad(name, 48), " native path errors: ", first(sprint(showerror, e), 200))
        println(out, join([origin, name, kind, "", "", "native errors", "", "", "", "", "", "", "", first(replace(sprint(showerror, e), '\n' => ' '), 160)], '\t')); flush(out)
        global ntot += 1; global nbad += 1
        continue
    end
    try
        tp = @elapsed Rp = manifold === nothing ? CR.Regridder(dT, sT; threaded=false) : CR.Regridder(manifold, dT, sT; threaded=false)
    catch e
        println(rpad(name, 48), " provider path errors: ", first(sprint(showerror, e), 300))
        println(out, join([origin, name, kind, "", "", "provider errors", "", "", "", "", "", @sprintf("%.3f", tn), "", first(replace(sprint(showerror, e), '\n' => ' '), 160)], '\t')); flush(out)
        global ntot += 1; global nbad += 1
        continue
    end
    rows = kind === :se_dst ? nothing : _pp(D)
    cols = kind === :se_src ? nothing : _pp(S)
    c = compare(Rn, Rp, rows, cols)
    verdict = c.ident ? "identical" : (c.pattern && c.rel < 1e-12 && c.rel_area < 1e-12) ? "within 1e-12" :
        (c.pattern && c.rel < 1e-8 && c.rel_area < 1e-8) ? "within 1e-8" : "DIFFERENT"
    all(c.outside_zero) || push!(note, "native rows/cols outside the unique cells not all zero")
    kind === :fv || push!(note, "SE side native, FV side provider")
    global ntot += 1
    c.ident ? (global nident += 1) : verdict == "DIFFERENT" ? (global nbad += 1) : (global nclose += 1)
    row = [origin, name, string(kind), c.size[1][1], c.size[1][2], verdict, c.pattern, @sprintf("%.2e", c.rel), @sprintf("%.2e", c.rel_area),
        c.outside_zero, nnz(Rp.intersections), @sprintf("%.3f", tn), @sprintf("%.3f", tp), join(note, "; ")]
    println(rpad(name, 48), " ", rpad(verdict, 13), " pattern=", c.pattern, " rel=", @sprintf("%.1e", c.rel), " areas=", @sprintf("%.1e", c.rel_area),
        " t=", @sprintf("%.2f", tn), "s/", @sprintf("%.2f", tp), "s", isempty(note) ? "" : "  [" * join(note, "; ") * "]"); flush(stdout)
    println(out, join(row, '\t')); flush(out)
end
close(out)
println("\n", ntot, " pairs: ", nident, " identical, ", nclose, " within tolerance, ", nbad, " different or erroring")
