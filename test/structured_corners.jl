module CornersRegression
using Test
import GridInterface as GrI
import GridInterface.Structured as S
using GeometryOps.UnitSpherical: UnitSphericalPoint

@testset "Corners preserve side identities and winding through degeneracies" begin
    a, b, c, d = (0.,0.), (2.,0.), (2.,1.), (0.,1.)
    cases = (
        ((a,b,c,d), (1,2,3,4), [a,b,c,d]),
        ((a,a,c,d), (2,3,4,0), [a,c,d]),
        ((a,b,c,c), (1,2,4,0), [a,b,c]),
        ((a,b,c,a), (1,2,3,0), [a,b,c]),
        ((a,a,a,d), (3,4,0,0), [a,d]),
        ((a,a,c,c), (2,4,0,0), [a,c]),
        ((a,a,a,a), (0,0,0,0), Tuple{Float64,Float64}[]),
        # Repeated nonadjacent vertices do not make adjacent sides degenerate.
        ((a,b,a,d), (1,2,3,4), [a,b,a,d]),
        # Nearness is not transitive: only the closing side survives here.
        (((0.,0.),(.75e-12,0.),(1.5e-12,0.),(2.25e-12,0.)), (4,0,0,0), [(2.25e-12,0.)]),
    )
    for (points, sides, expected) in cases
        r = S.Corners(points)
        @test r.pts === points
        @test r.sides == sides
        @test length(r) == length(expected)
        @test collect(r) == expected
        @test_throws BoundsError r[0]
        @test_throws BoundsError r[length(r)+1]
    end
    # A collapsed south-pole side keeps the two meridian sides in their order.
    pole = UnitSphericalPoint(0.,0.,-1.)
    east = UnitSphericalPoint(1.,0.,0.)
    north = UnitSphericalPoint(0.,1.,0.)
    r = S.Corners((pole,pole,east,north))
    @test r.sides == (2,3,4,0)
    @test collect(r) == [pole,east,north]
end
end
