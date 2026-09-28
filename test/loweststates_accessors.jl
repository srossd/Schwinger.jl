# §3.1: `vacuumof` / `quasiparticle` normalize the ordinary-vs-soliton and single-vs-momentum-list
# nesting of a `loweststates` result, so callers stop indexing `res[2][j][1]` by hand.
using Schwinger, Test

@testset "loweststates accessors" begin
    latI = Lattice(Inf; F = 1, a = 0.5, m = 1.0)
    HI   = Hamiltonian(latI, MPSKitBackend())

    # ordinary, single momentum: res[1] vacuum, res[2] a QP state
    res1 = loweststates(HI, 2; bonddim = 16, energy_tol = 1e-8)
    @test vacuumof(res1) === res1[1]
    @test quasiparticle(res1, 2) === res1[2]
    @test_throws ArgumentError vacuumof(res1; which = :second)          # single vacuum
    @test_throws ArgumentError quasiparticle(res1, 2; kind = :antiparticle)  # no partner
    @test_throws ArgumentError quasiparticle(res1, 1)                   # level 1 is the vacuum

    # ordinary, momentum list: res[2] is a vector of QPs (reuse the vacuum just solved)
    res2 = loweststates(HI, 2; momentum = [-0.3, 0.3], groundstate = res1[1], bonddim = 16)
    @test quasiparticle(res2, 2; momentum = 1) === res2[2][1]
    @test quasiparticle(res2, 2; momentum = 2) === res2[2][2]
    @test_throws ArgumentError quasiparticle(res2, 2; momentum = 3)

    # solitons at θ = π: res[1] = (v1, v2), res[2] = (soliton, antisoliton)
    latP = Lattice(Inf; F = 1, a = 0.5, m = 0.5, θ2π = 0.5)
    HP   = Hamiltonian(latP, MPSKitBackend())
    resS = loweststates(HP, 2; solitons = true, bonddim = 16, energy_tol = 1e-8)
    @test resS[1] isa Tuple && resS[2] isa Tuple
    @test vacuumof(resS)               === resS[1][1]
    @test vacuumof(resS; which = :second) === resS[1][2]
    @test quasiparticle(resS, 2)                    === resS[2][1]   # soliton
    @test quasiparticle(resS, 2; kind = :antiparticle) === resS[2][2]  # antisoliton
    @test quasiparticle(resS, 2; kind = :soliton)    === resS[2][1]   # alias
end
