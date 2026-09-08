using Test
using Schwinger

# Passing a precomputed ground state to `loweststates` must reproduce exactly the states obtained
# by re-solving from scratch (the cached path only skips the GS solve; the excitation solve on top
# is unchanged). This lets a cached GS be reused to find excitations at several momenta cheaply.

@testset "MPSKit infinite: cached GS reproduces excitations" begin
    latinf = Lattice(Inf; F = 1, a = 0.2, m = 1.0)
    Hinf   = Hamiltonian(latinf; backend = :MPSKit)

    # Baseline: solve GS + excitation from scratch.
    scratch = loweststates(Hinf, 2; bonddim = 24, energy_tol = 1e-9, momentum = 0.5)

    # Cache the GS once, then reuse it for the excitation — no second VUMPS solve.
    gs    = groundstate(Hinf; bonddim = 24, energy_tol = 1e-9)
    cached = loweststates(Hinf, 2; groundstate = gs, momentum = 0.5)

    # The cached state is returned verbatim as state 1.
    @test cached[1] === gs
    @test isapprox(energy(scratch[1]), energy(cached[1]); rtol = 1e-8)
    # Excitation built on the reused GS matches the from-scratch excitation.
    @test isapprox(energy(scratch[2]), energy(cached[2]); rtol = 1e-6)

    # The point of caching: reuse the same GS at a different momentum without re-solving.
    other = loweststates(Hinf, 2; groundstate = gs, momentum = 0.25)
    @test isfinite(energy(other[2]))
    # A distinct momentum gives a distinct excitation energy (dispersion is nontrivial).
    @test !isapprox(energy(other[2]), energy(cached[2]); rtol = 1e-3)
end

@testset "ITensor finite: cached GS reproduces excited states" begin
    lat = Lattice(6; F = 1, m = 1.0)
    H   = Hamiltonian(lat; backend = :ITensors)

    scratch = loweststates(H, 2; energy_tol = 1e-10)
    gs      = groundstate(H; energy_tol = 1e-10)
    cached  = loweststates(H, 2; groundstate = gs, energy_tol = 1e-10)

    @test cached[1] === gs
    @test isapprox(energy(scratch[1]), energy(cached[1]); rtol = 1e-6)
    # First excited state, found by DMRG penalizing the (cached) GS, matches the scratch run.
    @test isapprox(energy(scratch[2]), energy(cached[2]); rtol = 1e-5)
end

@testset "cached GS validation" begin
    Hinf = Hamiltonian(Lattice(Inf; F = 1, m = 1.0); backend = :MPSKit)
    Hfin = Hamiltonian(Lattice(6; F = 1, m = 1.0); backend = :MPSKit)

    # Wrong type.
    @test_throws ArgumentError loweststates(Hinf, 2; groundstate = 42)

    # Finite GS handed to an infinite solve (and vice versa) is rejected.
    gs_fin = groundstate(Hfin; bonddim = 16)
    @test_throws ArgumentError loweststates(Hinf, 2; groundstate = gs_fin)

    gs_inf = groundstate(Hinf; bonddim = 16, energy_tol = 1e-8)
    @test_throws ArgumentError loweststates(Hfin, 2; groundstate = gs_inf)
end
