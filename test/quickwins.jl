# Smoke + light-correctness tests for the "quick win" additions:
#   2.3 generic chargecurrents/energycurrents/momentumdensities (ED, ITensors)
#   2.7 dispersion / groupvelocity
#   2.9 rehost / quench
#   2.8 normalize / normalize!
#   2.10 standard_densities observable menu
using Schwinger, LinearAlgebra, Test

@testset "quick wins" begin
    N = 8
    lat = Lattice(N; F = 1, q = 2, a = 0.25, m = 0.5)     # mprime = 0

    @testset "2.3 generic current profiles (ED, ITensors)" begin
        for be in (EDBackend(), ITensorsBackend())
            H  = Hamiltonian(lat, be)
            gs = groundstate(H; (be isa EDBackend ? (;) : (; energy_tol = 1e-10))...)

            jc = chargecurrents(gs)
            ec = energycurrents(gs)
            md = momentumdensities(gs)

            @test length(jc) == N - 1
            @test length(ec) == N - 2          # sites 2..N-1
            @test length(md) == N - 2          # sites 1..N-2
            @test eltype(jc) <: Real

            # wraps the validated per-bond/per-site operators — confirm identical to a manual loop
            kw = be isa EDBackend ?
                 (; L_max = H.L_max, universe = H.universe, charge = 0) :
                 (; L_max = H.L_max, universe = H.universe)
            jc_ref = [real(expectation(ChargeCurrent(lat, b; backend = be, kw...), gs)) for b in 1:N-1]
            @test jc ≈ jc_ref atol=1e-12

            # parity-symmetric ground state ⇒ currents & momentum vanish
            @test maximum(abs, jc) < 1e-6
            @test maximum(abs, ec) < 1e-6
            @test abs(totalmomentum(gs)) < 1e-6
            @test totalmomentum(gs) ≈ sum(md) atol=1e-12
        end
    end

    @testset "2.9 rehost / quench" begin
        H0 = Hamiltonian(lat, EDBackend())
        gs = groundstate(H0)
        lat1 = Lattice(N; F = 1, q = 2, a = 0.25, m = 2.0)   # heavier mass
        H1 = Hamiltonian(lat1, EDBackend())

        q = rehost(gs, H1)
        @test q isa EDState
        @test q.hamiltonian === H1
        @test q.coeffs === gs.coeffs                 # data shared, not copied
        @test quench(gs, H1).hamiltonian === H1      # alias
        # energy under the new H differs from the old ground-state energy
        @test !isapprox(energy(q), energy(gs); atol = 1e-6)
        # <ψ|H1|ψ> ≥ groundstate energy of H1 (variational)
        @test energy(q) > energy(groundstate(H1)) - 1e-8

        # ED sector-mismatch guard
        latbig = Lattice(N + 2; F = 1, q = 2, a = 0.25, m = 0.5)
        @test_throws ArgumentError rehost(gs, Hamiltonian(latbig, EDBackend()))
    end

    @testset "2.8 normalize / normalize!" begin
        H  = Hamiltonian(lat, EDBackend())
        gs = groundstate(H)
        # act returns an unnormalized O|ψ>
        Oψ = act(Mass(lat; backend = EDBackend(), bare = false), gs)
        @test !isapprox(norm(Oψ.coeffs), 1.0; atol = 1e-6)
        n1 = normalize(Oψ)
        @test norm(n1.coeffs) ≈ 1.0 atol=1e-12
        @test norm(Oψ.coeffs) != 1.0                 # normalize made a fresh state
        normalize!(Oψ)
        @test norm(Oψ.coeffs) ≈ 1.0 atol=1e-12
    end

    @testset "2.10 standard_densities" begin
        d = standard_densities([:charge, :current, :energy])
        @test d isa Dict{String,Function}
        @test Set(keys(d)) == Set(["charge", "current", "energy"])
        @test_throws ArgumentError standard_densities([:not_a_density])

        H  = Hamiltonian(lat, EDBackend())
        gs = groundstate(H)
        # each callback produces a vector on the state
        @test d["charge"](gs, 0.0) isa AbstractVector
        @test length(d["current"](gs, 0.0)) == N - 1
        # plugs into evolve
        _, obs = evolve(gs, 0.05; nsteps = 2, observable = standard_densities([:charge, :current]))
        @test length(obs.charge) == 2 && length(obs.current) == 2   # one row per step
    end

    @testset "2.8 savestate / loadstate" begin
        tmp = tempname() * ".jld2"
        # ED round-trip: coeffs, defects, charge sector
        H  = Hamiltonian(lat, EDBackend())
        gs = groundstate(H)
        savestate(tmp, gs)
        gs2 = loadstate(tmp, H)
        @test gs2 isa EDState
        @test gs2.coeffs ≈ gs.coeffs atol=1e-14
        @test gs2.net_charge == gs.net_charge
        @test energy(gs2) ≈ energy(gs) atol=1e-12

        # backend-mismatch guard
        @test_throws ArgumentError loadstate(tmp, Hamiltonian(lat, ITensorsBackend()))

        # ITensors round-trip (compare observables, not the psi object)
        HI  = Hamiltonian(lat, ITensorsBackend())
        gsI = groundstate(HI; energy_tol = 1e-10)
        tmpI = tempname() * ".jld2"
        savestate(tmpI, gsI)
        gsI2 = loadstate(tmpI, HI)
        @test gsI2 isa ITensorState
        @test real.(vec(charges(gsI2))) ≈ real.(vec(charges(gsI))) atol=1e-8
        @test energy(gsI2) ≈ energy(gsI) atol=1e-8

        rm(tmp; force = true); rm(tmpI; force = true)
    end

    @testset "2.7 dispersion / groupvelocity (infinite MPSKit)" begin
        latI = Lattice(Inf; F = 1, a = 0.5, m = 1.0)
        HI   = Hamiltonian(latI, MPSKitBackend())
        ps   = [-0.4, 0.0, 0.4]
        Es   = dispersion(HI, ps; bonddim = 24, energy_tol = 1e-8)
        @test length(Es) == length(ps)
        @test all(>(0), Es)                          # gapped meson band
        @test Es[2] ≤ Es[1] + 1e-6 && Es[2] ≤ Es[3] + 1e-6   # minimum at k = 0
        vg0 = groupvelocity(HI, 0.0; dp = 0.2, bonddim = 24, energy_tol = 1e-8)
        @test abs(vg0) < 5e-2                         # zero group velocity at the band minimum
    end
end
