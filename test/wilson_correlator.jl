using LinearAlgebra
using ProgressMeter

@testset "Wilson loop correlator" begin
    @showprogress desc = "Wilson loop correlator" for L = 1:.2:2
        lat = Lattice(10; F = 1, L = 2, periodic=true);
        gs = groundstate(Hamiltonian(lat, ITensorsBackend()));
        gse = energy(gs);

        gsW = act(WilsonLoop(lat; backend = ITensorsBackend()), gs);
        _, obs = evolve(gsW, 2; observable = Dict("correlator" => (ψ, t) -> dot(gs, exp(1im*gse*t)*act(WilsonLoop(lat, true; backend = ITensorsBackend()), ψ))), nsteps=20)

        obs.exact = map(exp, -√(π)/2*lat.L .* (1 .- map(exp, -1im .* obs.time ./ √(π))))

        @test obs.correlator ≈ obs.exact rtol=1E-2
    end
end