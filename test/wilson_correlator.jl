using LinearAlgebra
using ProgressMeter

@testset "Wilson loop correlator" begin
    @showprogress desc = "Time evolution" for L = 1:.2:2
        lat = SchwingerLattice{10,1}(L = 2, periodic=true);
        gs = groundstate(MPOHamiltonian(lat));
        gse = energy(gs);

        gsW = act(MPOWilsonLoop(lat), gs);
        _, obs = evolve(gsW, 2; observable = Dict("correlator" => (ψ, t) -> dot(gs, exp(1im*gse*t)*act(MPOWilsonLoop(lat, true), ψ))), nsteps=20)

        obs.exact = map(exp, -√(π)/2*lat.L .* (1 .- map(exp, -1im .* obs.time ./ √(π))))

        @test obs.correlator ≈ obs.exact rtol=1E-2
    end
end