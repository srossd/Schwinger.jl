# §3.4: opt-in `pad` for boundary-limited detector sweeps returns a site-aligned length-N vector
# with NaN at the sites the operator is undefined on (energycurrents: sites 1,N; momentumdensities:
# sites N-1,N), so whole-lattice sweeps/plots don't need manual range guards.
using Schwinger, Test

@testset "detector pad option" begin
    N = 8
    lat = Lattice(N; F = 1, q = 2, a = 0.25, m = 0.5)
    for be in (EDBackend(), ITensorsBackend(), MPSKitBackend())
        gs = groundstate(Hamiltonian(lat, be); (be isa EDBackend ? (;) : (; energy_tol = 1e-10))...)

        ec  = energycurrents(gs)
        ecp = energycurrents(gs; pad = true)
        @test length(ec)  == N - 2
        @test length(ecp) == N
        @test isnan(ecp[1]) && isnan(ecp[N])
        @test ecp[2:N-1] == ec                     # interior identical, NaN only at the ends

        md  = momentumdensities(gs)
        mdp = momentumdensities(gs; pad = true)
        @test length(md)  == N - 2
        @test length(mdp) == N
        @test isnan(mdp[N-1]) && isnan(mdp[N])
        @test mdp[1:N-2] == md
    end
end
