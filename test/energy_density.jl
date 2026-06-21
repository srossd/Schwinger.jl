using Statistics

# energy_densities is an energy *per unit length* (electric on link n, mass on site n,
# hopping + hopping-mass on bond (n,n+1), divided by a), so Σ energy_densities · a must
# equal the total energy of the state, for every backend. (a ≠ 1 cases exercise the /a.)
@testset "Energy densities sum to total energy" begin
    for (N, a, m, θ, mp) in ((6, 1.0, 0.6, 0.0, 0.0), (6, 0.5, 0.5, 0.3, 0.4), (8, 0.7, 0.4, 0.2, 0.0))
        lat = Lattice(N; F = 1, a = a, m = m, θ2π = θ, mprime = mp)
        for backend in (EDBackend(), ITensorsBackend(), MPSKitBackend())
            gs = groundstate(Hamiltonian(lat; backend = backend))
            @test sum(energy_densities(gs)) * lat.a ≈ energy(gs) rtol = 1e-8 atol = 1e-8
        end
    end
end
