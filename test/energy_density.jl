using Statistics

# The per-site energies (electric on link n, mass on site n, hopping + hopping-mass on
# bond (n,n+1)) must sum to the total energy of the state, for every backend.
@testset "Energy densities sum to total energy" begin
    for (N, m, θ, mp) in ((6, 0.6, 0.0, 0.0), (6, 0.5, 0.3, 0.4), (8, 0.4, 0.2, 0.0))
        lat = Lattice(N; F = 1, m = m, θ2π = θ, mprime = mp)
        for backend in (EDBackend(), ITensorsBackend(), MPSKitBackend())
            gs = groundstate(Hamiltonian(lat; backend = backend))
            @test sum(energy_densities(gs)) ≈ energy(gs) rtol = 1e-8 atol = 1e-8
        end
    end
end
