using LinearAlgebra

# A zero-length Wilson line (start == finish) is χ†ₙχₙ, i.e. the on-site number
# operator. Check that every backend realises it as that projector — ⟨W⟩ = ⟨n⟩ and
# ‖W|gs⟩‖² = ⟨n⟩ (since n² = n) — independently of the `conjugate` flag.
@testset "Zero-length Wilson line" begin
    lat = Lattice(4; F = 1, m = 0.7)
    n = 2
    for B in (EDBackend(), ITensorsBackend(), MPSKitBackend())
        gs  = groundstate(Hamiltonian(lat; backend = B))
        occ = real(occupations(gs)[n, 1])
        for conj in (false, true)
            ψ = act(WilsonLine(lat, conj, 1, n, n; backend = B), gs)
            @test real(dot(gs, ψ)) ≈ occ      rtol = 1e-6
            @test real(dot(ψ, ψ)) ≈ occ       rtol = 1e-6
        end
    end
end
