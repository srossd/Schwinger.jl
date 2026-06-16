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

# All three backends share ED's imaginary-hopping gauge: ITensors/MPSKit carry the
# factor i^(finish-start) so that ⟨gs|W(i,j)|gs⟩ (independent of the ground-state
# phase) agrees across backends, for both `conjugate` values.
@testset "Wilson line backend phase agreement" begin
    lat = Lattice(4; F = 1, m = 0.7)
    gsED = groundstate(Hamiltonian(lat; backend = EDBackend()))
    gsIT = groundstate(Hamiltonian(lat; backend = ITensorsBackend()))
    gsMP = groundstate(Hamiltonian(lat; backend = MPSKitBackend()))
    w(gs, B, c, i, j) = dot(gs, act(WilsonLine(lat, c, 1, i, j; backend = B), gs))
    for (i, j) in ((1,3), (2,4), (1,4)), c in (false, true)   # lengths 2 and 3
        e = w(gsED, EDBackend(), c, i, j)
        @test w(gsIT, ITensorsBackend(), c, i, j) ≈ e  rtol = 1e-5 atol = 1e-8
        @test w(gsMP, MPSKitBackend(),  c, i, j) ≈ e  rtol = 1e-5 atol = 1e-8
    end
end
