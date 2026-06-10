using LinearAlgebra

# The single staggered fermion field φ̃_s. We check two backend-independent properties:
#  * ‖φ̃_s|gs⟩‖² = ⟨n_s⟩  (the Jordan–Wigner string is unitary);
#  * the gauge-invariant Δz=0 light-cone matrix element agrees across backends, which
#    pins the relative phase between fields on different sites (each backend's hopping
#    convention requires a matching string phase).
@testset "Fermion field" begin
    lat = Lattice(4; F=1, a=0.5, m=0.7)
    c = 2 * (Int(lat.N) ÷ 4)

    # M(0) = Σ (-1)^{f+s} [⟨φ̃_{c+f}h|φ̃_{c+s}h⟩ − ⟨φ̃_{c+f}0|φ̃_{c+s}0⟩]
    function fieldcheck(backend)
        gs, h = backend isa EDBackend ? loweststates(Hamiltonian(lat, backend), 2) :
                                        loweststates(Hamiltonian(lat, backend), 2; energy_tol=1e-10)
        occ = occupations(gs)
        for s in (c, c+1)
            ψ = act(FermionField(lat, s; backend=backend), gs)
            @test real(dot(ψ, ψ)) ≈ occ[s, 1] rtol=1e-5
        end
        amp(ket) = [dot(act(FermionField(lat, c+f; backend=backend), ket),
                        act(FermionField(lat, c+s; backend=backend), ket)) for f in 0:1, s in 0:1]
        Ah, A0 = amp(h), amp(gs)
        return sum((-1)^(f+s) * (Ah[f+1, s+1] - A0[f+1, s+1]) for f in 0:1, s in 0:1)
    end

    M0 = Dict(b => fieldcheck(b) for b in (EDBackend(), ITensorsBackend(), MPSKitBackend()))
    @test M0[EDBackend()]       ≈ M0[MPSKitBackend()] rtol=1e-3 atol=1e-4
    @test M0[ITensorsBackend()] ≈ M0[MPSKitBackend()] rtol=1e-3 atol=1e-4
end
