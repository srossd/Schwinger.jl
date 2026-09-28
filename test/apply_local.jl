# §2.1 (core): `apply_local` applies a single-site operator, editing only the window tensor at
# `site` and leaving the infinite wings untouched — so it works on a wavepacket `WindowMPS`
# (which `act(::MPSKitOperator, …)` cannot) and on a bare `FiniteMPS` window. Validated against
# `occupations`: with N_k the site number operator, ⟨ψ|N_k|ψ⟩/‖ψ‖² == occupations(ψ)[k].
using Schwinger, LinearAlgebra, Test
using TensorKit, MPSKit

# site number operator on physical space P (matches the convention in `occupations`)
function numop_on(P, site::Int, q::Int)
    op = zeros(ComplexF64, P ← P)
    block(op, U1Irrep(isodd(site) ? 0 : q)) .= 1.0
    return op
end

@testset "apply_local" begin
    q = 2
    lat = Lattice(Inf; F = 1, q = q, a = 0.25, m = 0.5)
    # one infinite VUMPS solve, reused for both window flavors (finite MPSKit groundstate is slow)
    gs, qpst = loweststates(Hamiltonian(lat, MPSKitBackend()), 2;
                            bonddim = 10, maxiters = 300, momentum = 0.3)

    # window = true  → WindowMPS (explicit infinite wings); window = false → bare FiniteMPS window
    for windowed in (true, false)
        w = wavepacket(qpst, 20; support = 5:16, gauge = :symmetric, window = windowed)
        @test w.psi isa (windowed ? WindowMPS : FiniteMPS)
        occ  = occupations(w)[:, 1]
        nrm2 = real(dot(w, w))
        for k in (6, 9, 12)
            P  = TensorKit.space(w.psi.AC[k], 2)
            Nw = apply_local(w, numop_on(P, k, q), k)
            @test Nw.psi isa (windowed ? WindowMPS : FiniteMPS)      # backing type preserved
            @test length(Nw.psi) == length(w.psi)                   # bond structure untouched
            @test real(dot(w, Nw)) / nrm2 ≈ occ[k] atol = 1e-7      # ⟨N_k⟩
            @test real(dot(Nw, Nw)) / nrm2 ≈ occ[k] atol = 1e-7     # N_k projector: ⟨N²⟩ = ⟨N⟩
        end
        # error paths
        @test_throws ArgumentError apply_local(w, numop_on(TensorKit.space(w.psi.AC[1], 2), 1, q), 0)
        @test_throws ArgumentError apply_local(w, numop_on(TensorKit.space(w.psi.AC[1], 2), 1, q), 21)
    end
end
