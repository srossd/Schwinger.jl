using LinearAlgebra, TensorKit, MPSKit

# The reflection-symmetric gauge (arXiv:2012.07243, eq. C17) fixes the residual Bloch
# gauge freedom on the quasiparticle B tensors by minimising the symmetric objective
#   f(B) = Σ_i ‖Σ_s B[i]^s ⊗ (AL[i]^s)*‖² + ‖Σ_s B[i]^s ⊗ (AR[i]^s)*‖².
@testset "Reflection-symmetric QP gauge" begin
    lat = Lattice(Inf; F = 1, m = 0.5)
    # m = 0.5 is gapped and converges fast; converge well so the multi-particle domain
    # wall (a center-matrix pseudo-inverse) is not corrupted by stray Schmidt values.
    gs, qpst = loweststates(Hamiltonian(lat; backend = :MPSKit), 2;
                            bonddim = 10, maxiters = 300, momentum = 0.3)
    qp = qpst.psi
    L  = length(qp)
    AL = [qp.left_gs.AL[i]  for i in 1:L]
    AR = [qp.right_gs.AR[i] for i in 1:L]

    function term(B, A, i)                        # ‖Σ_s B[i]^s ⊗ (A[i]^s)*‖²
        @tensor t[-1 -2; -3 -4 -5] := A[i][-1 1; -2] * conj(B[i][-3 1; -4 -5])
        return real(dot(t, t))
    end
    fL(B) = sum(term(B, AL, i) for i in 1:L)
    fR(B) = sum(term(B, AR, i) for i in 1:L)

    Bleft = [qp[i] for i in 1:L]
    Bsym  = reflection_symmetric_gauge(qp)

    # Optimality: the symmetric gauge lowers the objective well below either canonical gauge.
    @test fL(Bsym) + fR(Bsym) < 0.5 * (fL(Bleft) + fR(Bleft))

    # Gauge-orbit invariance: the symmetric tensors are a property of the physical state,
    # so re-gauging the right-gauge representation of the *same* state must agree.
    Bsym2 = reflection_symmetric_gauge(convert(MPSKit.RightGaugedQP, qp))
    @test maximum(norm(Bsym[i] - Bsym2[i]) for i in 1:L) < 1e-7

    # End-to-end: the wavepacket builds a valid (finite, nonzero norm) state in either
    # gauge. (The two localized states differ by gauge — that is precisely the freedom
    # being fixed — so their norms need not agree.)
    wsym  = wavepacket(qpst, 20; support = 5:16, gauge = :symmetric)
    wleft = wavepacket(qpst, 20; support = 5:16, gauge = :left)
    @test isfinite(norm(wsym.psi))  && norm(wsym.psi)  > 0
    @test isfinite(norm(wleft.psi)) && norm(wleft.psi) > 0

    # energy_densities works on the wavepacket WindowMPS (full per-site energy via
    # link_expectation + local terms).
    ed = energy_densities(wsym)
    @test length(ed) == 20
    @test all(isfinite, ed)

    # Multi-particle: the mixed-canonical background centres EVERY QP identically. Build a
    # two-particle wavepacket (same QP at two disjoint supports) and check each excitation
    # is centred on its support and the inter-particle gap returns to the vacuum.
    W2   = 48
    s1, s2 = 9:16, 33:40
    w2   = wavepacket([qpst, qpst], W2; supports = [s1, s2], gauge = :symmetric)
    ed2  = energy_densities(w2)
    vac  = ed2[3]                                    # far-field site is vacuum
    exc  = ed2 .- vac
    centroid(s) = sum(k * max(exc[k], 0) for k in s) / sum(max(exc[k], 0) for k in s)
    @test abs(centroid(s1) - (first(s1) + last(s1)) / 2) < 1.0   # centred on support 1
    @test abs(centroid(s2) - (first(s2) + last(s2)) / 2) < 1.0   # centred on support 2
    @test maximum(abs, exc[22:26]) < 1e-2                        # gap returns to vacuum
    # the two QPs are treated identically → matching profiles
    @test abs(maximum(exc[s1]) - maximum(exc[s2])) < 1e-3
    @test maximum(abs, exc[s1] .- exc[s2]) < 1e-3
end
