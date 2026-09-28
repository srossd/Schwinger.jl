# Tests for the matter momentum density  p_n = (-i/4a)(χ†_n U_n U_{n+1} χ_{n+2} − h.c.),  T⁰¹,
# and its total  P = Σ_n p_n.  Built from the length-2 Wilson line WilsonLine(n → n+2).
using Schwinger, LinearAlgebra, Test           # ITensors reached via Schwinger.ITensors (not a test dep)

@testset "Momentum density" begin
    N    = 12
    lat  = Lattice(N; F = 1, q = 2, a = 0.25, m = 0.5)         # mprime = 0
    c    = N ÷ 2
    θp   = [(c-2 ≤ n < c+2) ? 1.0 : 0.0 for n in 1:N]           # parity-breaking θ-step
    latθ = Lattice(N; F = 1, q = 2, a = 0.25, m = 0.5, θ2π = θp)
    D    = 64

    # (1) ground-state parity: the momentum density vanishes on the reflection-symmetric vacuum.
    for be in (EDBackend(), MPSKitBackend())
        gs = groundstate(Hamiltonian(lat, be); (be isa EDBackend ? (;) : (; energy_tol = 1e-10))...)
        @test maximum(abs(real(expectation(MomentumDensity(lat, n; backend = be), gs))) for n in 1:N-2) < 1e-6
    end

    # (2) θ-quench moving state (nonzero momentum): ED == MPSKit operator, window 3-site contraction
    #     == per-site operator, and ⟨p_n⟩ is real (Hermitian).  Build operators on lattice(state)=latθ
    #     (momentum density is θ-independent; only the operator/state lattices must match).
    gsED = groundstate(Hamiltonian(lat, EDBackend()))
    psED = evolve(EDState(Hamiltonian(latθ, EDBackend()), gsED.coeffs, gsED.defects, gsED.net_charge),
                  0.4; nsteps = 20)[1]
    pED  = [expectation(EDMomentumDensity(latθ, n), psED) for n in 1:N-2]
    @test maximum(abs(imag(p)) for p in pED) < 1e-10                       # Hermitian ⇒ real
    @test sum(abs, real.(pED)) > 1e-2                                      # genuinely nonzero

    gsMK   = groundstate(Hamiltonian(lat, MPSKitBackend()); energy_tol = 1e-10)
    psMK   = evolve(MPSKitState(Hamiltonian(latθ, MPSKitBackend()), gsMK.psi), 0.4;
                    nsteps = 20, two_site = true, maxlinkdim = D)[1]
    pMKop  = [real(expectation(MPSKitMomentumDensity(latθ, n), psMK)) for n in 1:N-2]
    pMKwin = momentumdensities(psMK)                                       # local 3-site contraction
    @test maximum(abs(pMKwin[n] - pMKop[n]) for n in 1:N-2) < 1e-8         # window == operator
    @test isapprox(sum(pMKwin), totalmomentum(psMK); atol = 1e-10)
    @test maximum(abs(real(pED[n]) - pMKop[n]) for n in 1:N-2) < 1e-3      # ED == MPSKit (truncation)

    # (3) ITensors via an evolution-free momentum boost U = exp(iκ Σ_n n N_n) (diagonal ⇒ gauge-
    #     independent, parity-breaking): applied identically to the ED and ITensors ground states,
    #     ED-boost and ITensors-boost momentum densities must agree.  (θ2π-vector TDVP OOMs in ITensors.)
    κb   = 0.3
    phED = [exp(im*κb*sum(n*bs.occupations[n,1] for n in 1:N)) for bs in Schwinger._edbasis(gsED)]
    bED  = EDState(gsED.hamiltonian, gsED.coeffs .* phED, gsED.defects, gsED.net_charge)
    pEDb = [real(expectation(EDMomentumDensity(lat, n), bED)) for n in 1:N-2]

    gsIT = groundstate(Hamiltonian(lat, ITensorsBackend()))
    sIT  = Schwinger.get_sites(gsIT.hamiltonian)     # same site indices the MPS was built from
    psib = copy(gsIT.psi)
    for n in 1:N
        s = sIT[n]
        projUp = 0.5*Schwinger.ITensors.op("Id", s) + Schwinger.ITensors.op("Sz", s)  # N_n = Sz+1/2 (Up=occupied)
        g = Schwinger.ITensors.op("Id", s) + (exp(im*κb*n) - 1)*projUp                # boost gate exp(iκ n N_n)
        psib[n] = Schwinger.ITensors.noprime(g * psib[n])
    end
    pITb = [real(expectation(ITensorMomentumDensity(lat, n),
                             ITensorState(gsIT.hamiltonian, psib, gsIT.defects))) for n in 1:N-2]
    @test sum(abs, pEDb) > 1e-3                                            # boost produced momentum
    @test maximum(abs(pEDb[n] - pITb[n]) for n in 1:N-2) < 1e-5            # ED == ITensors

    # (4) infinite lattice: a quasiparticle wavepacket carries the momentum it was built at.  The
    #     total momentum is odd in p, has the right sign, and matches the free lattice dispersion
    #     ½a·sin(2ap) — hence ≈ p at small a·p (continuum).  (Small a so a·p is small; modest bond dim.)
    latinf = Lattice(Inf; F = 1, a = 0.2, m = 1.0)
    Hinf   = Hamiltonian(latinf; backend = :MPSKit)
    res    = loweststates(Hinf, 2; bonddim = 24, energy_tol = 1e-9, momentum = [-0.5, 0.5])
    qm, qp = res[2]
    Wi, sup = 64, 17:48
    x0 = (first(sup) + last(sup)) / 2
    Pp = totalmomentum(wavepacket(qp, Wi; support = sup, sigma = 8.0, center = x0))
    Pm = totalmomentum(wavepacket(qm, Wi; support = sup, sigma = 8.0, center = x0))
    ap = 0.2 * 0.5
    @test Pp > 0 && Pm < 0                                  # right / left movers
    @test isapprox(Pp, -Pm; rtol = 0.03)                   # odd in p
    @test isapprox(Pp, sin(2ap)/(2*0.2); rtol = 0.06)      # free lattice dispersion ½a·sin(2ap)
    @test isapprox(Pp, 0.5; rtol = 0.06)                   # ≈ physical momentum (small a·p)
end
