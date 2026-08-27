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
#
# For |finish-start| ≥ 2 the Wilson line is the fermion bilinear χ†ᵢχⱼ, carrying the
# Jordan–Wigner σᶻ string ((-1)^n on the intervening sites). ED works in the occupation
# basis (occupied → −1 is unambiguous); the MPSKit string negates the site's occupied
# *charge sector*, which is staggered (charge 0 on odd sites, +q on even sites). The
# longer lines below straddle both parities of intervening site, so a wrong staggered
# sign in any single backend breaks agreement here.
@testset "Wilson line backend phase agreement" begin
    for Nsites in (4, 6)
        lat = Lattice(Nsites; F = 1, m = 0.7)
        gsED = groundstate(Hamiltonian(lat; backend = EDBackend()))
        gsIT = groundstate(Hamiltonian(lat; backend = ITensorsBackend()))
        gsMP = groundstate(Hamiltonian(lat; backend = MPSKitBackend()))
        w(gs, B, c, i, j) = dot(gs, act(WilsonLine(lat, c, 1, i, j; backend = B), gs))
        pairs = Nsites == 4 ? ((1,3), (2,4), (1,4)) :          # lengths 2, 3
                              ((1,3), (2,4), (1,4), (1,5), (2,6), (1,6))  # up to length 5
        for (i, j) in pairs, c in (false, true)
            e = w(gsED, EDBackend(), c, i, j)
            @test w(gsIT, ITensorsBackend(), c, i, j) ≈ e  rtol = 1e-5 atol = 1e-8
            @test w(gsMP, MPSKitBackend(),  c, i, j) ≈ e  rtol = 1e-5 atol = 1e-8
        end
    end
end

# On a periodic (translation-invariant) lattice the two staggered-sublattice Wilson lines of
# equal length n are related by the parity × charge-conjugation symmetry of the vacuum:
#
#     ⟨W_{0,n}⟩ = (-1)^{n-1} conj(⟨W_{1,n+1}⟩)   (bare χ†_i χ_j, no phase convention),
#
# where W_{0,n} starts on sublattice 0 (odd site 1) and W_{1,n+1} on sublattice 1 (even site 2),
# both of length n. The package's WilsonLine carries the factor im^(finish-start) = i^n (the same
# for both lines), which turns the bare (-1)^{n-1}conj into a plain −conj on the raw operator:
#
#     ⟨W(1,1+n)⟩ = −conj(⟨W(2,2+n)⟩)             (raw package operator).
#
# i.e. the real parts are opposite and the imaginary parts equal. Only ED and ITensors support
# periodic lattices (MPSKit does not); we also check they agree, incl. sign/phase.
@testset "Wilson line sublattice identity (periodic)" begin
    lat  = Lattice(8; F = 1, m = 0.7, periodic = true)
    gsED = groundstate(Hamiltonian(lat; backend = EDBackend()))
    gsIT = groundstate(Hamiltonian(lat; backend = ITensorsBackend()))
    for (gs, B) in ((gsED, EDBackend()), (gsIT, ITensorsBackend()))
        w(i, j) = dot(gs, act(WilsonLine(lat, false, 1, i, j; backend = B), gs))
        for n in 1:4
            A  = w(1, 1 + n)                       # W_{0,n}: sublattice-0 start (odd site 1)
            Bv = w(2, 2 + n)                       # W_{1,n+1}: sublattice-1 start (even site 2)
            # rtol 1e-5 (not 1e-6): the identity is exact in the true vacuum; the residual is
            # ground-state numerical precision, amplified in the relative error of this ~1e-2 value.
            @test A ≈ -conj(Bv) rtol = 1e-5 atol = 1e-8                        # raw operator
        end
    end
    # cross-backend agreement on the periodic lattice (ED vs ITensors), including sign/phase
    for (i, j) in ((1,3), (2,4), (1,4), (1,5), (2,6)), c in (false, true)
        e = dot(gsED, act(WilsonLine(lat, c, 1, i, j; backend = EDBackend()),      gsED))
        t = dot(gsIT, act(WilsonLine(lat, c, 1, i, j; backend = ITensorsBackend()), gsIT))
        @test t ≈ e rtol = 1e-5 atol = 1e-8
    end
end

# On an infinite lattice `WilsonLine` returns an `InfiniteMPOHamiltonian` whose expectation value
# is the two-site-translation-averaged Wilson line — the length-ℓ = finish-start fermion bilinear
# χ†χ (with JW string) averaged over the two sublattice start sites of the unit cell,
# (⟨W_{0,ℓ}⟩ + ⟨W_{1,ℓ}⟩)/2. Check it against a large finite lattice's mid-lattice sublattice
# average, and that it is purely imaginary (a consequence of the parity×charge-conjugation
# sublattice identity W_0 = -conj(W_1)).
@testset "Wilson line infinite lattice (translation-averaged)" begin
    lat_inf = Lattice(Inf; F = 1, m = 0.7)
    lat_fin = Lattice(20;  F = 1, m = 0.7)
    gs_inf = groundstate(Hamiltonian(lat_inf, MPSKitBackend()); bonddim = 40)
    gs_fin = groundstate(Hamiltonian(lat_fin, MPSKitBackend()); bonddim = 40)
    nrm = real(dot(gs_fin, gs_fin))
    wfin(i, j) = dot(gs_fin, act(WilsonLine(lat_fin, false, 1, i, j; backend = MPSKitBackend()), gs_fin)) / nrm

    for ℓ in (1, 2, 3)
        infavg = expectation(WilsonLine(lat_inf, false, 1, 1, 1 + ℓ; backend = MPSKitBackend()), gs_inf)
        finavg = (wfin(9, 9 + ℓ) + wfin(10, 10 + ℓ)) / 2
        @test infavg ≈ finavg rtol = 1e-2 atol = 1e-4     # matches bulk of the finite lattice
        @test abs(real(infavg)) < 1e-6                    # purely imaginary (sublattice identity)
    end

    # a finite `finish` must be passed explicitly (absolute position is meaningless when N = ∞)
    @test_throws ArgumentError WilsonLine(lat_inf; backend = MPSKitBackend())            # no finish
    @test_throws ArgumentError WilsonLine(lat_inf, false, 1, 1, 1; backend = MPSKitBackend())  # ℓ = 0
end
