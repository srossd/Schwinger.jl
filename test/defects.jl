using ProgressMeter
using LinearAlgebra

# Static defect charges: a `DefectCharge` must reproduce the equivalent θ2π step.
#  * ED / ITensors realise it as that θ step (no extra site);
#  * MPSKit inserts a genuine extra lattice site that is invisible to observables.
# We check (i) every backend's defect build against the "direct θ" alternative,
# (ii) the three backends against each other, and (iii) observable invisibility.

@testset "Defect charges" begin
    # (N, m, θ2π_base, defects)
    cases = [
        (8, 0.5, 0.0, [DefectCharge(5, 1)]),
        (8, 0.5, 0.2, [DefectCharge(3, 1)]),
        (8, 1.0, 0.0, [DefectCharge(3, 1), DefectCharge(6, -1)]),
        (6, 1.0, 0.0, [DefectCharge(4, 2)]),
    ]

    @showprogress desc="Defects" for (N, m, θ0, defs) in cases
        a = 0.7
        lat = Lattice(N; F=1, a=a, m=m, θ2π=θ0)

        # "direct θ" reference: θ2π[n] += charge for n ≥ link, built by hand (no defect API)
        θ = fill(Float64(θ0), N)
        for d in defs, n in d.link:N
            θ[n] += d.charge
        end
        latθ   = Lattice(N; F=1, a=a, m=m, θ2π=θ)
        ed_ref = groundstate(Hamiltonian(latθ, EDBackend()))

        # defect API on every backend
        gss = [
            (:ED,      groundstate(Hamiltonian(lat, EDBackend();        defects=defs))),
            (:ITensor, groundstate(Hamiltonian(lat, ITensorsBackend(); defects=defs), energy_tol=1e-10)),
            (:MPSKit,  groundstate(Hamiltonian(lat, MPSKitBackend();    defects=defs); energy_tol=1e-10)),
        ]

        for (name, gs) in gss
            # uniform across backends: original (unaltered) lattice + the defect list
            @test collect(Float64, lattice(gs).θ2π[1:N]) ≈ fill(Float64(θ0), N)
            @test defects(gs) == defs
            # matches the direct-θ reference (this also compares the backends to each other)
            @test energy(gs)         ≈ energy(ed_ref)         rtol=1e-6
            @test occupations(gs)    ≈ occupations(ed_ref)    rtol=1e-3 atol=1e-4
            @test electricfields(gs) ≈ electricfields(ed_ref) rtol=1e-3 atol=1e-4
            @test scalardensities(gs) ≈ scalardensities(ed_ref) rtol=1e-3 atol=1e-4
            # the inserted defect site(s) / zero-width links are invisible to observables
            @test size(occupations(gs), 1) == N
            @test length(charges(gs))       == N
            @test length(electricfields(gs)) == N
        end

        # the "direct θ" build also agrees across backends (exercises non-uniform θ)
        @test energy(groundstate(Hamiltonian(latθ, ITensorsBackend()); energy_tol=1e-10)) ≈ energy(ed_ref) rtol=1e-6
        @test energy(groundstate(Hamiltonian(latθ, MPSKitBackend());   energy_tol=1e-10)) ≈ energy(ed_ref) rtol=1e-6
    end

    # translate_defect: move a static charge, keeping the matter wavefunction fixed.
    @testset "translation" begin
        N = 8; a = 0.7; m = 0.5
        lat = Lattice(N; F=1, a=a, m=m)
        d0 = DefectCharge(3, 1); ℓ = 6

        ed   = groundstate(Hamiltonian(lat, EDBackend();     defects=[d0]))
        it   = groundstate(Hamiltonian(lat, ITensorsBackend(); defects=[d0]), energy_tol=1e-10)
        mps  = groundstate(Hamiltonian(lat, MPSKitBackend(); defects=[d0]); energy_tol=1e-10)
        edt  = translate_defect(ed,  d0, ℓ)
        itt  = translate_defect(it,  d0, ℓ)
        mpst = translate_defect(mps, d0, ℓ)

        for tr in (edt, itt, mpst)
            @test defects(tr) == [DefectCharge(ℓ, 1)]              # defect relocated
        end
        # matter wavefunction unchanged by the move
        @test occupations(edt)  ≈ occupations(ed)   atol=1e-10
        @test occupations(itt)  ≈ occupations(it)   atol=1e-8
        @test occupations(mpst) ≈ occupations(mps)  atol=1e-6     # MPS move is exact
        @test length(mpst.psi) == N + 1                           # site count preserved
        # all backends agree after translation (matter + relocated flux)
        @test occupations(mpst)    ≈ occupations(edt)    rtol=1e-3 atol=1e-4
        @test electricfields(mpst) ≈ electricfields(edt) rtol=1e-3 atol=1e-4
        @test electricfields(itt)  ≈ electricfields(edt) rtol=1e-3 atol=1e-4
        # the move genuinely changed the electric field (relocated the unit of flux)
        @test maximum(abs.(electricfields(edt) .- electricfields(ed))) > 0.5
        # MPSKit fuse/split round-trip recovers the original state exactly
        mpsrt = translate_defect(mpst, DefectCharge(ℓ, 1), 3)
        @test abs(dot(mpsrt.psi, mps.psi)) / (norm(mpsrt.psi) * norm(mps.psi)) ≈ 1 atol=1e-8
        @test occupations(mpsrt) ≈ occupations(mps) atol=1e-8
    end
end

# Static-defect time evolution uses the "absorb" trick: each inert 1-D defect site is
# fused into its matter neighbour (so it can't cap the TDVP/DMRG bond dimension), the
# pure-matter chain is evolved, and the defect sites are split back out. We check the
# fused Hamiltonian/state against the defect ones and the dynamics against ED.
@testset "Defect evolution (absorb)" begin
    N = 8; a = 0.7
    lat  = Lattice(N; F=1, a=a, m=0.5); Q = lat.q
    defs = [DefectCharge(4, Q)]

    Hd = Hamiltonian(lat, MPSKitBackend(); defects=defs)
    gd = groundstate(Hd; energy_tol=1e-10)

    # (1) the fused (defect-absorbed) LEMPO reproduces the defect Hamiltonian's energy
    flempo = Schwinger._fused_defect_lempo(lat, defs, 0)
    fop    = MPSKitOperator(lat, flempo, 0, DefectCharge[])
    fgd    = Schwinger._fuse_defect_mps(gd.psi, lat, defs)
    @test length(fgd) == N                                     # defect site absorbed away
    @test real(energy(MPSKitState(fop, fgd, DefectCharge[]))) ≈ real(energy(gd)) rtol=1e-8

    # (2) fuse ∘ split is the identity on the state
    back = Schwinger._split_defect_mps(fgd, lat, defs)
    @test abs(dot(gd.psi, back)) / (norm(gd.psi) * norm(back)) ≈ 1 atol=1e-8

    # (3) real-time evolution of a non-eigenstate defect state conserves energy and the
    #     MPSKit (absorb) dynamics agree with exact diagonalisation
    lat2 = Lattice(N; F=1, a=a, m=1.0)
    ψ0 = MPSKitState(Hamiltonian(lat2, MPSKitBackend(); defects=defs), gd.psi, defs)
    ψt, _ = evolve(ψ0, 0.4; nsteps=20, two_site=true, maxlinkdim=64)
    @test length(ψt.psi) == N + 1                             # defect site restored
    @test real(energy(ψt)) ≈ real(energy(ψ0)) rtol=1e-6       # unitary ⇒ energy conserved

    ed0 = EDState(Hamiltonian(lat2, EDBackend(); defects=defs),
                  groundstate(Hamiltonian(lat, EDBackend(); defects=defs)).coeffs, defs)
    edt, _ = evolve(ed0, 0.4; nsteps=20)
    @test occupations(ψt) ≈ occupations(edt) rtol=1e-3 atol=1e-4
end
