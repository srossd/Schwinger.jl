# Tests for the conserved lattice currents ChargeCurrent (j¹) and EnergyCurrent (T⁰¹ = 𝒥).
using Schwinger, LinearAlgebra, Test

@testset "Charge/energy currents" begin
    N = 12
    lat = Lattice(N; F = 1, q = 2, a = 0.25, m = 0)          # mprime = 0 required
    H   = Hamiltonian(lat, EDBackend())
    kw  = (; L_max = H.L_max, universe = H.universe, charge = 0)

    # --- ED: exact energy-current operator identity  i[H,h_n] = 𝒥_n − 𝒥_{n+1} ---
    hfull(m) = EDGaugeKinetic(lat, m; bare = false, kw...) +
               EDHopping(lat, m; bare = false, kw...) +
               0.5*EDMass(lat, m; bare = false, kw...) + 0.5*EDMass(lat, m+1; bare = false, kw...)
    Hm = H.matrix
    for n in 3:N-3
        lhs = im*(Hm*hfull(n).matrix - hfull(n).matrix*Hm)
        rhs = EDEnergyCurrent(lat, n; kw...).matrix - EDEnergyCurrent(lat, n+1; kw...).matrix
        @test maximum(abs, lhs - rhs) < 1e-9
    end

    # --- ED: charge- & energy-current continuity vs a two-step time derivative, on a quench ---
    # The current is the flux in the continuity equation, so measuring the current at t must match the
    # time derivative of the corresponding density, estimated from a symmetric difference of two
    # evolve steps (t ± dt):  ρ̇_n = -(j¹_n - j¹_{n-1})  and  ḣ_n = 𝒥_n - 𝒥_{n+1}.
    psi0 = groundstate(H)
    c = N ÷ 2
    psiq = insert_defect(insert_defect(psi0, c-2, +1), c+2, -1)
    Hq   = Hamiltonian(lat, EDBackend(); defects = defects(psiq))
    psit = evolve(EDState(Hq, psiq.coeffs, psiq.defects, psiq.net_charge), 0.5; nsteps = 25)[1]
    dt = 0.005
    sp = evolve(psit, dt)[1]; sm = evolve(psit, -dt)[1]           # t ± dt (once, reused below)

    # charge:  ρ̇_n = -(j¹_n - j¹_{n-1})
    j = [real(expectation(EDChargeCurrent(lat, b; kw...), psit)) for b in 1:N-1]
    rhodot = (real.(vec(charges(sp))) .- real.(vec(charges(sm)))) ./ (2dt)
    div = zeros(N); div[1] = j[1]; for n in 2:N-1; div[n] = j[n]-j[n-1]; end; div[N] = -j[N-1]
    @test norm(rhodot .+ div) / norm(rhodot) < 2e-3

    # energy:  ḣ_n = 𝒥_n - 𝒥_{n+1}   (h_n = the full local energy density hfull(n)).
    # Evolve a non-eigenstate (the m = 1 ground state) under the SAME defect-free H that 𝒥/hfull are
    # built from, so the continuity holds exactly. (The charge quench above evolves under Hq, whose
    # static defect charges shift the electric-field energy via Gauss law — fine for charge continuity
    # since they commute with fermion number, but the defect-free 𝒥 would not capture that here.)
    psiM = groundstate(Hamiltonian(Lattice(N; F = 1, q = 2, a = 0.25, m = 1.0), EDBackend()))
    pe   = evolve(EDState(H, psiM.coeffs, psiM.defects, psiM.net_charge), 0.3; nsteps = 15)[1]
    ep = evolve(pe, dt)[1]; em = evolve(pe, -dt)[1]
    hdot = [(real(expectation(hfull(n), ep)) - real(expectation(hfull(n), em))) / (2dt) for n in 3:N-4]
    Jdiv = [real(expectation(EDEnergyCurrent(lat, n;   kw...), pe)) -
            real(expectation(EDEnergyCurrent(lat, n+1; kw...), pe)) for n in 3:N-4]
    @test norm(hdot .- Jdiv) / norm(Jdiv) < 5e-3

    # --- MPSKit: charge-current continuity ρ̇ = -(j_n - j_{n-1}) via two evolve steps (θ-string quench) ---
    # Exercises the memory-light `chargecurrents` observable (local 2-site contraction, not the operator).
    let be = MPSKitBackend(), dtm = 0.01, D = 48
        θp   = [(c-2 ≤ n < c+2) ? 1.0 : 0.0 for n in 1:N]
        latθ = Lattice(N; F = 1, q = 2, a = 0.25, m = 0, θ2π = θp)
        gsm  = groundstate(Hamiltonian(lat, be); energy_tol = 1e-10)
        stm  = evolve(MPSKitState(Hamiltonian(latθ, be), gsm.psi), 0.4; nsteps = 8, two_site = true, maxlinkdim = D)[1]
        jm   = chargecurrents(stm)                               # local 2-site path (memory-light)
        ρp   = real.(vec(charges(evolve(stm,  dtm; two_site = true, maxlinkdim = D)[1])))
        ρm   = real.(vec(charges(evolve(stm, -dtm; two_site = true, maxlinkdim = D)[1])))
        ρdot = (ρp .- ρm) ./ (2dtm)
        divj = zeros(N); divj[1] = jm[1]; for n in 2:N-1; divj[n] = jm[n]-jm[n-1]; end; divj[N] = -jm[N-1]
        @test norm(ρdot .+ divj) / norm(ρdot) < 3e-2             # MPSKit truncation ⇒ looser tol
    end

    # --- all backends: currents vanish on the (parity-symmetric) ground state ---
    for be in (EDBackend(), ITensorsBackend(), MPSKitBackend())
        gs = groundstate(Hamiltonian(lat, be); (be isa EDBackend ? (;) : (; energy_tol = 1e-10))...)
        @test maximum(abs(real(expectation(ChargeCurrent(lat, b; backend = be), gs))) for b in 2:N-2) < 1e-6
        @test maximum(abs(real(expectation(EnergyCurrent(lat, s; backend = be), gs))) for s in 3:N-2) < 1e-6
    end
end
