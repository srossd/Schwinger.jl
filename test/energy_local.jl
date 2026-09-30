# Request 1: the fast local 3-site `energycurrents` (contract_mpo_expval3) must equal the reference
# per-site `MPSKitEnergyCurrent` operator it replaces, on a finite lattice.
using Schwinger, LinearAlgebra, Test

@testset "MPSKit energycurrents (local 3-site path)" begin
    N   = 10
    lat = Lattice(N; F = 1, q = 2, a = 0.25, m = 0.7)     # mprime = 0
    H   = Hamiltonian(lat, MPSKitBackend())
    gs  = groundstate(H; energy_tol = 1e-10)
    u   = H.universe

    # local 3-site contraction == the per-site MPSKitEnergyCurrent operator.
    # Tolerance 1e-6 (not machine): the two paths contract the SAME operator but in different orders,
    # and the DMRG ground state's local observables converge to ~1e-7 (energy converges faster than
    # densities — see the loweststates docstring), so the agreement floors at ~1e-7. A wrong
    # coefficient/sign in the local operator would show as an O(1) discrepancy, not 1e-7.
    Jloc = energycurrents(gs)
    Jop  = [real(expectation(MPSKitEnergyCurrent(lat, s; universe = u), gs)) for s in 2:N-1]
    @test length(Jloc) == N - 2
    @test norm(Jloc - Jop) / max(norm(Jop), 1e-12) < 1e-6

    # pad option: site-aligned length-N with NaN at the two boundary sites
    ecp = energycurrents(gs; pad = true)
    @test length(ecp) == N && isnan(ecp[1]) && isnan(ecp[N]) && ecp[2:N-1] == Jloc

    # energy_densities :site sum rule still holds (per-site path)
    @test sum(energy_densities(gs)) * lat.a ≈ energy(gs) atol = 1e-8
end
