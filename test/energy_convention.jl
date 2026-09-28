# §1.1: `energy_densities(...; convention=:bond)` is the density that partners `EnergyCurrent` in
# the lattice continuity ∂_t h_n = 𝒥_n − 𝒥_{n+1}; the default (:site) is NOT (it differs by a
# lattice total derivative). Also checks the :site sum rule is preserved.
using Schwinger, LinearAlgebra, Test

@testset "energy_density conventions" begin
    N   = 12
    lat = Lattice(N; F = 1, q = 2, a = 0.25, m = 1.0)     # mprime = 0
    H   = Hamiltonian(lat, EDBackend())

    # --- :site sum rule preserved (Σ site density · a = total energy) ---
    gs = groundstate(H)
    @test sum(energy_densities(gs)) * lat.a ≈ energy(gs) atol=1e-9
    @test length(energy_densities(gs))              == N
    @test length(energy_densities(gs; convention = :bond)) == N - 1

    # --- continuity: evolve a non-eigenstate, finite-difference each density's ḣ_n ---
    psiM = groundstate(Hamiltonian(Lattice(N; F = 1, q = 2, a = 0.25, m = 0.3), EDBackend()))
    pe   = evolve(EDState(H, psiM.coeffs, psiM.defects, psiM.net_charge), 0.3; nsteps = 15)[1]
    dt   = 0.005
    ep = evolve(pe, dt)[1]; em = evolve(pe, -dt)[1]
    kw = (; L_max = H.L_max, universe = H.universe, charge = 0)

    interior = 3:N-4
    Jdiv = [real(expectation(EDEnergyCurrent(lat, n; kw...), pe)) -
            real(expectation(EDEnergyCurrent(lat, n+1; kw...), pe)) for n in interior]

    # bond-centered (extensive): ∂_t h_n should equal 𝒥_n − 𝒥_{n+1}
    hb_p = energy_densities(ep; convention = :bond)
    hb_m = energy_densities(em; convention = :bond)
    hbdot = [(hb_p[n] - hb_m[n]) / (2dt) for n in interior]
    @test norm(hbdot .- Jdiv) / norm(Jdiv) < 1e-2          # probe measured 2.8e-3

    # site-centered (extensive = a·density): does NOT satisfy the same continuity
    hs_p = lat.a .* energy_densities(ep)
    hs_m = lat.a .* energy_densities(em)
    hsdot = [(hs_p[n] - hs_m[n]) / (2dt) for n in interior]
    @test norm(hsdot .- Jdiv) / norm(Jdiv) > 1e-1          # probe measured ≈ 1.0 (the trap)

    # boundary / range guards
    @test_throws ArgumentError energy_density(gs, N; convention = :bond)   # bond N-1 is the last
    @test_throws ArgumentError energy_densities(gs; convention = :nonsense)
end
