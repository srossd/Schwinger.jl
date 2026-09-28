# §1.3: guard against a silent convention split in `pseudoscalardensities`.
#   • the MPSKitState override (bond-contraction path) must equal the generic per-site operator path
#   • ED / ITensors / MPSKit must agree on the unique finite-lattice ground state (same physical state)
# Boundaries differ BY DESIGN (the generic per-site path wraps `before = N` for site 1; the override
# uses half-weight open boundaries), so the interior 2..N-1 is compared.
using Schwinger, LinearAlgebra, Test

@testset "pseudoscalardensities conventions agree" begin
    N   = 10
    lat = Lattice(N; F = 1, q = 2, a = 0.25, m = 0.7)
    int = 2:N-1

    gsM = groundstate(Hamiltonian(lat, MPSKitBackend()); energy_tol = 1e-10)
    gsE = groundstate(Hamiltonian(lat, EDBackend()))
    gsI = groundstate(Hamiltonian(lat, ITensorsBackend()); energy_tol = 1e-10)

    pd_override = pseudoscalardensities(gsM)                 # MPSKitState method (bond contraction)
    pd_persite  = [pseudoscalardensity(gsM, s) for s in 1:N] # generic per-site operator path

    # the two MPSKit methods must be the same quantity (docstring claims machine precision)
    @test pd_override[int] ≈ pd_persite[int] atol=1e-10

    # cross-backend: same physical ground state ⇒ same interior profile
    @test pseudoscalardensities(gsE)[int] ≈ pd_persite[int] atol=1e-6
    @test pseudoscalardensities(gsI)[int] ≈ pd_persite[int] atol=1e-6

    # sanity: the profile is non-trivial (so the agreement is meaningful, not all-zeros)
    @test norm(pd_persite[int]) > 1e-2
end
