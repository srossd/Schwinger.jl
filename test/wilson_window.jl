# Applying a bounded-support operator to a window via window_lattice + apply_local(state, op).
# Validation on a vacuum window with a length-1 Wilson line W(i→i+1) = χ†_i U χ_{i+1} (nonzero on the
# vacuum): (a) same-sublattice translation invariance ⟨W(9→10)⟩ == ⟨W(11→12)⟩, and (b) the odd/even
# start average equals the infinite Wilson-line automaton expectation (documented to be exactly that
# two-site-translation average).
using Schwinger, LinearAlgebra, Test
using MPSKit, TensorKit

@testset "apply_local(state, op) on a window (Wilson line)" begin
    lat = Lattice(Inf; F = 1, q = 2, a = 0.25, m = 0.5)
    H   = Hamiltonian(lat, MPSKitBackend())
    gsi = groundstate(H; bonddim = 16)
    W   = 20
    vac = MPSKitState(H, WindowMPS(gsi.psi, W))
    nrm2 = real(dot(vac, vac))

    latW = window_lattice(vac)
    @test latW isa Lattice
    @test Int(latW.N) == W && latW.q == lat.q && latW.F == lat.F && latW.a == lat.a

    wl(i, j) = WilsonLine(latW, false, 1, i, j; backend = MPSKitBackend())
    meas(i, j) = dot(vac, apply_local(vac, wl(i, j))) / nrm2      # ⟨vac| W(i→j) |vac⟩

    W9  = apply_local(vac, wl(9, 10))
    @test W9.psi isa WindowMPS && length(W9.psi) == W             # wings preserved, bonds intact

    v9  = dot(vac, W9) / nrm2                                     # odd start (nonzero: hopping bilinear)
    v11 = meas(11, 12)                                            # odd start, one unit cell over
    v10 = meas(10, 11)                                            # even start
    @test abs(v9) > 1e-3                                          # a genuine (nonzero) signal
    @test isapprox(v9, v11; rtol = 1e-4)                         # same-sublattice translation invariance

    # odd/even-start average == infinite Wilson-line automaton expectation on the vacuum
    val_inf = expectation(WilsonLine(lat, false, 1, 1, 2; backend = MPSKitBackend()), gsi)
    @test isapprox((v9 + v10) / 2, val_inf; atol = 1e-8, rtol = 1e-4)

    # guard: an infinite (whole-lattice) operator is rejected — must be built over window_lattice
    @test_throws ArgumentError apply_local(vac, H)
end
