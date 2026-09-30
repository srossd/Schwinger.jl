# apply_wilsonline / apply_local(::MPSKitOperator): apply a bounded-support Wilson line to a window.
# Validation: ⟨vac| W(i→i+ℓ) |vac⟩ measured on a vacuum window, averaged over the odd- and even-start
# sublattices, must equal the infinite Wilson-line automaton expectation on the vacuum (which the
# code documents as exactly that two-site-translation average).
using Schwinger, LinearAlgebra, Test
using MPSKit, TensorKit

@testset "apply_wilsonline on a window" begin
    lat = Lattice(Inf; F = 1, q = 2, a = 0.25, m = 0.5)
    H   = Hamiltonian(lat, MPSKitBackend())
    gsi = groundstate(H; bonddim = 16)
    W   = 20
    vac = MPSKitState(H, WindowMPS(gsi.psi, W))
    nrm2 = real(dot(vac, vac))

    ℓ = 2
    Wodd  = apply_wilsonline(vac, 9,  9 + ℓ)     # odd start
    Weven = apply_wilsonline(vac, 10, 10 + ℓ)    # even start
    @test Wodd.psi  isa WindowMPS && length(Wodd.psi)  == W   # wings preserved, bonds intact
    @test Weven.psi isa WindowMPS && length(Weven.psi) == W

    val_odd  = dot(vac, Wodd)  / nrm2
    val_even = dot(vac, Weven) / nrm2
    val_window_avg = (val_odd + val_even) / 2

    # infinite Wilson-line automaton expectation on the vacuum = the odd/even-start average
    val_inf = expectation(WilsonLine(lat, false, 1, 1, 1 + ℓ; backend = MPSKitBackend()), gsi)

    @test isapprox(val_window_avg, val_inf; rtol = 1e-4)
    @test abs(imag(val_window_avg)) < 1e-6            # neutral bilinear ⇒ real on the vacuum

    # a length-1 (nearest-neighbour) line also works and stays a window
    W1 = apply_wilsonline(vac, 9, 10)
    @test W1.psi isa WindowMPS

    # guard: an infinite-lattice operator is rejected by apply_local(state, op)
    @test_throws ArgumentError apply_local(vac, H)   # H is the whole-lattice (infinite) operator
end
