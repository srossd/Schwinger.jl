# Request 2: infinite ground state -> vacuum WindowMPS -> act (apply_local) -> evolve with adaptive
# window growth. Also checks window `energycurrents` (now supported) and the normalize/WindowMPS fix.
using Schwinger, LinearAlgebra, Test
using MPSKit, TensorKit

@testset "infinite gs -> WindowMPS -> act -> evolve+grow" begin
    lat = Lattice(Inf; F = 1, q = 2, a = 0.25, m = 0.5)
    H   = Hamiltonian(lat, MPSKitBackend())
    gsi = groundstate(H; bonddim = 16)
    W   = 20
    wvac = MPSKitState(H, WindowMPS(gsi.psi, W))          # vacuum window straight from the infinite gs
    @test wvac.psi isa WindowMPS

    # energycurrents on a window now works; a translation-invariant vacuum carries ~zero energy flow
    Jw = energycurrents(wvac)
    @test length(Jw) == W - 2
    @test maximum(abs, Jw) < 1e-4
    @test length(energy_densities(wvac)) == W

    # acting with a global (whole-lattice) MPSKitOperator on a window gives a clear error that
    # points to apply_local (use H, which is a validly-constructed infinite-lattice operator).
    err = try
        act(H, wvac); ""
    catch e; sprint(showerror, e) end
    @test occursin("apply_local", err)

    # apply_local a single-site operator, normalize, then evolve with adaptive growth.
    k  = 10
    P  = TensorKit.space(wvac.psi.AC[k], 2)
    op = zeros(ComplexF64, P ← P); block(op, U1Irrep(iseven(k) ? 2 : 0)) .= 1.0
    exc = normalize(apply_local(wvac, op, k))
    @test exc.psi isa WindowMPS                          # normalize must NOT downgrade the window (bug fix)

    cond = window_growth_condition(nsites = 6, threshold = 1e-3, growth = 2)
    ev, _ = evolve(exc, 0.2; nsteps = 2, two_site = true, maxlinkdim = 32, grow = cond)
    @test ev.psi isa WindowMPS                           # evolve+grow runs and stays a window
    @test length(ev.psi) ≥ W
    @test (length(ev.psi) - W) % 2 == 0                  # any growth is in whole unit cells
end
