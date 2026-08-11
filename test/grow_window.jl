# Adaptive window growth: `grow_window` splices exact vacuum sites onto a WindowMPS
# wavepacket (preserving the evolved content), `boundary_energy_excess` measures how much
# the energy density near each boundary exceeds the wing vacuum, and
# `window_growth_condition` turns that into an `evolve` growth callback.
using LinearAlgebra, MPSKit

@testset "Adaptive window growth" begin
    lat = Lattice(Inf; F = 1, m = 0.5)
    # m = 0.5 is gapped and converges fast (same setup as the QP-gauge tests).
    gs, qpst = loweststates(Hamiltonian(lat; backend = :MPSKit), 2;
                            bonddim = 10, maxiters = 300, momentum = 0.3)
    Uc = length(qpst.psi.left_gs)                            # background unit cell (= 2)

    @testset "grow_window preserves the physical state" begin
        w  = wavepacket(qpst, 20; support = 5:16, gauge = :symmetric, window = true)
        @test w.psi isa WindowMPS
        ed = energy_densities(w)

        left, right = 8, 8
        w2 = grow_window(w; left = left, right = right)
        @test w2.psi isa WindowMPS
        @test length(w2.psi) == 20 + left + right

        # Isometric vacuum padding leaves the norm unchanged...
        @test norm(w2.psi) ≈ norm(w.psi) rtol = 1e-10
        # ...and the interior energy densities are identical, shifted by `left`.
        ed2 = energy_densities(w2)
        @test maximum(abs(ed[k] - ed2[k + left]) for k in 2:19) < 1e-10

        # The spliced sites are exact vacuum: deep in the (generous) left padding the energy
        # density is a translation-invariant plateau — same-parity sites two apart agree far
        # below any packet-tail signal. (Sites 2..5 are ≥4 sites from the shifted packet.)
        @test abs(ed2[2] - ed2[4]) < 1e-5      # even sites
        @test abs(ed2[3] - ed2[5]) < 1e-5      # odd sites

        # Growing by zero is the identity (returns the state unchanged).
        @test grow_window(w; left = 0, right = 0) === w
    end

    @testset "grow_window argument checking" begin
        w  = wavepacket(qpst, 20; support = 5:16, gauge = :symmetric, window = true)
        wf = wavepacket(qpst, 20; support = 5:16, gauge = :symmetric, window = false)

        @test_throws ArgumentError grow_window(wf; left = 2)          # not a WindowMPS
        @test_throws ArgumentError grow_window(w; left = 1)           # not a multiple of Uc
        @test_throws ArgumentError grow_window(w; right = Uc + 1)     # not a multiple of Uc
        @test_throws ArgumentError grow_window(w; left = -Uc)         # negative
    end

    @testset "boundary_energy_excess and default condition" begin
        # A packet parked far from both edges has near-zero boundary excess; one parked at the
        # left edge has a large left excess and a negligible right excess.
        centered = wavepacket(qpst, 40; support = 17:24, gauge = :symmetric, window = true)
        leftnear = wavepacket(qpst, 40; support = 3:10,  gauge = :symmetric, window = true)

        cl, cr = boundary_energy_excess(centered; nsites = 6)
        @test cl < 1e-4 && cr < 1e-4

        ll, lr = boundary_energy_excess(leftnear; nsites = 6)
        @test ll > 1e-2                                             # packet touches the left edge
        @test lr < 1e-4                                            # right edge is undisturbed

        # Only a WindowMPS has explicit infinite wings to supply the vacuum reference; a bare
        # FiniteMPS window has no wings, so the boundary check is undefined there.
        wf = wavepacket(qpst, 20; support = 5:16, gauge = :symmetric, window = false)
        @test_throws ArgumentError boundary_energy_excess(wf; nsites = 4)

        # The default callback grows a side exactly when its excess clears the threshold.
        cond = window_growth_condition(nsites = 6, threshold = 1e-3, growth = Uc)
        @test cond(centered) == (0, 0)
        @test cond(leftnear) == (Uc, 0)
        # A huge threshold never triggers.
        @test window_growth_condition(threshold = 1e6)(leftnear) == (0, 0)
    end

    @testset "evolve grows the window adaptively" begin
        w = wavepacket(qpst, 20; support = 5:16, gauge = :symmetric, window = true)
        cond = window_growth_condition(nsites = 6, threshold = 1e-2, growth = Uc)
        wf, _ = evolve(w, 0.1; nsteps = 2, grow = cond)
        @test wf.psi isa WindowMPS
        @test length(wf.psi) > 20                                  # the window grew...
        @test (length(wf.psi) - 20) % Uc == 0                      # ...in whole unit cells

        # `grow` on a non-window state is an error (no infinite wings to draw vacuum from).
        wfin = wavepacket(qpst, 20; support = 5:16, gauge = :symmetric, window = false)
        @test_throws ArgumentError evolve(wfin, 0.1; nsteps = 1, grow = true)
    end
end
