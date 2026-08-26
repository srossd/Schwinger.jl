using Test
using Schwinger
using MPSKit, MPSKitLEMPO, TensorKit, LinearAlgebra

# ⟨vac|𝒯⁰⁰₀|vac⟩ — the cell-0 vacuum expectation of the FiniteLEMPO operator, computed
# with the same machinery `qp_matrix_element` uses for its connected subtraction constant.
_vacuum_cell(O, vac) =
    MPSKitLEMPO._vacuum_expectation(O, MPSKitLEMPO.QPMEContext(vac.psi, length(vac.psi)))

# The `emt` cell operator's translates by the two-site unit cell sum to H, so its expectation
# values must agree with `energy_density`: on the vacuum ⟨𝒯⁰⁰₀⟩ is the two-site-cell energy
# 2a·energy_density(vac), and on a sharp-momentum quasiparticle the forward diagonal is E(p).
# The last case has mprime ≠ 0 (staggered hopping-mass), folded into `emt`'s hopping term.
@testset "emt vacuum cell energy == 2a·energy_density" begin
    for (a, m, θ, mp) in ((1.0, 0.6, 0.0, 0.0), (0.5, 0.5, 0.3, 0.0), (0.7, 1.0, -0.2, 0.4))
        lat = Lattice(Inf; F = 1, a = a, m = m, θ2π = θ, mprime = mp)
        Hop = Hamiltonian(lat; backend = :MPSKit)
        vac = groundstate(Hop; bonddim = 40)

        cell_from_ed = real(energy_density(vac)) * 2 * lat.a
        cell_total   = real(_vacuum_cell(emt(lat),                vac))
        cell_matter  = real(_vacuum_cell(emt(lat; part = :matter), vac))
        cell_mass    = real(_vacuum_cell(emt(lat; part = :mass),   vac))
        cell_hop     = real(_vacuum_cell(emt(lat; part = :hopping), vac))

        @test cell_total ≈ cell_from_ed rtol = 1e-6
        # mass + hopping + electric split stays exactly additive
        @test cell_mass + cell_hop ≈ cell_matter rtol = 1e-8
    end
end

@testset "emt forward oracle == E(p)" begin
    lat = Lattice(Inf; F = 1, a = 0.6, m = 0.7, θ2π = 0.0, mprime = 0.3)
    Hop = Hamiltonian(lat; backend = :MPSKit)
    for pmom in (0.0, 0.4)
        qp = loweststates(Hop, 2; momentum = pmom, bonddim = 40)[2]   # [1]=vacuum, [2]=lightest QP
        Ep = real(energy(qp))
        me = real(qp_matrix_element(qp.psi, qp.psi, emt(lat))) / real(dot(qp.psi, qp.psi))
        @test me ≈ Ep rtol = 1e-4
    end
end
