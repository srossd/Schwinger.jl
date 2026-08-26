# =============================================================================
# Cell-0 energy density 𝒯⁰⁰₀ of the Schwinger model as a `FiniteLEMPOHamiltonian`.
#
# The k-local matter terms (mass, hopping) live in the `FiniteMPOHamiltonian` part;
# the electric energy (a/2)(L+θ)², diagonal in the accumulated bond charge, lives in
# the `link_fcts`.  The operator is supported on sites 1..3 (the hop on bond 2 reaches
# site 3) and its translates by the two-site unit cell sum to H, so it is exactly
# normalized (Σ_cells 𝒯⁰⁰ = H ⇒ Z₀₀ = 1, i.e. G(0) = 1).  Use with
# `MPSKitLEMPO.qp_matrix_element` on infinite sharp-momentum quasiparticle states.
# =============================================================================

"""
    emt(lat::Lattice; part=:total)

Cell-0 energy density `𝒯⁰⁰₀` of the Schwinger model as a `FiniteLEMPOHamiltonian`, for
exact quasiparticle EMT matrix elements `⟨ϕ(p′)|𝒯⁰⁰₀|ϕ(p)⟩` (via
`MPSKitLEMPO.qp_matrix_element`). In 1+1d the elastic EMT matrix element is rank-1
(`𝔐^μν = 2P^μP^ν G(t)`), so `T⁰⁰` alone suffices and is exactly normalized
(`Σ_cells 𝒯⁰⁰ = H ⇒ Z₀₀ = 1`, i.e. `G(0)=1` with no renormalization constant):

    G(t) = |qp_matrix_element(bra.qp, ket.qp, emt(lat); normalize=true)| / E(p),  t = -4p².

Cell 0 owns sites `{1,2}`, links `{1,2}`, bond `{1,2}`; support is sites `1..3` (the hop
on bond 2 reaches site 3). `𝒯⁰⁰₀ = mass(1,2) + electric(1,2) + hopping(1,2)`.

`part`:
- `:total`   — mass + hopping + electric (default; `Σ_cells = H`).
- `:matter`  — mass + hopping (the fermion piece; pure kinetic at `m=0`, `mprime=0`).
- `:mass`    — staggered mass only (`m_lat (-1)^s χ†χ`; identically 0 at `m=0`).
- `:hopping` — the bond hopping term `(1/2a + (-1)^(ℓ+1) mprime)(χ†χ + h.c.)`: the kinetic
  hop plus the staggered hopping-mass `mprime` (pure kinetic at `mprime=0`).
Electric is `total − matter`, so the mass/hopping/electric split stays exactly additive.

The vacuum expectation `⟨vac|𝒯⁰⁰₀|vac⟩` equals the two-site-cell energy `2a · energy_density(vac)`
and the forward diagonal `qp_matrix_element(ϕ,ϕ,emt(lat))/⟨ϕ|ϕ⟩` equals `E(p)`, matching
[`energy_densities`](@ref)/[`energy_density`](@ref) (including the `mprime` hopping-mass).
"""
function emt(lat::Lattice; part::Symbol = :total)
    part in (:total, :matter, :mass, :hopping) ||
        throw(ArgumentError("part must be :total, :matter, :mass, or :hopping"))
    lat.F == 1 || throw(ArgumentError("emt currently supports F = 1"))
    q = lat.q
    even = U1Space(0 => 1, q => 1)
    odd  = U1Space(-q => 1, 0 => 1)
    P(s) = isodd(s) ? odd : even

    function massop(s)                                   # m_lat (-1)^s χ†χ
        ml = lat.mlat[mod1(s, length(lat.mlat))][1]
        O = TensorKit.id(ComplexF64, P(s))
        for (c, b) in blocks(O)
            b .*= ml * (c.charge == 0 ? -0.5 : 0.5)
        end
        return O
    end

    function hopop(s0)                                   # coeff Σ_{±q} χ†_{s0} χ_{s0+1} + h.c.
        # kinetic hopping (1/2a) plus the staggered hopping-mass (-1)^(s0+1) m′ on bond s0:
        # both are the same χ†χ + h.c. transport, so they share one coefficient (as in
        # `_energy_densities_window`); m′ ≡ 0 leaves the pure kinetic term.
        coeff = 1 / (2 * lat.a) +
                (-1)^(s0 + 1) * lat.mprime[mod1(s0, length(lat.mprime))][1]
        O = nothing
        for qs in (q, -q)
            OL = coeff * ones(ComplexF64, P(s0) ← P(s0) ⊗ U1Space(qs => 1))
            OR = ones(ComplexF64, U1Space(qs => 1) ⊗ P(s0 + 1) ← P(s0 + 1))
            @tensor o[-1 -2; -3 -4] := OL[-1; -3 1] * OR[1 -2; -4]
            O = O === nothing ? o : O + o
        end
        return O
    end

    # electric energy on link ℓ: (a/2)(L_ℓ + θ/2π)², diagonal in the accumulated bond charge
    Fel(ℓ) = (c::U1Irrep -> (lat.a / 2) *
              (c.charge + Float64(lat.θ2π[mod1(ℓ, length(lat.θ2π))]))^2)

    massterms = [(1,) => massop(1), (2,) => massop(2)]
    hopterms  = [(1, 2) => hopop(1), (2, 3) => hopop(2)]
    terms = part === :mass    ? massterms :
            part === :hopping ? hopterms  :
                                vcat(massterms, hopterms)          # :total, :matter
    mpo = FiniteMPOHamiltonian([P(1), P(2), P(3)], terms)
    links = part === :total ? Union{Missing,Function}[Fel(1), Fel(2), missing] :
                              Union{Missing,Function}[missing, missing, missing]
    return FiniteLEMPOHamiltonian(mpo, links)
end
