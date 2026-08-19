# =============================================================================
# Conserved lattice currents: the vector (charge) current j¹ and the energy
# current T⁰¹ = 𝒥.  Both use the SAME sign convention: the current through a
# point is  -i[H, (conserved charge in the region to the LEFT of that point)],
# which is the flux flowing to the RIGHT (+x), positive.  Both are gauge-
# independent (built from each backend's own hopping/mass operators) and exactly
# conserved on the lattice by construction.
#
# The two densities live on dual sublattices of the staggered chain:
#   • charge density ρ_n on SITES, charge current j¹_n on BONDS;
#   • energy density h_m on BONDS,  energy current 𝒥_n on SITES.
#
# Charge current through bond (n, n+1) — flux rightward past the bond:
#     j¹_n = -i [H, Q_{≤n}] = -i q [Hop(n), N_n],
# from ∂_t ρ_n = -(j¹_n - j¹_{n-1}) (charge in from the left bond minus out the
# right bond).  Q_{≤n} = q Σ_{m≤n} N_m is the charge left of the bond; only the
# bond-n hopping fails to commute with it.  N_n = Σ_f χ†_n χ_n is the site number
# operator (here (-1)^n·Mass(n; bare); the additive constant drops in the
# commutator).  Requires mprime = 0.
#
# Energy current through site n — flux rightward past the site:
#     𝒥_n = -i [H, E_{<n}] = -i [b_n, b_{n-1}],   b_m = Hop(m) + ½[Mass(m)+Mass(m+1)],
# from ∂_t h_n = 𝒥_n - 𝒥_{n+1} (energy in from the left site minus out the right
# site) — the exact analogue of the charge relation.  E_{<n} = Σ_{m<n} h_m is the
# energy left of the site.  The electric part of h_m contributes nothing to the
# commutator (no lattice Poynting vector), so 𝒥 is a plain product of local
# operators, reproducing −(i/4a²)(χ†_{n-1}U_{n-1}U_n χ_{n+1} − h.c.) +
# (m_lat(-1)^n/2)(j¹_n + j¹_{n-1}).
# =============================================================================

_comm(A::SchwingerOperator, B::SchwingerOperator) = A * B + (-1.0) * (B * A)

# MPSKitHopping(lattice, bond) returns a vector of per-flavor/direction operators;
# ED/ITensor return a single operator. Collapse to a single operator either way.
_assingle(op::SchwingerOperator) = op
_assingle(ops::AbstractVector) = sum(ops)

# MPSKit cannot multiply/add across MPO types: `MPSKitHopping` is a plain `FiniteMPO{TensorMap}`
# while `MPSKitMass` is a `FiniteMPOHamiltonian{JordanMPOTensorMap}`, and mixing the two (or the
# converted `FiniteMPO{JordanMPOTensorMap}`) fails in the space fuser. So we rebuild the per-site
# staggered mass / number operator directly as a plain `FiniteMPO{TensorMap}` — one trivial-bond
# MPO per flavor (identity elsewhere), summed over flavors — matching the hopping's element type so
# every `*`/`+` in the commutators stays within one MPO type. Mirrors `MPSKitWilsonLine`'s
# zero-length number-operator idiom and `MPSKitMass`'s diagonal ±0.5 blocks.
function _mpskit_mass_fmpo(lattice::Lattice, site::Int; universe::Int = 0, bare::Bool = true)
    isinf(lattice.N) && throw(ArgumentError("per-site mass FiniteMPO requires a finite lattice"))
    N, F = Int(lattice.N), lattice.F
    _, universe = process_L_max_universe(lattice, nothing, universe)
    sp = collect(get_mpskit_spaces(lattice))
    idop(P) = isomorphism(ComplexF64, U1Space(0 => 1) ⊗ P, P ⊗ U1Space(0 => 1))
    function massop(P, coeff)                              # diagonal: empty −½, occupied +½
        T = zeros(ComplexF64, U1Space(0 => 1) ⊗ P ← P ⊗ U1Space(0 => 1))
        block(T, U1Irrep(0)) .= coeff * -0.5
        block(T, U1Irrep(isodd(site) ? -lattice.q : lattice.q)) .= coeff * 0.5
        return T
    end
    MPOT = typeof(idop(sp[1]))
    fmpo = nothing
    for f in 1:F
        idx = (site - 1)*F + f
        coeff = ComplexF64(bare ? 1 : lattice.mlat[site][f])
        ts = MPOT[k == idx ? massop(sp[k], coeff) : idop(sp[k]) for k in 1:N*F]
        term = MPSKit.FiniteMPO(ts)
        fmpo = fmpo === nothing ? term : fmpo + term
    end
    return MPSKitOperator(lattice, fmpo, universe)
end

# ---- ED ----
"""
`EDChargeCurrent(lattice, bond)`

Charge (vector) current `j¹` through `bond` `(n,n+1)`: the rightward charge flux
`j¹_n = -i[H, Q_{≤n}] = -i q [Hop(n), N_n]`. Requires `mprime = 0`.
"""
function EDChargeCurrent(lattice::Lattice, bond::Int; L_max::Union{Nothing,Int} = nothing,
                         universe::Int = 0, charge::Int = 0)
    Hn = _assingle(EDHopping(lattice, bond; L_max = L_max, universe = universe, bare = false, charge = charge))
    Nn = EDMass(lattice, bond;    L_max = L_max, universe = universe, bare = true,  charge = charge)
    return (-im * lattice.q * (-1)^bond) * _comm(Hn, Nn)
end

"""
`EDEnergyCurrent(lattice, site)`

Energy current `𝒥 = T⁰¹` through `site` `n`: the rightward energy flux
`𝒥_n = -i[H, E_{<n}] = -i [b_n, b_{n-1}]`, `b_m = Hop(m) + ½[Mass(m)+Mass(m+1)]`. Requires `mprime = 0`.
"""
function EDEnergyCurrent(lattice::Lattice, site::Int; L_max::Union{Nothing,Int} = nothing,
                         universe::Int = 0, charge::Int = 0)
    b(m) = _assingle(EDHopping(lattice, m; L_max = L_max, universe = universe, bare = false, charge = charge)) +
           0.5 * EDMass(lattice, m;     L_max = L_max, universe = universe, bare = false, charge = charge) +
           0.5 * EDMass(lattice, m + 1; L_max = L_max, universe = universe, bare = false, charge = charge)
    return -im * _comm(b(site), b(site - 1))
end

# ---- ITensors ----
"""
`ITensorChargeCurrent(lattice, bond)`

Charge (vector) current `j¹` through `bond` `(n,n+1)`: the rightward charge flux
`j¹_n = -i[H, Q_{≤n}] = -i q [Hop(n), N_n]`. Requires `mprime = 0`.
"""
function ITensorChargeCurrent(lattice::Lattice, bond::Int; L_max::Union{Nothing,Int} = nothing, universe::Int = 0)
    Hn = _assingle(ITensorHopping(lattice, bond; L_max = L_max, universe = universe, bare = false))
    Nn = ITensorMass(lattice, bond;    L_max = L_max, universe = universe, bare = true)
    return (-im * lattice.q * (-1)^bond) * _comm(Hn, Nn)
end

"""
`ITensorEnergyCurrent(lattice, site)`

Energy current `𝒥 = T⁰¹` through `site` `n`: the rightward energy flux
`𝒥_n = -i[H, E_{<n}] = -i [b_n, b_{n-1}]`, `b_m = Hop(m) + ½[Mass(m)+Mass(m+1)]`. Requires `mprime = 0`.
"""
function ITensorEnergyCurrent(lattice::Lattice, site::Int; L_max::Union{Nothing,Int} = nothing, universe::Int = 0)
    b(m) = _assingle(ITensorHopping(lattice, m; L_max = L_max, universe = universe, bare = false)) +
           0.5 * ITensorMass(lattice, m;     L_max = L_max, universe = universe, bare = false) +
           0.5 * ITensorMass(lattice, m + 1; L_max = L_max, universe = universe, bare = false)
    return -im * _comm(b(site), b(site - 1))
end

# ---- MPSKit ----
"""
`MPSKitChargeCurrent(lattice, bond)`

Charge (vector) current `j¹` through `bond` `(n,n+1)`: the rightward charge flux
`j¹_n = -i[H, Q_{≤n}] = -i q [Hop(n), N_n]`. Requires `mprime = 0`. (See `chargecurrents` for a
memory-light profile that also works on a wavepacket window.)
"""
function MPSKitChargeCurrent(lattice::Lattice, bond::Int; universe::Int = 0)
    Hn = _assingle(MPSKitHopping(lattice, bond; universe = universe, bare = false))
    Nn = _mpskit_mass_fmpo(lattice, bond;       universe = universe, bare = true)
    return (-im * lattice.q * (-1)^bond) * _comm(Hn, Nn)
end

"""
`MPSKitEnergyCurrent(lattice, site)`

Energy current `𝒥 = T⁰¹` through `site` `n`: the rightward energy flux
`𝒥_n = -i[H, E_{<n}] = -i [b_n, b_{n-1}]`, `b_m = Hop(m) + ½[Mass(m)+Mass(m+1)]`. Requires `mprime = 0`.
"""
function MPSKitEnergyCurrent(lattice::Lattice, site::Int; universe::Int = 0)
    b(m) = _assingle(MPSKitHopping(lattice, m; universe = universe, bare = false)) +
           0.5 * _mpskit_mass_fmpo(lattice, m;     universe = universe, bare = false) +
           0.5 * _mpskit_mass_fmpo(lattice, m + 1; universe = universe, bare = false)
    return -im * _comm(b(site), b(site - 1))
end

# Charge current on window bond ℓ (between MPS sites ℓ, ℓ+1), as a *local 2-site operator tensor*
# contracted with `contract_mpo_expval2`.  This is required for a `WindowMPS`/`FiniteMPS` window:
# its Hamiltonian is the infinite lattice, so `expectation_value(ψ, FiniteMPO)` fails (there is no
# `FiniteMPO * WindowMPS`) — only local `contract_mpo_expval1/2` contractions work, exactly as
# `_energy_densities_window` does for the hopping term.  The tensor
#     j¹_ℓ = (q/2a)·i · Σ_{qs=±q} sgn(qs)·(χ transport qs across the bond)
# was validated bond-by-bond against the per-bond `MPSKitChargeCurrent` operator on a finite
# θ-quench state (agrees to machine precision; the transport tensor already carries the (−1)^ℓ
# staggering, so no extra sign is needed).  `openT`/`closeT` mirror `_energy_densities_window`'s
# `hopop`, but combined antisymmetrically (sgn(qs)) with an `i` prefactor to make the current.
function _window_charge_current(ψ, q::Int, a::Float64, ℓ::Int)
    Pℓ  = TensorKit.space(ψ.AC[ℓ],   2)                   # physical space at site ℓ
    Pℓ1 = TensorKit.space(ψ.AC[ℓ+1], 2)                   # physical space at site ℓ+1
    raw = nothing
    for (s, qs) in ((1.0, q), (-1.0, -q))
        openT  = ones(ComplexF64, U1Space(0 => 1)  ⊗ Pℓ  ← Pℓ  ⊗ U1Space(qs => 1))
        closeT = ones(ComplexF64, U1Space(qs => 1) ⊗ Pℓ1 ← Pℓ1 ⊗ U1Space(0 => 1))
        @tensor t[-1 -2; -3 -4] := openT[1, -1; -3, 2] * closeT[2, -2; -4, 1]
        raw = raw === nothing ? s * t : raw + s * t
    end
    op = (q / (2a)) * im * raw
    return real(MPSKit.contract_mpo_expval2(ψ.AC[ℓ], ψ.AR[ℓ+1], op))
end

# Bare hopping-mass (pseudoscalar) bilinear on bond ℓ, as a *local 2-site operator tensor* contracted
# with `contract_mpo_expval2` — the memory-light companion of `_window_charge_current`, mirroring
# `_energy_densities_window`'s `hopop`.  This is the SAME symmetric transport χ†_ℓ χ_{ℓ+1} + h.c.
# that the charge current uses, but combined SYMMETRICALLY (no sgn(qs), no i) and with the (−1)^{ℓ+1}
# staggering that `matrices_hoppingmass` gives the bare mprime term (odd site_ind → +, even → −).
# Returns the bare `MPSKitHoppingMass(lat, ℓ)` bond expectation p(ℓ); `pseudoscalardensity` combines
# neighbouring bonds as (1/a)(p(site)/2 + p(before)/2).  Validated bond-by-bond against the per-site
# operator path (`pseudoscalardensity`) on a finite θ-quench state (machine precision).
function _window_pseudoscalar_bond(ψ, q::Int, ℓ::Int)
    Pℓ  = TensorKit.space(ψ.AC[ℓ],   2)
    Pℓ1 = TensorKit.space(ψ.AC[ℓ+1], 2)
    raw = nothing
    for qs in (q, -q)
        openT  = ones(ComplexF64, U1Space(0 => 1)  ⊗ Pℓ  ← Pℓ  ⊗ U1Space(qs => 1))
        closeT = ones(ComplexF64, U1Space(qs => 1) ⊗ Pℓ1 ← Pℓ1 ⊗ U1Space(0 => 1))
        @tensor t[-1 -2; -3 -4] := openT[1, -1; -3, 2] * closeT[2, -2; -4, 1]
        raw = raw === nothing ? t : raw + t
    end
    op = ((-1)^(ℓ + 1)) * raw
    return real(MPSKit.contract_mpo_expval2(ψ.AC[ℓ], ψ.AR[ℓ+1], op))
end

# Uniform bulk pseudoscalar density of an infinite (translation-invariant) vacuum, e.g. a wavepacket
# wing (`left_gs`/`right_gs`).  Every site is flanked by one odd and one even bond, so the density is
# uniform = (p_odd + p_even)/(2a); averaging both sublattice bonds makes this robust to the wing's
# absolute parity offset.  Used to replace the outermost window sites, which miss the bond into the
# wing (exactly as `_energy_densities_window` replaces them with the wing vacuum energy density).
function _vacuum_pseudoscalar_density(vac::MPSKitState, q::Int, a::Float64)
    ψ = vac.psi                                          # an InfiniteMPS wing vacuum
    return (_window_pseudoscalar_bond(ψ, q, 1) + _window_pseudoscalar_bond(ψ, q, 2)) / (2a)
end

"""
`pseudoscalardensities(state::MPSKitState)`

The pseudoscalar density P = ⟨ψ̄ iγ⁵ψ⟩ on each site as a profile over the lattice. For `F = 1` with no
defects this uses a local 2-site contraction (`contract_mpo_expval2`, O(1) memory per bond) of the
bare hopping-mass bilinear — the memory-light analogue of `chargecurrents`, avoiding the O(N²) memory
of materialising the N per-site `MPSKitHoppingMass` operators. It matches the per-site
`pseudoscalardensity` to machine precision, and works on a wavepacket window (`WindowMPS`) as well as a
finite lattice.  Otherwise it falls back to the generic per-site path.

On a `WindowMPS`, the two outermost sites miss the bond into the (infinite-vacuum) wing, which would
leave them at half the true density; they are replaced with the wing's uniform vacuum pseudoscalar
density so the profile connects smoothly to the background — mirroring `_energy_densities_window`.

This is the operator that appears in the explicit-mass term of the axial Ward identity,
∂_μ j₅^μ = (q/π)E + q·m_lat·P (lattice; continuum reading 2m·⟨ψ̄iγ⁵ψ⟩).
"""
function pseudoscalardensities(state::MPSKitState)
    lat = lattice(state)
    if lat.F == 1 && isempty(state.defects) && (isfinite(lat.N) || _isfinitewindow(state))
        ψ = state.psi; W = length(ψ)
        nrm2 = ψ isa MPSKit.InfiniteMPS ? 1.0 : real(dot(ψ, ψ))   # windows can drift from norm 1
        p = [_window_pseudoscalar_bond(ψ, lat.q, ℓ) / nrm2 for ℓ in 1:W-1]   # bare bond bilinears
        # pseudoscalardensity(site) = (1/a)·½(left bond + right bond); a boundary site keeps only its
        # one existing bond (open end), so it is half-weight — then repaired to the wing value below.
        pds = [_averaged_bond(ℓ -> p[ℓ], site, W, false) / lat.a for site in 1:W]
        if ψ isa WindowMPS                                        # connect boundaries to the wings
            lv, rv = _window_vacua(state)
            pds[1]   = _vacuum_pseudoscalar_density(lv, lat.q, lat.a)
            pds[end] = _vacuum_pseudoscalar_density(rv, lat.q, lat.a)
        end
        return pds
    else
        N = isinf(lat.N) ? 2 : Int(lat.N)
        return pseudoscalardensity.(Ref(state), 1:N)
    end
end

"""
`chargecurrents(state)` / `energycurrents(state)`

The vector current j¹ on each bond (`1..N-1`) and the energy current 𝒥 on each interior site
(`2..N-1`), as a profile over the lattice — observable convenience wrappers around the
`ChargeCurrent`/`EnergyCurrent` operators.  For `F = 1` (no defects) `chargecurrents` measures the
current as a local 2-site operator tensor via `contract_mpo_expval2` (O(1) memory per bond, no
full-length MPO built): this both avoids materialising the N−1 per-bond `FiniteMPO` operators for a
2-site observable (O(N²) memory) and lets it run on a **wavepacket window** (`WindowMPS`/`FiniteMPS`
cut from an infinite background, e.g. an evolving soliton), where the per-bond `FiniteMPO` operator
cannot be applied (`expectation_value(ψ, mpo)` would need a `FiniteMPO * WindowMPS` product MPSKit
does not define).  It matches the `MPSKitChargeCurrent` operator to machine precision.  For
`F > 1`/defects it falls back to the per-bond operator (finite lattice only).  `energycurrents` on a
window is not yet supported (its `𝒥_n = -i[b_n, b_{n-1}]` is a 3-site operator).
"""
function chargecurrents(state::MPSKitState)
    lat = lattice(state)
    # Measure the current as a LOCAL 2-site operator (`contract_mpo_expval2`) whenever possible: this
    # works on any MPS exposing `AC`/`AR` (a finite `FiniteMPS` or a wavepacket `WindowMPS`), costs
    # O(1) memory per bond, and avoids materialising N−1 full-length `FiniteMPO` operators for a
    # 2-site observable (which was O(N²) memory). It matches the `MPSKitChargeCurrent` operator to
    # machine precision. Restricted to F = 1 with no defects (MPS bond = lattice bond); otherwise fall
    # back to the per-bond operator (finite lattice only — a window has no such operator).
    if lat.F == 1 && isempty(state.defects)
        ψ = state.psi; nrm2 = real(dot(ψ, ψ))
        return [_window_charge_current(ψ, lat.q, lat.a, ℓ) / nrm2 for ℓ in 1:length(ψ)-1]
    elseif _isfinitewindow(state)
        throw(ArgumentError("chargecurrents on a window supports F = 1 with no defects"))
    else
        u = state.hamiltonian.universe
        return [real(expectation(MPSKitChargeCurrent(lat, b; universe = u), state)) for b in 1:Int(lat.N)-1]
    end
end

"""
`energycurrents(state)`

The energy current `𝒥 = T⁰¹` on each interior site (`2..N-1`), as a profile over the lattice, via
the per-site `EnergyCurrent` operator.  See [`chargecurrents`](@ref) for the companion charge
current.  Not supported on a wavepacket window (`𝒥_n = -i[b_n, b_{n-1}]` is a 3-site operator);
measure on a finite lattice, or use `energy_densities` on the window.
"""
function energycurrents(state::MPSKitState)
    lat = lattice(state)
    if _isfinitewindow(state)
        # 𝒥_n = -i[b_n, b_{n-1}] spans sites n-1,n,n+1 — a genuine 3-site operator that the local
        # `contract_mpo_expval2` machinery (all that works on a WindowMPS) cannot contract, and the
        # `FiniteMPO` operator cannot be applied to a window (no `FiniteMPO * WindowMPS`). Not yet
        # supported on windows; measure the energy current on a finite lattice instead.
        throw(ArgumentError("energycurrents on a window is not yet supported (3-site operator); " *
                            "use a finite lattice, or energy_densities for the window"))
    else
        u = state.hamiltonian.universe
        return [real(expectation(MPSKitEnergyCurrent(lat, s; universe = u), state)) for s in 2:Int(lat.N)-1]
    end
end
