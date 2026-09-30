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

# Pad a natural-range detector profile `v` (covering sites `first_site : first_site+length(v)-1`) to a
# site-aligned length-`N` vector, filling the boundary sites the operator is undefined on with `NaN`.
# Opt-in convenience for whole-lattice sweeps/plots that want a value at every site index (see the
# `pad` keyword on `energycurrents`/`momentumdensities`).
function _pad_profile(v::AbstractVector, first_site::Int, N::Int)
    out = fill(NaN, N)
    out[first_site:first_site + length(v) - 1] .= v
    return out
end

"""
`energycurrents(state; pad = false)`

The energy current `𝒥 = T⁰¹` on each interior site (`2..N-1`), as a profile over the lattice.  See
[`chargecurrents`](@ref) for the companion charge current.

For `F = 1` (no defects) this measures `𝒥_n = -i[b_n, b_{n-1}]` as a **local 3-site contraction**
(`_contract_mpo_expval3`, O(1) memory / O(D³) per site) instead of building the N−1 full-length
per-site `FiniteMPO` operators and contracting each over the whole chain (which was O(N²)).  This
both speeds up the finite-lattice profile dramatically and lets it run on a **wavepacket window**
(`WindowMPS`/`FiniteMPS`), where the per-site `FiniteMPO` operator cannot be applied (no
`FiniteMPO * WindowMPS`).  It matches the `MPSKitEnergyCurrent` operator to machine precision.
Requires `mprime = 0` (as does the operator).  For `F > 1`/defects it falls back to the per-site
operator (finite lattice only).

`𝒥` is only defined on interior sites, so by default this returns a length-`N-2` vector (sites
`2..N-1`). Pass `pad = true` to instead get a site-aligned length-`N` vector with `NaN` at the two
boundary sites — convenient for sweeping/plotting a detector across the whole lattice without manual
range guards.
"""
function energycurrents(state::MPSKitState; pad::Bool = false)
    lat = lattice(state)
    if lat.F == 1 && isempty(state.defects) && !lat.flavor_sym && (isfinite(lat.N) || _isfinitewindow(state))
        ψ = state.psi; W = length(ψ)
        nrm2 = ψ isa MPSKit.InfiniteMPS ? 1.0 : real(dot(ψ, ψ))   # windows can drift from norm 1
        vals = [_window_energy_current(ψ, lat, n) / nrm2 for n in 2:W-1]
        return pad ? _pad_profile(vals, 2, W) : vals
    elseif _isfinitewindow(state)
        throw(ArgumentError("energycurrents on a window supports F = 1 with no defects"))
    else
        u = state.hamiltonian.universe; N = Int(lat.N)
        vals = [real(expectation(MPSKitEnergyCurrent(lat, s; universe = u), state)) for s in 2:N-1]
        return pad ? _pad_profile(vals, 2, N) : vals
    end
end

# =============================================================================
# Momentum density  p_n = T⁰¹  (matter/canonical momentum density)
#
#     p_n = (-i/4a) (χ†_n U_n U_{n+1} χ_{n+2} − h.c.)
#
# The total momentum is P = Σ_n p_n (the 1/a is already carried by each p_n).  The
# gauge-invariant bilinear χ†_n U_n U_{n+1} χ_{n+2} — a length-2 fermion hop across two
# links, on the SAME staggered sublattice as site n — is exactly the length-2 Wilson
# line `WilsonLine(n → n+2)` (Jordan–Wigner σᶻ string on site n+1 + the two gauge links
# U_n U_{n+1}).  So the operator reuses the validated `WilsonLine` on every backend, and
# `p_n = (-i/4a)(W − W†)` is manifestly Hermitian.
#
# This is DISTINCT from `EnergyCurrent` (also written T⁰¹, but the *energy* current):
# that one is the length-2 hop on the OTHER sublattice (n−1 → n+1), carries a 1/4a²
# coefficient, and adds the (m_lat) mass pieces.  `MomentumDensity` is the pure matter
# momentum density: same-sublattice hop, coefficient 1/4a, no mass term.  For a bound
# state (meson quasiparticle) Σ_n p_n is the total = centre-of-mass momentum.
#
# SIGN CONVENTION.  `WilsonLine` carries the backend-shared phase i^(finish−start) = i² =
# −1 (ED's imaginary-hopping gauge), a real overall factor that is common to W and W†.
# The physical sign of p_n is fixed by requiring a right-mover (lattice momentum k > 0) to
# have ⟨P⟩ > 0; that calibration gives the prefactor `_MOMDENS_SIGN` below (validated in
# test/momentum_density.jl against the sharp-momentum QP dispersion 𝒫(k) ≈ k/a).
# =============================================================================

const _MOMDENS_SIGN = +1   # calibrated so a k>0 wavepacket has ⟨P⟩>0 (see test/momentum_density.jl)

# W† for a length-ℓ Wilson line is the `conjugate = true` line of the same endpoints (the
# real phase i^ℓ is common to both), so p_n = (-i/4a)(W − W†) = (-i/4a)(W_false − W_true).
_momdens_prefactor(lattice::Lattice) = _MOMDENS_SIGN * (-im) / (4 * lattice.a)

"""
`EDMomentumDensity(lattice, site)`

Momentum density `p_n = T⁰¹ = (-i/4a)(χ†_n U_n U_{n+1} χ_{n+2} − h.c.)` at `site` `n`, built from
the length-2 Wilson line `WilsonLine(n → n+2)`.  The total momentum is `Σ_n` of these.  Requires a
charge-neutral sector (the underlying `EDWilsonLine` is built at `in_charge = 0`).
"""
function EDMomentumDensity(lattice::Lattice, site::Int; L_max::Union{Nothing,Int} = nothing, universe::Int = 0)
    Wf = EDWilsonLine(lattice, false, 1, site, site + 2; L_max = L_max, universe = universe)
    Wc = EDWilsonLine(lattice, true,  1, site, site + 2; L_max = L_max, universe = universe)
    return _momdens_prefactor(lattice) * (Wf + (-1.0) * Wc)
end

"""
`ITensorMomentumDensity(lattice, site)`

Momentum density `p_n = T⁰¹ = (-i/4a)(χ†_n U_n U_{n+1} χ_{n+2} − h.c.)` at `site` `n` (ITensors),
built from the length-2 Wilson line `WilsonLine(n → n+2)`.
"""
function ITensorMomentumDensity(lattice::Lattice, site::Int; L_max::Union{Nothing,Int} = nothing, universe::Int = 0)
    Wf = ITensorWilsonLine(lattice, false, 1, site, site + 2; L_max = L_max, universe = universe)
    Wc = ITensorWilsonLine(lattice, true,  1, site, site + 2; L_max = L_max, universe = universe)
    return _momdens_prefactor(lattice) * (Wf + (-1.0) * Wc)
end

"""
`MPSKitMomentumDensity(lattice, site)`

Momentum density `p_n = T⁰¹ = (-i/4a)(χ†_n U_n U_{n+1} χ_{n+2} − h.c.)` at `site` `n` (MPSKit),
built from the length-2 Wilson line `WilsonLine(n → n+2)`.  Finite lattice only — on a wavepacket
window use [`momentumdensities`](@ref), which contracts the local 3-site operator directly.

See also [`momentumdensities`](@ref), [`totalmomentum`](@ref), and [`EnergyCurrent`](@ref) (the
distinct *energy* current on the other sublattice).
"""
function MPSKitMomentumDensity(lattice::Lattice, site::Int; universe::Int = 0)
    isinf(lattice.N) && throw(ArgumentError("per-site MPSKitMomentumDensity requires a finite lattice; " *
                                            "use `momentumdensities` on a window"))
    Wf = MPSKitWilsonLine(lattice, false, 1, site, site + 2; universe = universe)
    Wc = MPSKitWilsonLine(lattice, true,  1, site, site + 2; universe = universe)
    return _momdens_prefactor(lattice) * (Wf + (-1.0) * Wc)
end

# Local 3-site expectation ⟨A1 A2 A3 | O | A1 A2 A3⟩ for standard single-physical-leg MPS
# tensors (the 3-site analogue of MPSKit.contract_mpo_expval2). Needed for a wavepacket window,
# whose infinite-lattice Hamiltonian rules out `expectation_value(ψ, FiniteMPO)` — only these
# local contractions work (see `_window_charge_current`). O[bra1 bra2 bra3; ket1 ket2 ket3].
function _contract_mpo_expval3(A1, A2, A3, O, A1b = A1, A2b = A2, A3b = A3)
    return @plansor conj(A1b[1 2; 3]) * conj(A2b[3 4; 5]) * conj(A3b[5 6; 7]) *
                    O[2 4 6; 8 9 10] * A1[1 8; 11] * A2[11 9; 12] * A3[12 10; 7]
end

# --- shared local tensor builders (F = 1) --------------------------------------------------------
# On-site mass operator (bare = false, physical mass) as a 1-site TensorMap on the physical space of
# `site`: diagonal −½ (empty) / +½ (occupied), scaled by `mlat[site]`. Identical to the closure used
# by `_energy_densities_window` and to `_mpskit_mass_fmpo`'s per-site block.
function _mpskit_massop(lat::Lattice, site::Int)
    sp = get_mpskit_spaces(lat); P = sp[mod1(site, length(sp))]
    mop = zeros(ComplexF64, P ← P)
    block(mop, U1Irrep(0)) .= -0.5
    block(mop, U1Irrep(isodd(site) ? -lat.q : lat.q)) .= 0.5
    return lat.mlat[mod1(site, length(lat.mlat))][1] * mop
end

# Symmetric hopping transport χ†_ℓ χ_{ℓ+1} + h.c. on bond ℓ as a 2-site TensorMap (the `raw`
# transport shared with `_energy_densities_window`). The kinetic hopping OPERATOR is (1/2a)·this.
function _mpskit_hop_transport(lat::Lattice, ℓ::Int)
    sp = get_mpskit_spaces(lat); q = lat.q
    Pℓ = sp[mod1(ℓ, length(sp))]; Pℓ1 = sp[mod1(ℓ + 1, length(sp))]
    raw = nothing
    for qs in (q, -q)
        openT  = ones(ComplexF64, U1Space(0 => 1)  ⊗ Pℓ  ← Pℓ  ⊗ U1Space(qs => 1))
        closeT = ones(ComplexF64, U1Space(qs => 1) ⊗ Pℓ1 ← Pℓ1 ⊗ U1Space(0 => 1))
        @tensor t[-1 -2; -3 -4] := openT[1, -1; -3, 2] * closeT[2, -2; -4, 1]
        raw = raw === nothing ? t : raw + t
    end
    return raw
end

# Energy current 𝒥_n = -i[b_n, b_{n-1}] as a local 3-site contraction on sites (n-1, n, n+1), where
# b_m = Hop(m) + ½Mass(m) + ½Mass(m+1), Hop(m) = (1/2a)·(symmetric transport). There is NO electric
# term in b_m (no lattice Poynting vector), so 𝒥 is a pure product of local operators. The two bond
# operators are built as dense 3-site operators (via ⊗ with single-site identities), commuted, and
# contracted with `_contract_mpo_expval3`. This is O(1) memory / O(D³) per site — the memory-light,
# window-capable analogue of the per-site `MPSKitEnergyCurrent` operator, which it matches to machine
# precision (validated in test/energy_currents_local.jl). Requires mprime = 0 (as does the operator).
function _window_energy_current(ψ, lat::Lattice, n::Int)
    a = lat.a
    sp = get_mpskit_spaces(lat)
    idm = TensorKit.id(sp[mod1(n - 1, length(sp))])
    id0 = TensorKit.id(sp[mod1(n,     length(sp))])
    idp = TensorKit.id(sp[mod1(n + 1, length(sp))])
    Mm = _mpskit_massop(lat, n - 1); M0 = _mpskit_massop(lat, n); Mp = _mpskit_massop(lat, n + 1)
    hopn  = (1 / (2a)) * _mpskit_hop_transport(lat, n)       # kinetic Hop on (n, n+1)
    hopnm = (1 / (2a)) * _mpskit_hop_transport(lat, n - 1)   # kinetic Hop on (n-1, n)
    bn = hopn  + 0.5 * (M0 ⊗ idp) + 0.5 * (id0 ⊗ Mp)         # b_n     on sites (n, n+1)
    bm = hopnm + 0.5 * (Mm ⊗ id0) + 0.5 * (idm ⊗ M0)         # b_{n-1} on sites (n-1, n)
    Bn = idm ⊗ bn                                            # embed b_n     on triple (n-1,n,n+1)
    Bm = bm ⊗ idp                                            # embed b_{n-1} on triple (n-1,n,n+1)
    J = (-im) * (Bn * Bm - Bm * Bn)
    return real(_contract_mpo_expval3(ψ.AC[n - 1], ψ.AR[n], ψ.AR[n + 1], J))
end

# Momentum density on window sites (n, n+1, n+2) as a *local 3-site operator tensor* contracted
# with `_contract_mpo_expval3`.  Mirrors `_window_charge_current` (open/close charge transport,
# antisymmetric in the transport direction, `i` prefactor) but spans TWO bonds, so the middle
# site carries the Jordan–Wigner σᶻ string via `_wilson_jw_passthrough` (exactly as the length-2
# `MPSKitWilsonLine` does).  Validated site-by-site against `MPSKitMomentumDensity` on a finite
# lattice (test/momentum_density.jl).
function _window_momentum_density(ψ, lat::Lattice, n::Int)
    q = lat.q
    Pn  = TensorKit.space(ψ.AC[n],   2)
    Pn1 = TensorKit.space(ψ.AC[n+1], 2)
    Pn2 = TensorKit.space(ψ.AC[n+2], 2)
    raw = nothing
    for (s, qs) in ((1.0, q), (-1.0, -q))
        openT  = ones(ComplexF64, U1Space(0 => 1)  ⊗ Pn  ← Pn  ⊗ U1Space(qs => 1))
        midT   = _wilson_jw_passthrough(Pn1, qs)          # σᶻ string + charge carry on site n+1
        closeT = ones(ComplexF64, U1Space(qs => 1) ⊗ Pn2 ← Pn2 ⊗ U1Space(0 => 1))
        @tensor t[-1 -2 -3; -4 -5 -6] := openT[1, -1; -4, 2] * midT[2, -2; -5, 3] *
                                         closeT[3, -3; -6, 1]
        raw = raw === nothing ? s * t : raw + s * t
    end
    # `WilsonLine(n → n+2)` carries the backend-shared phase i^(finish−start) = i² = −1 that puts
    # the fermion bilinear in ED's imaginary-hopping convention; this hand-built transport tensor is
    # in the raw real-hopping gauge and lacks it.  Multiply it in so the window path equals the
    # per-site `MPSKitMomentumDensity` operator exactly (verified to ~1e-10 in the finite-lattice
    # validation, where without it window = −operator).
    op = _momdens_prefactor(lat) * (im^2) * raw
    return real(_contract_mpo_expval3(ψ.AC[n], ψ.AR[n+1], ψ.AR[n+2], op))
end

"""
`momentumdensities(state::MPSKitState)`

The momentum density `p_n = T⁰¹` on each site `n = 1 … N−2` as a profile over the lattice; its sum
is the total momentum `⟨P⟩` (see [`totalmomentum`](@ref)).  For `F = 1` with no defects this uses a
local 3-site contraction (O(1) memory per site) of the length-2 charge-transport bilinear — the
memory-light analogue of [`chargecurrents`](@ref) — and so also works on a wavepacket window
(`WindowMPS`/`FiniteMPS`), where the per-site `FiniteMPO` operator cannot be applied.  It matches
`MPSKitMomentumDensity` to machine precision.  Otherwise it falls back to the per-site operator
(finite lattice only).

`p_n` is defined on sites `1..N-2`, so by default this returns a length-`N-2` vector. Pass
`pad = true` for a site-aligned length-`N` vector with `NaN` at the two trailing boundary sites
(convenient for whole-lattice sweeps/plots).
"""
function momentumdensities(state::MPSKitState; pad::Bool = false)
    lat = lattice(state)
    if lat.F == 1 && isempty(state.defects) && (isfinite(lat.N) || _isfinitewindow(state))
        ψ = state.psi
        W = length(ψ)
        nrm2 = ψ isa MPSKit.InfiniteMPS ? 1.0 : real(dot(ψ, ψ))
        vals = [_window_momentum_density(ψ, lat, n) / nrm2 for n in 1:W-2]
        return pad ? _pad_profile(vals, 1, W) : vals
    elseif _isfinitewindow(state)
        throw(ArgumentError("momentumdensities on a window supports F = 1 with no defects"))
    else
        u = state.hamiltonian.universe; N = Int(lat.N)
        vals = [real(expectation(MPSKitMomentumDensity(lat, n; universe = u), state)) for n in 1:N-2]
        return pad ? _pad_profile(vals, 1, N) : vals
    end
end

"""
`totalmomentum(state::MPSKitState)`

The total momentum `⟨P⟩ = Σ_n ⟨p_n⟩` of `state`, the sum of the [`momentumdensities`](@ref).  For a
single-quasiparticle wavepacket this is the centre-of-mass momentum; at small lattice momentum it
equals the physical momentum `p` the wavepacket was built at (up to O(k³) lattice-dispersion
corrections).

!!! note "This is the mean only"
    `⟨P⟩` is the clean, vacuum-free first moment.  Do **not** read the *distribution* over total
    momentum off the second moment `⟨(Σ_n p_n)²⟩`: summing a local density over a window imports the
    vacuum's momentum fluctuations (extensive in the window size) and `Σ_n p_n` does not commute
    with `H`, so `Var(P)` is dominated by a state-independent floor.  Obtain the distribution instead
    as the crystal-momentum distribution pushed through the dispersion `𝒫(κ)` (the diagonal of this
    operator on sharp-momentum states).

See also [`momentumdensities`](@ref), [`MomentumDensity`](@ref).
"""
totalmomentum(state::MPSKitState) = sum(momentumdensities(state))

# =============================================================================
# Generic cross-backend current / energy-current / momentum profiles (ED, ITensors)
#
# The MPSKit methods above use a memory-light local contraction (and support windows); ED and
# ITensors have no such shortcut, so these simply loop the validated per-bond / per-site operator
# (`ChargeCurrent`/`EnergyCurrent`/`MomentumDensity`) over the lattice. Their purpose is coverage
# parity — so `chargecurrents`/`energycurrents`/`momentumdensities` return a whole-lattice profile
# on *every* backend, which matters for ED↔MPSKit cross-checks. Finite lattice only.
#
# Ranges match the operators' domains of validity: charge current on bonds 1..N-1, energy current
# on interior sites 2..N-1, momentum density on sites 1..N-2 (see the operator docstrings).
# =============================================================================

"""
`chargecurrents(state::EDState)` / `chargecurrents(state::ITensorState)`

The vector current j¹ on each bond `1..N-1`, by looping the per-bond `ChargeCurrent` operator —
the ED/ITensors companion of the MPSKit [`chargecurrents`](@ref). Finite lattice only.
"""
function chargecurrents(state::EDState)
    lat = lattice(state)
    isinf(lat.N) && throw(ArgumentError("chargecurrents requires a finite lattice"))
    kw = (; L_max = state.hamiltonian.L_max, universe = state.hamiltonian.universe, charge = state.net_charge)
    return [real(expectation(EDChargeCurrent(lat, b; kw...), state)) for b in 1:Int(lat.N)-1]
end
function chargecurrents(state::ITensorState)
    lat = lattice(state)
    isinf(lat.N) && throw(ArgumentError("chargecurrents requires a finite lattice"))
    kw = (; L_max = state.hamiltonian.L_max, universe = state.hamiltonian.universe)
    return [real(expectation(ITensorChargeCurrent(lat, b; kw...), state)) for b in 1:Int(lat.N)-1]
end

"""
`energycurrents(state::EDState)` / `energycurrents(state::ITensorState)`

The energy current 𝒥 = T⁰¹ on each interior site `2..N-1`, by looping the per-site `EnergyCurrent`
operator — the ED/ITensors companion of the MPSKit [`energycurrents`](@ref). Finite lattice only.
Pass `pad = true` for a site-aligned length-`N` vector with `NaN` at the two boundary sites.
"""
function energycurrents(state::EDState; pad::Bool = false)
    lat = lattice(state)
    isinf(lat.N) && throw(ArgumentError("energycurrents requires a finite lattice"))
    N = Int(lat.N)
    kw = (; L_max = state.hamiltonian.L_max, universe = state.hamiltonian.universe, charge = state.net_charge)
    vals = [real(expectation(EDEnergyCurrent(lat, s; kw...), state)) for s in 2:N-1]
    return pad ? _pad_profile(vals, 2, N) : vals
end
function energycurrents(state::ITensorState; pad::Bool = false)
    lat = lattice(state)
    isinf(lat.N) && throw(ArgumentError("energycurrents requires a finite lattice"))
    N = Int(lat.N)
    kw = (; L_max = state.hamiltonian.L_max, universe = state.hamiltonian.universe)
    vals = [real(expectation(ITensorEnergyCurrent(lat, s; kw...), state)) for s in 2:N-1]
    return pad ? _pad_profile(vals, 2, N) : vals
end

"""
`momentumdensities(state::EDState)` / `momentumdensities(state::ITensorState)`

The momentum density p_n = T⁰¹ on each site `1..N-2`, by looping the per-site `MomentumDensity`
operator — the ED/ITensors companion of the MPSKit [`momentumdensities`](@ref). Its sum is the
total momentum ⟨P⟩ (see [`totalmomentum`](@ref)). Finite lattice only; a charge-neutral sector is
required (the underlying Wilson line is built at `in_charge = 0`). Pass `pad = true` for a
site-aligned length-`N` vector with `NaN` at the two trailing boundary sites.
"""
function momentumdensities(state::EDState; pad::Bool = false)
    lat = lattice(state)
    isinf(lat.N) && throw(ArgumentError("momentumdensities requires a finite lattice"))
    N = Int(lat.N)
    kw = (; L_max = state.hamiltonian.L_max, universe = state.hamiltonian.universe)
    vals = [real(expectation(EDMomentumDensity(lat, n; kw...), state)) for n in 1:N-2]
    return pad ? _pad_profile(vals, 1, N) : vals
end
function momentumdensities(state::ITensorState; pad::Bool = false)
    lat = lattice(state)
    isinf(lat.N) && throw(ArgumentError("momentumdensities requires a finite lattice"))
    N = Int(lat.N)
    kw = (; L_max = state.hamiltonian.L_max, universe = state.hamiltonian.universe)
    vals = [real(expectation(ITensorMomentumDensity(lat, n; kw...), state)) for n in 1:N-2]
    return pad ? _pad_profile(vals, 1, N) : vals
end

"""
`totalmomentum(state::EDState)` / `totalmomentum(state::ITensorState)`

The total momentum ⟨P⟩ = Σ_n ⟨p_n⟩, the sum of the [`momentumdensities`](@ref). Finite lattice only.
"""
totalmomentum(state::Union{EDState,ITensorState}) = sum(momentumdensities(state))
