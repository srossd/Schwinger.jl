# =============================================================================
# SU(F) flavor symmetry for F equal-mass flavors  (MPSKit backend, opt-in)
# =============================================================================
# When `Lattice(...; flavor_sym = true)` is set, the F flavors of a staggered site
# are gauged into a single physical site whose Hilbert space is the fermionic Fock
# space `Λ*(ℂ^F) = ⊕ₖ Λᵏ`, organized into SU(F) flavor multiplets (antisymmetric
# irreps Λᵏ, the k-fermion sectors).  The symmetry group is `U(1) × SU(F)`; the
# staggered unit cell is 2 physical sites (any F), and the hopping MPO has fixed
# bond dimension 4.  This replaces the default representation of F flavors as F
# separate U(1) sites.  See `su_nf_flavor_design.md`.
#
# The creation operator `χ† : Λᵏ → Λᵏ⁺¹` (a fundamental tensor operator, VS = fund)
# has reduced matrix element rₖ = √(k+1) for every F (canonical fermions); because
# the physical space contains only the antisymmetric irreps, χ† auto-projects onto
# the fermionic channel.  All else (mass, JW string, electric term) is SU(F)-blind.
# =============================================================================

# Antisymmetric irrep Λᵏ of SU(Nf): Young diagram = single column of k boxes.
# (Uniform SUNIrrep for all Nf ≥ 2 — it gives sector-pure quasiparticle excitations once the
# trial bond space is a proper fiducial product grid; see `_flavorsym_bond_space`.)
_flavor_irrep(Nf::Int, k::Int) = SUNIrrep{Nf}(ntuple(i -> i <= k ? 1 : 0, Nf))
# Concrete U(1) × SU(Nf) product-sector type (use `typeof(...)`, not the SUNIrrep{Nf}
# UnionAll — the latter breaks `@tensor conj(...)` with a sector-type convert error).
_flavor_sector(Nf::Int) = U1Irrep ⊠ typeof(_flavor_irrep(Nf, 1))

# Physical space of a staggered site: k fermions ↦ (charge, Λᵏ).
#   even site: charge q·k ; odd site: charge q·(k − Nf)  (particle–hole staggered shift)
function _flavorsym_physical(Nf::Int, q::Int, even::Bool)
    G = _flavor_sector(Nf)
    off = even ? 0 : Nf
    return Vect[G]((q * (k - off), _flavor_irrep(Nf, k)) => 1 for k in 0:Nf)
end

# Physical spaces of the whole chain (unit cell 2 for infinite; N sites for finite).
function _flavorsym_spaces(lattice::Lattice)
    Nf, q = lattice.F, lattice.q
    Po = _flavorsym_physical(Nf, q, false)   # site 1 is odd
    Pe = _flavorsym_physical(Nf, q, true)
    isinf(lattice.N) ? [Po, Pe] : [isodd(n) ? Po : Pe for n in 1:Int(lattice.N)]
end

# link function of a bond: the electric energy (a/2)(L + θ/2π)² of the U(1) charge.
_flavorsym_link(a, θ) = (r -> (a / 2) * (r[1].charge + θ)^2)

# Matter MPO channel-matrices (mass + hopping) and the electric link functions.
# The hopping bond is 4-dimensional (start / χ†-forward / χ-backward / done); the
# flavor index rides the fundamental auxiliary leg, so this is F-independent.
function _flavorsym_matter_and_links(lattice::Lattice; universe::Int = 0)
    Nf, q, a = lattice.F, lattice.q, lattice.a
    G  = _flavor_sector(Nf)
    Po = _flavorsym_physical(Nf, q, false)
    Pe = _flavorsym_physical(Nf, q, true)
    VS = Vect[G]((q, _flavor_irrep(Nf, 1)) => 1)          # creation carries (charge q, fund)

    nf_e(c) = c ÷ q                                       # #fermions from the U(1) charge
    nf_o(c) = c ÷ q + Nf

    mkc(P, nf) = (c = zeros(Float64, VS ⊗ P ← P);         # χ† : reduced m.e. √(k+1)
                  for (cc, b) in blocks(c); b .= sqrt(nf(cc[1].charge)); end; c)
    ce, co = mkc(Pe, nf_e), mkc(Po, nf_o)
    mkF(P, nf) = (Fop = zeros(Float64, P ← P);            # Jordan–Wigner string (−1)^#ferm
                  for (cc, b) in blocks(Fop); b .= iseven(nf(cc[1].charge)) ? 1 : -1; end; Fop)
    Fe, Fo = mkF(Pe, nf_e), mkF(Po, nf_o)
    mkM(P, nf, s) = (M = zeros(Float64, P ← P);           # staggered mass (−1)ⁿ(#ferm − Nf/2)
                     for (cc, b) in blocks(M); b .= s * (nf(cc[1].charge) - Nf / 2); end; M)
    Me = insertleftunit(insertrightunit(mkM(Pe, nf_e, +1), 2), 1)
    Mo = insertleftunit(insertrightunit(mkM(Po, nf_o, -1), 2), 1)
    mkH(cc) = (@tensor cH[-1 -2; -3] := conj(cc[-1 -3; -2]); flip(cH, 1))
    cHe, cHo = mkH(ce), mkH(co)
    cRe, cHRe = insertrightunit(ce, 3), insertrightunit(cHe, 3)
    cRo, cHRo = insertrightunit(co, 3), insertrightunit(cHo, 3)
    mkL(cc, Fop) = (@tensor t[-1; -2 -3] := cc[-3 -1; 2] * Fop[2; -2]; insertleftunit(flip(t, 3), 1))
    cLe, cHLe = mkL(ce, Fe), mkL(cHe, Fe)
    cLo, cHLo = mkL(co, Fo), mkL(cHo, Fo)

    nsites = isinf(lattice.N) ? 2 : Int(lattice.N)
    Elt = Union{Missing, typeof(cRe), scalartype(cRe)}
    A = Vector{Matrix{Elt}}(undef, nsites)
    for n in 1:nsites
        mlat = lattice.mlat[n][1]                          # DKPZ-shifted, uniform across flavors
        W = Matrix{Elt}(missing, 4, 4)
        W[1, 1] = 1.0; W[4, 4] = 1.0
        W[2, 4] = 0.5 / a * (isodd(n) ? cRo  : cRe)
        W[3, 4] = 0.5 / a * (isodd(n) ? cHRo : cHRe)
        W[1, 2] =           (isodd(n) ? cHLo : cHLe)
        W[1, 3] =           (isodd(n) ? cLo  : cLe)
        W[1, 4] = mlat *    (isodd(n) ? Mo   : Me)
        A[n] = W
    end

    θ2πu = collect(Float64, lattice.θ2π) .+ universe
    link_fcts = Vector{Union{Missing,Function}}(
        [_flavorsym_link(a, θ2πu[mod1(n, length(θ2πu))]) for n in 1:nsites])
    return A, link_fcts
end

# Full SU(F)-symmetric LEMPO Hamiltonian.  Dispatched from `MPSKitHamiltonian`.
function _flavorsym_hamiltonian(lattice::Lattice; universe::Int = 0,
                                defects::Vector{DefectCharge} = DefectCharge[])
    isempty(defects) || throw(ArgumentError("flavor_sym does not support defects (yet)."))
    all(all(isapprox.(mp, 0)) for mp in lattice.mprime) ||
        throw(ArgumentError("flavor_sym does not support the mprime (m′) term."))
    for ml in lattice.mlat                                 # SU(F) needs equal flavor masses
        all(isapprox.(ml, ml[1])) ||
            throw(ArgumentError("flavor_sym requires equal masses across flavors at each site."))
    end

    A, link_fcts = _flavorsym_matter_and_links(lattice; universe = universe)
    if isinf(lattice.N)
        lempo = InfiniteLEMPOHamiltonian(InfiniteMPOHamiltonian(A), link_fcts)
    else
        A[1]   = A[1][1:1, :]
        A[end] = A[end][:, end:end]
        lempo = FiniteLEMPOHamiltonian(FiniteMPOHamiltonian(A), link_fcts)
    end
    return MPSKitOperator(lattice, lempo, universe, defects)
end

# ---- initial-state virtual spaces (used by `loweststates`) -------------------
# Fiducial trial bond space for the infinite MPS: a PRODUCT grid of a U(1) charge cutoff
# (net fermion number c ∈ -Lmax:Lmax, charge q·c — the same cutoff the U(1)-only path uses)
# with an SU(Nf) irrep cutoff (all irreps up to 2·Lmax boxes), each sector at multiplicity
# `dim`. The two axes are linked only by N-ality: a physical site is (charge q·k, Λᵏ) with
# k ≡ n-ality(Λᵏ) mod Nf, and both add under fusion, so every bond sector obeys
# c ≡ n-ality(R) mod Nf — but each irrep is paired with EVERY allowed charge in range, not
# just the one it first fuses to. (Pairing each irrep with a single charge, as fusing a fixed
# number of physical sites does, mimics the finite lattice's zero-flux boundary and yields a
# non-injective infinite trial state, which makes the QP excitation solver leak across flavor
# sectors.) Unreachable sectors are harmlessly pruned by VUMPS/DMRG.
function _flavorsym_bond_space(lattice::Lattice, dim::Int, Lmax::Int)
    Nf, q = lattice.F, lattice.q
    G = _flavor_sector(Nf)
    Girr = typeof(_flavor_irrep(Nf, 1))
    # SU(Nf) irreps up to `repcut` boxes, tagged by n-ality (box count mod Nf = fusion level
    # mod Nf), enumerated by repeatedly fusing the fundamental.
    Vfund = Vect[Girr](_flavor_irrep(Nf, 1) => 1)
    repcut = max(2 * Lmax, Nf)
    cur = Vect[Girr](_flavor_irrep(Nf, 0) => 1)          # Λ⁰ (trivial), 0 boxes
    nality = Dict{Girr,Int}(_flavor_irrep(Nf, 0) => 0)
    for level in 1:repcut
        cur = Vect[Girr](s => 1 for s in sectors(fuse(cur ⊗ Vfund)))
        for s in sectors(cur)
            get!(nality, s, level % Nf)
        end
    end
    pairs = Pair{G,Int}[]
    for (R, na) in nality, c in -Lmax:Lmax
        mod(c - na, Nf) == 0 && push!(pairs, G(q * c, R) => dim)
    end
    return Vect[G](pairs...)
end
# Right boundary: the total-charge sector in a given flavor irrep (default: singlet Λ⁰).
_flavorsym_right_space(lattice::Lattice, total_charge::Int, R = _flavor_irrep(lattice.F, 0)) =
    Vect[_flavor_sector(lattice.F)]((total_charge, R) => 1)

# Dispatchers shared with the default (U(1)-only) path.
_mpskit_bond_space(lattice::Lattice, dim::Int, Lmax::Int) =
    lattice.flavor_sym ? _flavorsym_bond_space(lattice, dim, Lmax) :
    U1Space([lattice.q * c => dim for c in -Lmax:Lmax])
_mpskit_right_space(lattice::Lattice, total_charge::Int) =
    lattice.flavor_sym ? _flavorsym_right_space(lattice, total_charge) :
    U1Space(total_charge => 1)
# Right boundary that pins the global flavor sector to irrep `R` (flavor_sym only): used to
# target the lowest finite excitation in a chosen flavor channel by sector-selective DMRG.
_mpskit_right_space(lattice::Lattice, total_charge::Int, R) =
    _flavorsym_right_space(lattice, total_charge, R)

# ---- observables & excitations -----------------------------------------------
# Total fermion-number operator on a staggered site (eigenvalue = #fermions = k on Λᵏ).
function _flavorsym_number_op(lattice::Lattice, site::Int)
    Nf, q = lattice.F, lattice.q
    even = iseven(site)
    P = _flavorsym_physical(Nf, q, even)
    nf(c) = even ? c ÷ q : c ÷ q + Nf
    Nop = zeros(Float64, P ← P)
    for (cc, b) in blocks(Nop); b .= nf(cc[1].charge); end
    return Nop
end

# Excitation sector for a neutral quasiparticle of flavor irrep `R`: (charge 0, R).
_flavorsym_excitation_sector(lattice::Lattice, R) = U1Irrep(0) ⊠ R

"""
    flavor_singlet(lattice)

The SU(F) flavor-singlet irrep (for `loweststates(...; flavor_irrep = flavor_singlet(lattice))`).
"""
flavor_singlet(lattice::Lattice) = _flavor_irrep(lattice.F, 0)

"""
    flavor_adjoint(lattice)

The SU(F) adjoint irrep (dimension F²−1). Mesons decompose as fund ⊗ fund̄ = singlet ⊕ adjoint,
so `flavor_singlet` and `flavor_adjoint` are the two meson flavor channels.
"""
flavor_adjoint(lattice::Lattice) =
    (F = lattice.F; SUNIrrep{F}((2, ntuple(_ -> 1, F - 2)..., 0)))
