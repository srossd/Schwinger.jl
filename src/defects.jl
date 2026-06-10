# =============================================================================
# Static defect charges
# =============================================================================

"""
`DefectCharge(link, charge)`

A static (external) charge of strength `charge` inserted into the lattice on the
left side of link `link`, i.e. between physical sites `link-1` and `link`.

Its only effect is on Gauss's law: it shifts the electric flux on every link
`n ≥ link` by `+charge`.  Equivalently, it adds `charge` to `θ2π[n]` for `n ≥ link`.

A `Hamiltonian` may be built with a list of `DefectCharge`s (keyword `defects`):
  * ED / ITensors backends realise them as that θ2π shift (no extra site);
  * the MPSKit backend inserts a genuine extra lattice site carrying the charge
    (a 1-dimensional `U1Space(charge => 1)`), inert under the matter Hamiltonian
    but contributing to the gauge/Gauss-law accumulation. The inserted site (and
    its zero-width link) are invisible to Schwinger.jl observables.
"""
struct DefectCharge
    link::Int
    charge::Int
end

# defects mapped to a per-link θ2π shift (length N); used by ED / ITensors
function _defect_theta_shift(lattice::Lattice, defects)
    N = Int(lattice.N)
    Δθ = zeros(Float64, N)
    for d in defects
        (2 ≤ d.link ≤ N) || throw(ArgumentError("defect link $(d.link) must lie in 2:$N"))
        for n in d.link:N
            Δθ[n] += d.charge
        end
    end
    return Δθ
end

# a copy of `lattice` with θ2π shifted by the defects (preserving everything else)
function _lattice_with_defects(lattice::Lattice, defects)
    N = Int(lattice.N)
    θ = collect(Float64, lattice.θ2π[1:N]) .+ _defect_theta_shift(lattice, defects)
    mlatmat   = permutedims(reduce(hcat, lattice.mlat))     # N×F
    mprimemat = permutedims(reduce(hcat, lattice.mprime))   # N×F
    return Lattice(N; F=lattice.F, q=lattice.q, periodic=lattice.periodic,
                   a=lattice.a, mlat=mlatmat, mprime=mprimemat, θ2π=θ)
end

# MPSKit chain index of matter (physical) site `j`, flavor `k`, accounting for
# the defect sites that have been inserted to its left.  This makes the inserted
# defect sites invisible to Schwinger.jl observables.
function _matter_index(lattice::Lattice, defects, j::Int, k::Int = 1)
    base = (j - 1) * lattice.F + k
    return base + count(d -> (d.link - 1) * lattice.F + 1 <= base, defects)
end

# MPSKit chain index of the defect-site for `defect` within `defects` (defects are
# inserted in ascending link order, so earlier ones shift later ones rightward).
function _defect_position(lattice::Lattice, defects, defect::DefectCharge)
    sorted = sort(collect(defects); by = x -> x.link)
    j = findfirst(==(defect), sorted)
    isnothing(j) && throw(ArgumentError("defect $defect not present"))
    return (defect.link - 1) * lattice.F + 1 + (j - 1)
end

# `defects` with `defect` replaced (relocated) by `repl`, or removed if `repl===nothing`
function _edit_defect(defects, defect::DefectCharge, repl::Union{DefectCharge,Nothing})
    i = findfirst(==(defect), defects)
    isnothing(i) && throw(ArgumentError("defect $defect not present"))
    out = copy(collect(defects))
    isnothing(repl) ? deleteat!(out, i) : (out[i] = repl)
    return out
end

# the FiniteMPS site tensors as a plain vector (center gauged to site 1)
_mps_tensors(psi::MPSKit.FiniteMPS) = [i == 1 ? psi.AC[1] : psi.AR[i] for i in 1:length(psi)]

# ---- moving a (1-D, entanglement-free) defect site within a FiniteMPS ----------
# Swap MPS sites `a` and `a+1` exactly: contract the pair, swap the two physical
# legs, and re-split by QR. Exact (no truncation) because one of the sites carries
# no entanglement; for non-abelian symmetries the leg swap is the braiding/F-move.
function _swap_mps_sites!(tensors, a::Int)
    L = tensors[a]; R = tensors[a + 1]
    @tensor T[vl, pL, pR; vr] := L[vl, pL; m] * R[m, pR; vr]
    # regroup so the RIGHT site's leg joins the left bond, then exact QR-split
    T2 = TensorKit.permute(T, ((1, 3), (2, 4)))            # (vl ⊗ pR) ← (pL ⊗ vr)
    Q, C = TensorKit.left_orth(T2)                          # Q:(vl⊗pR)←s, C: s←(pL⊗vr)
    tensors[a]     = Q                                      # now carries the right site's leg
    tensors[a + 1] = TensorKit.permute(C, ((1, 2), (3,)))   # (s ⊗ pL) ← vr : the left site's leg
    return tensors
end

function _move_mps_defect(psi::MPSKit.FiniteMPS, from::Int, to::Int)
    from == to && return psi
    tensors = _mps_tensors(psi)
    rng = from < to ? (from:(to - 1)) : ((from - 1):-1:to)
    for a in rng; _swap_mps_sites!(tensors, a); end
    return MPSKit.FiniteMPS(tensors; normalize = false)
end

# ---- inserting / removing a (1-D) defect site (non-gauge-invariant; the total
#      charge of the MPS changes by `charge`) -----------------------------------
# insert a charge-`charge` defect site at MPSKit index `pos`: append it at the right
# boundary (a 1-D charge-shifted isometry) then bump it left into place.
function _insert_mps_defect(psi::MPSKit.FiniteMPS, pos::Int, charge::Int)
    tensors = _mps_tensors(psi)
    bN = MPSKit.right_virtualspace(psi, length(psi))
    Pd = U1Space(charge => 1)
    push!(tensors, isometry(ComplexF64, bN ⊗ Pd, fuse(bN ⊗ Pd)))
    return _move_mps_defect(MPSKit.FiniteMPS(tensors; normalize = false), length(tensors), pos)
end

# remove the (1-D) defect site at MPSKit index `pos`
function _remove_mps_defect(psi::MPSKit.FiniteMPS, pos::Int)
    tensors = _mps_tensors(_move_mps_defect(psi, pos, length(psi)))   # move to the far right
    pop!(tensors)                                                     # drop the boundary defect
    return MPSKit.FiniteMPS(tensors; normalize = false)
end

# ---- absorbing a defect: fuse its (inert, 1-D) site into the matter site to its right,
#      giving a pure-matter chain with one enlarged site. Lets the bond dimension grow
#      freely during evolution/optimisation; split back out afterwards. ---------------

# fuse one MPO channel-matrix entry `(cL ⊗ P) ← (P ⊗ cR)` with the identity on the 1-D
# defect space (defect on the left of the fused leg); scalars/`missing` pass through.
function _fuse_mpo_entry(t, S)
    (t isa Number || ismissing(t)) && return t
    @tensor tf[cL, pf; pfp, cR] := conj(S[pd, p; pf]) * t[cL, p; pp, cR] * S[pd, pp; pfp]
    return tf
end

# the defect-absorbed finite LEMPO: each defect's charge is fused into the matter site at
# its `link` (its right neighbour), so the flux it adds starts on the correct link.
function _fused_defect_lempo(lattice::Lattice, defects, universe::Int)
    isinf(lattice.N) && throw(ArgumentError("fused defect LEMPO requires a finite lattice"))
    A, gauge = _matter_mpo_matrices(lattice; universe = universe)
    link_fcts = collect(Any, gauge.lempo.link_fcts)
    spaces = collect(get_mpskit_spaces(lattice))
    for d in sort(collect(defects); by = x -> x.link)
        (2 ≤ d.link ≤ Int(lattice.N)) ||
            throw(ArgumentError("defect link $(d.link) must lie in 2:$(Int(lattice.N))"))
        idx = (d.link - 1) * lattice.F + 1          # matter site to the right of the defect
        Pd  = U1Space(d.charge => 1); Pm = spaces[idx]
        S   = isometry(Pd ⊗ Pm, fuse(Pd ⊗ Pm))
        A[idx] = map(t -> _fuse_mpo_entry(t, S), A[idx])
        spaces[idx] = fuse(Pd ⊗ Pm)                 # in case another defect lands here
    end
    _trim_finite_mpo!(A, lattice.F)
    return FiniteLEMPOHamiltonian(FiniteMPOHamiltonian(A), link_fcts)
end

# fuse each defect's 1-D site into its right matter neighbour within a FiniteMPS
# (mirrors `_fused_defect_lempo`). Exact (the defect site is entanglement-free).
function _fuse_defect_mps(psi::MPSKit.FiniteMPS, lattice::Lattice, defects)
    isempty(defects) && return psi
    tensors = _mps_tensors(psi)
    for d in sort(collect(defects); by = x -> _defect_position(lattice, defects, x), rev = true)
        pos = _defect_position(lattice, defects, d)        # defect site; matter neighbour at pos+1
        Pd = TensorKit.space(tensors[pos], 2); Pm = TensorKit.space(tensors[pos + 1], 2)
        S  = isometry(Pd ⊗ Pm, fuse(Pd ⊗ Pm))
        @tensor T[vl, pd, pm; vr] := tensors[pos][vl, pd; m] * tensors[pos + 1][m, pm; vr]
        @tensor Fs[vl, pf; vr] := conj(S[pd, pm; pf]) * T[vl, pd, pm; vr]
        tensors = vcat(tensors[1:pos - 1], [Fs], tensors[pos + 2:end])
    end
    return MPSKit.FiniteMPS(tensors; normalize = false)
end

# split each fused site back into (defect, matter) — inverse of `_fuse_defect_mps`.
function _split_defect_mps(psi::MPSKit.FiniteMPS, lattice::Lattice, defects)
    isempty(defects) && return psi
    tensors = _mps_tensors(psi)
    for (i, d) in enumerate(sort(collect(defects); by = x -> x.link))
        midx = (d.link - 1) * lattice.F + 1                # matter index of the fused site
        pos  = midx + (i - 1)                              # shifted by earlier re-insertions
        Pd = U1Space(d.charge => 1); Pm = get_mpskit_spaces(lattice)[midx]
        S  = isometry(Pd ⊗ Pm, fuse(Pd ⊗ Pm))
        @tensor T[vl, pd, pm; vr] := S[pd, pm; pf] * tensors[pos][vl, pf; vr]
        Aq, Cc = TensorKit.left_orth(TensorKit.permute(T, ((1, 2), (3, 4))))   # (vl,pd)←s
        tensors = vcat(tensors[1:pos - 1], [Aq, TensorKit.permute(Cc, ((1, 2), (3,)))], tensors[pos + 1:end])
    end
    return MPSKit.FiniteMPS(tensors; normalize = false)
end
