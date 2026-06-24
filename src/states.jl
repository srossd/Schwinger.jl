"""
`SchwingerState`

Abstract type for Schwinger model states.
"""
abstract type SchwingerState end

"""
    validate_state_operator_compatibility(op, state)

Validate that an operator and state are compatible for expectation value or action.
Throws ArgumentError if L_max or universe don't match.
"""
function validate_state_operator_compatibility(op::SchwingerOperator, state::SchwingerState)
    if op.lattice != lattice(state)
        throw(ArgumentError("Operator and state are defined on different lattices"))
    end
    if !isa(op, MPSKitOperator) && op.L_max != state.hamiltonian.L_max
        throw(ArgumentError("Operator L_max $(op.L_max) does not match state L_max $(state.hamiltonian.L_max)"))
    end
    if op.universe != state.hamiltonian.universe
        throw(ArgumentError("Operator universe $(op.universe) does not match state universe $(state.hamiltonian.universe)"))
    end
    if isa(op, EDOperator) && isa(state, EDState) && op.in_charge != state.net_charge
        throw(ArgumentError("Operator in_charge $(op.in_charge) does not match state net_charge $(state.net_charge)"))
    end
end

"""
`BasisState(occupations, L0)`

A Schwinger model basis state.
"""
struct BasisState <: SchwingerState
    lattice::Lattice
    occupations::BitMatrix
    L₀::Int
    q::Int

    function BasisState(lattice::Lattice, occupations::BitMatrix, L₀::Int, q::Int)
        if size(occupations) != (lattice.N, lattice.F)
            throw(ArgumentError("occupations must be an NxF BitMatrix"))
        end
        if !(0 ≤ sum(occupations) ≤ lattice.N * lattice.F)
            throw(ArgumentError("number of occupied sites must lie in 0:N*F"))
        end
        if q < 1
            throw(ArgumentError("q must be ≥ 1"))
        end
        new(lattice, occupations, L₀, q)
    end
end

function nstates(N::Int, F::Int; L_max::Int = 3, fill::Int = N*F÷2)
    return (2*L_max+1)*binomial(N*F, fill)
end

function lattice(state::BasisState)
    return state.lattice
end

"""
`schwingerbasis(lattice; L_max, q, universe)`

Returns a basis of Schwinger model states.

# Arguments
- `lattice::Lattice`
- `L_max::Int`
- `q::Int`
- `universe::Int`
"""
@memoize function schwingerbasis(lattice::Lattice; L_max::Union{Nothing,Int} = nothing, universe::Int = 0,
                                 net_charge::Int = 0)
    N, F = Int(lattice.N), lattice.F
    if !isnothing(L_max) && L_max < 0
        throw(ArgumentError("L_max must be non-negative"))
    end
    if !isnothing(L_max) && L_max > 0 && !lattice.periodic
        throw(ArgumentError("L_max must be 0 for open boundary conditions"))
    end

    L_max = isnothing(L_max) ? (lattice.periodic ? 3 : 0) : L_max

    # matter charge c relates to occupation by c = q*(Σocc − N*F/2), so a sector of total
    # charge `net_charge` is the set of configurations with this many occupied sites:
    net_charge % lattice.q == 0 ||
        throw(ArgumentError("net_charge $net_charge must be a multiple of q = $(lattice.q)"))
    fill = N*F÷2 + net_charge ÷ lattice.q
    (0 ≤ fill ≤ N*F) ||
        throw(ArgumentError("net_charge $net_charge is out of range for an $N×$F lattice"))

    states = Vector{BasisState}(undef, nstates(Int(N), F; L_max = L_max, fill = fill))
    stateidx = 1
    for L₀ in ((-L_max:L_max) .* lattice.q) .+ universe
        for comb in combinations(1:N*F, fill)
            occupations = BitMatrix(undef, N, F)
            for idx in comb
                occupations[idx] = true
            end
            states[stateidx] = BasisState(lattice, occupations, L₀, lattice.q)
            stateidx += 1
        end
    end
    return states
end

@memoize function positionindex(lattice::Lattice; L_max::Union{Nothing,Int} = nothing, universe::Int = 0,
                                net_charge::Int = 0)
    states = schwingerbasis(lattice; L_max = L_max, universe = universe, net_charge = net_charge)
    return Dict{Tuple{BitMatrix,Int},Int}((states[i].occupations, states[i].L₀) => i for i in eachindex(states))
end

"""
`EDState(hamiltonian, coeffs)`

A Schwinger model state represented as a linear combination of basis states.
"""
struct EDState <: SchwingerState
    hamiltonian::EDOperator
    coeffs::Vector{ComplexF64}
    defects::Vector{DefectCharge}
    net_charge::Int

    function EDState(hamiltonian::EDOperator, coeffs::Vector{ComplexF64},
                     defects::Vector{DefectCharge} = hamiltonian.defects,
                     net_charge::Int = hamiltonian.in_charge)
        new(hamiltonian, coeffs, defects, net_charge)
    end
end

function lattice(state::EDState)
    return state.hamiltonian.lattice
end

# the basis an EDState's coefficients are expressed in (its matter charge sector)
function _edbasis(state::EDState)
    return schwingerbasis(lattice(state); L_max = state.hamiltonian.L_max,
                          universe = state.hamiltonian.universe, net_charge = state.net_charge)
end

function Base.:*(scalar::Number, state::EDState)
    return EDState(state.hamiltonian, scalar * state.coeffs, state.defects, state.net_charge)
end

function Base.:*(state::EDState, scalar::Number)
    return scalar * state
end

"""
`ITensorState(hamiltonian, psi)`

A Schwinger model MPS.
"""
struct ITensorState <: SchwingerState
    hamiltonian::ITensorOperator
    psi::ITensorMPS.MPS
    defects::Vector{DefectCharge}

    function ITensorState(hamiltonian::ITensorOperator, psi::ITensorMPS.MPS,
                          defects::Vector{DefectCharge} = hamiltonian.defects)
        new(hamiltonian, psi, defects)
    end
end

function lattice(state::ITensorState)
    return state.hamiltonian.lattice
end

function Base.:*(scalar::Number, state::ITensorState)
    return ITensorState(state.hamiltonian, scalar * state.psi, state.defects)
end

function Base.:*(state::ITensorState, scalar::Number)
    return scalar * state
end

# =============================================================================
# MPSKit Backend State
# =============================================================================

"""
`MPSKitState(hamiltonian, psi)`

A Schwinger model MPS using MPSKit.jl.
"""
struct MPSKitState <: SchwingerState
    hamiltonian::MPSKitOperator
    psi::MPSKit.AbstractMPS
    defects::Vector{DefectCharge}

    function MPSKitState(hamiltonian::MPSKitOperator, psi::MPSKit.AbstractMPS,
                         defects::Vector{DefectCharge} = hamiltonian.defects)
        if psi isa InfiniteMPS
            isinf(hamiltonian.lattice) || throw(ArgumentError("Hamiltonian and state must both be finite or infinite"))
        elseif psi isa WindowMPS
            isinf(hamiltonian.lattice) || throw(ArgumentError("WindowMPS requires an infinite-lattice Hamiltonian"))
        else  # FiniteMPS: a finite lattice, or a finite window on an infinite background
            isinf(hamiltonian.lattice) ||
                length(psi) == Int(hamiltonian.lattice.N) * hamiltonian.lattice.F + length(defects) ||
                throw(ArgumentError("State length does not match Hamiltonian lattice size"))
        end
        new(hamiltonian, psi, defects)
    end
end

function lattice(state::MPSKitState)
    return state.hamiltonian.lattice
end

# A finite window on an infinite background — either a `WindowMPS`, or a bare `FiniteMPS`
# paired with an infinite-lattice Hamiltonian (the two `wavepacket` representations). Such
# states are indexed over their actual length rather than the (infinite) unit cell.
_isfinitewindow(state::MPSKitState) = isinf(lattice(state).N) && !(state.psi isa MPSKit.InfiniteMPS)

"""
`defects(x)`

Return the list of static `DefectCharge`s carried by an operator or state. The
operator/state keeps the *original* (unshifted) lattice and this defect list,
uniformly across backends — ED/ITensors encode the defects as a θ2π step in the
stored matrix/MPO, MPSKit as an extra (hidden) lattice site.
"""
defects(op::EDOperator)      = op.defects
defects(op::ITensorOperator) = op.defects
defects(op::MPSKitOperator)  = op.defects
defects(state::Union{EDState,ITensorState,MPSKitState}) = state.defects
defects(::SchwingerState) = DefectCharge[]

"""
    translate_defect(state, defect::DefectCharge, new_link::Int)

Move the static `defect` (which must be in `defects(state)`) to link `new_link`,
returning a new state with the static charge relocated. The matter wavefunction is
unchanged — only the charge's position moves — so `occupations` are preserved while
`electricfields` shift. For ED/ITensors this merely updates the defect list (the
coefficients/MPS are untouched); for MPSKit the extra defect site is physically moved
within the MPS by an exact fuse/split, so the result is again a valid MPS.

The returned state keeps the original `hamiltonian`; build
`Hamiltonian(lattice(state); defects = defects(state))` if you need the matching
Hamiltonian for energy/evolution at the new charge location.
"""
function translate_defect(state::EDState, defect::DefectCharge, new_link::Int)
    return EDState(state.hamiltonian, state.coeffs, _edit_defect(state.defects, defect, DefectCharge(new_link, defect.charge)), state.net_charge)
end
function translate_defect(state::ITensorState, defect::DefectCharge, new_link::Int)
    return ITensorState(state.hamiltonian, state.psi, _edit_defect(state.defects, defect, DefectCharge(new_link, defect.charge)))
end
function translate_defect(state::MPSKitState, defect::DefectCharge, new_link::Int)
    newdefs = _edit_defect(state.defects, defect, DefectCharge(new_link, defect.charge))
    lat  = state.hamiltonian.lattice
    from = _defect_position(lat, state.defects, defect)
    to   = _defect_position(lat, newdefs, DefectCharge(new_link, defect.charge))
    return MPSKitState(state.hamiltonian, _move_mps_defect(state.psi, from, to), newdefs)
end

"""
    insert_defect(state, link, charge)

Insert a static `DefectCharge(link, charge)` into `state`, returning a new state.
This is the way to begin/end a Wilson line and is **non-gauge-invariant**: the total
charge changes by `charge`. ED/ITensors merely add the defect to the list (the
coefficients/MPS are untouched); MPSKit inserts the extra 1-D defect site at `link`
and bumps the virtual spaces of all tensors to its right.
"""
insert_defect(state::EDState, link::Int, charge::Int) =
    EDState(state.hamiltonian, state.coeffs, vcat(state.defects, DefectCharge(link, charge)), state.net_charge)
insert_defect(state::ITensorState, link::Int, charge::Int) =
    ITensorState(state.hamiltonian, state.psi, vcat(state.defects, DefectCharge(link, charge)))
function insert_defect(state::MPSKitState, link::Int, charge::Int)
    newdefs = vcat(state.defects, DefectCharge(link, charge))
    pos = _defect_position(state.hamiltonian.lattice, newdefs, DefectCharge(link, charge))
    return MPSKitState(state.hamiltonian, _insert_mps_defect(state.psi, pos, charge), newdefs)
end

"""
    remove_defect(state, defect::DefectCharge)

Remove `defect` (which must be present) from `state`; inverse of `insert_defect`.
"""
remove_defect(state::EDState, defect::DefectCharge) =
    EDState(state.hamiltonian, state.coeffs, _edit_defect(state.defects, defect, nothing), state.net_charge)
remove_defect(state::ITensorState, defect::DefectCharge) =
    ITensorState(state.hamiltonian, state.psi, _edit_defect(state.defects, defect, nothing))
function remove_defect(state::MPSKitState, defect::DefectCharge)
    pos = _defect_position(state.hamiltonian.lattice, state.defects, defect)
    return MPSKitState(state.hamiltonian, _remove_mps_defect(state.psi, pos), _edit_defect(state.defects, defect, nothing))
end

function Base.:*(scalar::Number, state::MPSKitState)
    return MPSKitState(state.hamiltonian, scalar * state.psi, state.defects)
end

function Base.:*(state::MPSKitState, scalar::Number)
    return scalar * state
end

struct MPSKitQPState <: SchwingerState
    hamiltonian::MPSKitOperator
    psi::MPSKit.QP

    function MPSKitQPState(hamiltonian::MPSKitOperator, psi::MPSKit.QP)
        isinf(hamiltonian.lattice) || throw(ArgumentError("Quasiparticle ansatz only supported for infinite lattices"))
        new(hamiltonian, psi)
    end
end

function lattice(state::MPSKitQPState)
    return state.hamiltonian.lattice
end

function Base.:*(scalar::Number, state::MPSKitQPState)
    return MPSKitQPState(state.hamiltonian, scalar * state.psi)
end

function Base.:*(state::MPSKitQPState, scalar::Number)
    return scalar * state
end

# =============================================================================
# Ground State Finding
# =============================================================================

"""
`loweststates(hamiltonian, nstates)`

Returns the lowest few eigenstates of the Schwinger model Hamiltonian.

# Arguments
- `hamiltonian::MPOOperator`: Schwinger model Hamiltonian.
- `nstates::Int`:: number of states to determine.
- `momentum::Union{Real,Nothing} = nothing`: for an infinite lattice (MPSKit backend), the
  physical momentum (in units of the coupling `g`, like the rest of the code) at which to
  build the excited (quasiparticle) states; defaults to `0`. Internally converted to MPSKit's
  dimensionless lattice momentum `a·p`. Only supported on infinite lattices — passing a
  non-`nothing` value for a finite lattice throws an `ArgumentError`.
"""
function loweststates(hamiltonian::ITensorOperator, nstates::Int;
    maxiters::Int = 500, initiallinkdim::Int = 4, maxlinkdim::Int = 600, energy_tol::Real = 1E-6, weight::Real = 100, outputlevel::Int = 0, minsweeps::Int = 5,
    momentum::Union{Real,Nothing} = nothing)

    isnothing(momentum) || throw(ArgumentError("Momentum-resolved excitations not supported for finite lattices"))
    H = hamiltonian.mpo
    N, F = hamiltonian.lattice.N, hamiltonian.lattice.F

    states = Vector{ITensorState}(undef, nstates)
    for idx in 1:nstates
        state = [n == N * F + 1 ? hamiltonian.L_max + 1 : isodd(floor((n-1)/F)) ? "Up" : "Dn" for n=1:(N * F + (hamiltonian.lattice.periodic ? 1 : 0))]
        psi = random_mps(get_sites(hamiltonian.lattice; L_max = hamiltonian.L_max), state; linkdims = initiallinkdim)

        sweeps = Sweeps(maxiters)
        maxdim!(sweeps, maxlinkdim)

        if outputlevel >= 1
            println("Running DMRG for state $idx")
        end
        _, psi = ITensorMPS.dmrg(H, [state.psi for state in states[1:idx-1]], psi, sweeps; weight = weight, observer = DMRGObserver(;energy_tol = energy_tol,minsweeps=minsweeps), outputlevel = outputlevel)
        if outputlevel >= 1
            println()
        end

        states[idx] = ITensorState(hamiltonian, psi)
    end

    return states
end

"""
`lowest_states(hamiltonian, nstates)`

Returns the lowest few eigenstates of the Schwinger model Hamiltonian.

# Arguments
- `hamiltonian::EDOperator`: Schwinger model Hamiltonian.
- `nstates::Int`:: number of states to determine.
"""
function loweststates(hamiltonian::EDOperator, nstates::Int; ncv::Int = max(20, 2*nstates + 1),
    momentum::Union{Real,Nothing} = nothing)
    isnothing(momentum) || throw(ArgumentError("Momentum-resolved excitations not supported for finite lattices"))
    _, vecs = if nstates + 2 < size(hamiltonian.matrix)[1]
        try
            Arpack.eigs(hamiltonian.matrix; nev=nstates, ncv = ncv, which=:SR)
        catch e
            if e isa Arpack.XYAUPD_Exception
                @info "Arpack.eigs failed to converge, increasing ncv to $(2*ncv)"
                return loweststates(hamiltonian, nstates; ncv = 2*ncv)
            else
                rethrow(e)
            end
        end
    else
        eigs = eigen(Matrix(hamiltonian.matrix))
        (eigs.values[1:nstates], eigs.vectors[:,1:nstates])
    end
    return [EDState(hamiltonian, vecs[:,n]) for n=1:nstates]
end


"""
`_reflect_vacuum(psi, n)`

Spatial-parity reflection of a θ=π vacuum `InfiniteMPS`, returning the other degenerate vacuum.

At θ = π the model has two degenerate vacua with background bond charge `n` and `n+1`, related by
parity. On the lattice the reflection is a one-site translation (sublattice swap) together with
the charge map (virtual `q → -q + (2n+1)`, physical `q → -q`). Crucially it is a **unitary**
operation — the tensor data is copied, NOT complex-conjugated. (Charge conjugation `C` is the same
charge map but antiunitary, i.e. with complex conjugation; that breaks the k → -k symmetry of the
soliton, giving a lopsided dispersion, so we use the unitary reflection instead.)

Implemented by remapping fusion channels: for each native fusion tree of the target tensor, the
data block of the matching source channel — identified by its coupled (flux) charge, which is
independent of the dual-leg convention — is copied in.

Building `v2 = P(v1)` instead of solving the second vacuum independently fixes the two vacua's
relative phase (removing the soliton's phase ambiguity), gives a symmetric soliton dispersion with
its minimum at k=0, and avoids a second VUMPS solve.
"""
function _reflect_vacuum(psi::MPSKit.InfiniteMPS, n::Int)
    Uc    = length(psi)
    shift = 2n + 1                        # virtual: q → -q + shift ;  physical: q → -q
    function reflect(A)
        vLs, phs, vRs = TensorKit.space(A, 1), TensorKit.space(A, 2), TensorKit.space(A, 3)
        newvL = U1Space((U1Irrep(-c.charge + shift) => TensorKit.dim(vLs, c) for c in TensorKit.sectors(vLs))...)
        newph = U1Space((U1Irrep(-c.charge)         => TensorKit.dim(phs, c) for c in TensorKit.sectors(phs))...)
        # vR is the (dual) right-virtual leg, so its block dim lives at the conjugate charge —
        # map c → c + shift here (not -c + shift as for the nondual codomain legs).
        newvR = U1Space((U1Irrep(c.charge + shift) => TensorKit.dim(vRs, c) for c in TensorKit.sectors(vRs))...)
        B = zeros(ComplexF64, (newvL ⊗ newph) ← newvR)
        targets = Dict{NTuple{3,Int},Any}()
        for (g1, g2) in TensorKit.fusiontrees(B)
            targets[(g1.coupled.charge, g1.uncoupled[1].charge, g1.uncoupled[2].charge)] = (g1, g2)
        end
        for (f1, f2) in TensorKit.fusiontrees(A)
            key = (-f1.coupled.charge + shift, -f1.uncoupled[1].charge + shift, -f1.uncoupled[2].charge)
            haskey(targets, key) || error("_reflect_vacuum: unmatched fusion channel $key")
            g1, g2 = targets[key]
            B[g1, g2] = copy(A[f1, f2])   # unitary reflection: copy, do NOT conjugate
        end
        return B
    end
    # one-site translation (sublattice swap) composed with the per-site reflection
    return MPSKit.InfiniteMPS([reflect(psi.AL[mod1(i + 1, Uc)]) for i in 1:Uc])
end


"""
`lowest_states(hamiltonian, nstates)`

Returns the lowest few eigenstates of the Schwinger model Hamiltonian using MPSKit.

# Arguments
- `hamiltonian::MPSKitOperator`: Schwinger model Hamiltonian.
- `nstates::Int`:: number of states to determine.

# Keywords
- `momentum`: physical momentum (units of the coupling g) for the excitations on an infinite
  lattice. May be a single value or a **list** of momenta. With a single value, each excited
  entry of the result is one quasiparticle; with a list, each excited entry is a list of
  quasiparticles at those momenta (so they can be combined into multi-particle states).
- `svdcut::Bool` (default `true`): on an infinite lattice, adapt the VUMPS bond dimension each
  iteration via a two-site SVD update (`VUMPSSvdCut`) truncated at relative tolerance `cutoff`,
  growing from `initiallinkdim`. Set `false` to run VUMPS at the fixed `initiallinkdim` instead.
- `solitons::Bool`: on an infinite θ=π lattice, return θ=π domain-wall (soliton) states.
  `result[1]` is the vacuum pair `(v1, v2)`; for `k ≥ 2`, `result[k]` is the `(soliton,
  antisoliton)` pair of the (k-1)-th band — or, for a momentum list, a list of such pairs (one
  per momentum), all built on the same vacua so any pair shares backgrounds.
"""
function loweststates(hamiltonian::MPSKitOperator, nstates::Int;
    maxiters::Int = 500, initiallinkdim::Int = 10, bonddim::Union{Nothing,Int} = nothing, initial_Lmax::Int = 3, energy_tol::Real = 1E-6, cutoff::Real = 1E-10, weight::Real = 100., verbose::Bool = false, momentum::Union{Real, Nothing, AbstractVector} = nothing,
    solitons::Bool = false, attenuation::Real = 1e-3, svdcut::Bool = true)

    initiallinkdim = something(bonddim, initiallinkdim)   # `bonddim` is an alias

    H = hamiltonian.lempo
    spaces = get_mpskit_spaces(hamiltonian.lattice)   # defect sites are added in the finite branch
    total_defect = sum(d.charge for d in hamiltonian.defects; init = 0)

    # VUMPS finalizer: with `svdcut`, after each iteration adapt the virtual spaces with a two-site
    # SVD update (`VUMPSSvdCut`), truncating singular values at relative tolerance `cutoff` — so the
    # bond dimension grows from the small `initiallinkdim` to whatever the state needs. With
    # `svdcut = false`, the finalizer is a no-op and VUMPS runs at the fixed `initiallinkdim`.
    vumps_finalize = svdcut ?
        ((it, st, op, ev) -> MPSKit.changebonds(st, op, MPSKit.VUMPSSvdCut(; trscheme = trunctol(; rtol = cutoff)), ev)) :
        ((it, st, op, ev) -> (st, ev))

    # Soliton sector at θ = π (mod 2π): the model has two degenerate vacua v1, v2 whose
    # background electric field sits near U1Irrep(n) and U1Irrep(n+1), with n = -1/2 - θ/2π.
    # A soliton is a domain wall from v1 to v2 (a quasiparticle built on the two vacua as
    # left/right backgrounds); the antisoliton is the reverse wall v2 → v1, built on the
    # SAME vacuum objects so the two share backgrounds (e.g. for a soliton/antisoliton pair).
    #
    # `result[1]` is the vacuum pair `(v1, v2)`; for k ≥ 2, `result[k]` is the pair
    # `(soliton, antisoliton)` of the (k-1)-th excitation.
    if solitons
        isinf(hamiltonian.lattice.N) ||
            throw(ArgumentError("solitons are only available on an infinite lattice"))
        θ2π = hamiltonian.lattice.θ2π[1]
        isapprox(mod(θ2π, 1.0), 0.5; atol = 1e-8) ||
            throw(ArgumentError("solitons require θ = π (mod 2π), i.e. θ2π ≡ 1/2 (mod 1); got θ2π = $θ2π"))

        n  = round(Int, -0.5 - θ2π)
        Uc = length(spaces)
        ψ₀ = MPSKit.InfiniteMPS(spaces, [U1Space([q => initiallinkdim for q in hamiltonian.lattice.q*(-initial_Lmax:initial_Lmax)]) for _ in 1:Uc])
        alg = MPSKit.VUMPS(; maxiter = maxiters, tol = energy_tol, verbosity = verbose ? 1 : 0, finalize = vumps_finalize)
        # Solve the first vacuum (biased to background charge n by `attenuateLinks`), then obtain
        # the second as its parity reflection. Using P rather than a second independent VUMPS solve
        # fixes the two vacua's relative phase (removing the soliton's phase ambiguity), gives a
        # symmetric soliton dispersion (min at k=0), and is cheaper. See `_reflect_vacuum`.
        ψ1, envs1, _ = MPSKit.find_groundstate(attenuateLinks(ψ₀, fill(U1Irrep(n), Uc), attenuation), H, alg)
        ψ2 = _reflect_vacuum(ψ1, n)
        envs2 = MPSKit.environments(ψ2, H)

        v1, v2 = MPSKitState(hamiltonian, ψ1), MPSKitState(hamiltonian, ψ2)
        # The two θ=π vacua must differ (opposite background fields, ±1/2). If `attenuateLinks`
        # failed to separate them, VUMPS lands on the same vacuum and the "soliton" is spurious.
        isapprox(electricfields(v1), electricfields(v2); atol = 1e-3) &&
            @warn "soliton vacua have nearly identical electric fields — attenuation did not select \
                   the two distinct θ=π vacua, so the soliton is likely spurious. Try a smaller \
                   `attenuation` (smaller = stronger)."

        result = Vector{Any}(undef, nstates)
        result[1] = (v1, v2)
        if nstates > 1
            isnothing(momentum) && (momentum = 0.0)
            islist = momentum isa AbstractVector
            moms   = islist ? collect(momentum) : [momentum]
            alg_qp = MPSKit.QuasiparticleAnsatz()
            # For each requested momentum, solve the soliton band (domain wall v1 → v2) and the
            # antisoliton band (the reverse v2 → v1). All are built on the SAME vacua, so any
            # (soliton, antisoliton) pair — even at different momenta — shares backgrounds.
            bands = map(moms) do mom
                k_lat = hamiltonian.lattice.a * mom
                _, sols  = MPSKit.excitations(H, alg_qp, k_lat, ψ1, envs1, ψ2, envs2; num = nstates - 1)
                _, antis = MPSKit.excitations(H, alg_qp, k_lat, ψ2, envs2, ψ1, envs1; num = nstates - 1)
                (sols, antis)
            end
            for k in 2:nstates
                # one (soliton, antisoliton) pair per requested momentum
                pairs = [(MPSKitQPState(hamiltonian, bands[mi][1][k-1]),
                          MPSKitQPState(hamiltonian, bands[mi][2][k-1])) for mi in eachindex(moms)]
                result[k] = islist ? pairs : pairs[1]
            end
        end
        return result
    end

    states = Vector{Any}(undef, nstates)   # entries may be a QP, or (for a momentum list) a list of QPs
    if isinf(hamiltonian.lattice.N)
        ψ₀ = MPSKit.InfiniteMPS(spaces, [U1Space([q => initiallinkdim for q in hamiltonian.lattice.q*(-initial_Lmax:initial_Lmax)]) for _ in 1:length(spaces)]) #TODO: add attenuation near theta = pi
        if abs(hamiltonian.lattice.θ2π[1] - 0.5) < 0.1
            ψ₀ = attenuateLinks(ψ₀, hamiltonian.lattice.θ2π[1] < 0.5 ? [U1Irrep(0), U1Irrep(0)] : [U1Irrep(-1), U1Irrep(-1)], 0.01)
        end
        alg = MPSKit.VUMPS(; maxiter=maxiters, tol=energy_tol, verbosity=verbose ? 1 : 0, finalize = vumps_finalize)
        ψ, envs, _ = MPSKit.find_groundstate(ψ₀, H, alg)
        states[1] = MPSKitState(hamiltonian, ψ)
        
        if nstates > 1
            isnothing(momentum) && (momentum = 0.0)
            islist = momentum isa AbstractVector
            moms   = islist ? collect(momentum) : [momentum]
            algqp  = MPSKit.QuasiparticleAnsatz()
            # `momentum` is the physical momentum (units of the coupling g); MPSKit's excitation
            # momentum is the dimensionless lattice momentum a·p (phase per site). A list of
            # momenta returns, per excited band, a list of quasiparticles at those momenta.
            bands = map(mom -> MPSKit.excitations(H, algqp, hamiltonian.lattice.a * mom, ψ, envs; num = nstates - 1)[2], moms)
            for n in 2:nstates
                qps = [MPSKitQPState(hamiltonian, bands[mi][n-1]) for mi in eachindex(moms)]
                states[n] = islist ? qps : qps[1]
            end
        end
    else
        # With static defects, optimise in the absorbed representation (defect sites fused
        # into their matter neighbours) so the bond dimension can grow, then split back out.
        if isempty(hamiltonian.defects)
            Huse, spc, split = H, collect(spaces), identity
        else
            lat = hamiltonian.lattice
            Huse = _fused_defect_lempo(lat, hamiltonian.defects, hamiltonian.universe)
            spc = collect(get_mpskit_spaces(lat))
            for d in hamiltonian.defects
                idx = (d.link - 1) * lat.F + 1
                spc[idx] = fuse(U1Space(d.charge => 1) ⊗ spc[idx])
            end
            split = ψ -> _split_defect_mps(ψ, lat, hamiltonian.defects)
        end
        ψ₀ = MPSKit.FiniteMPS(rand, ComplexF64, spc, U1Space([q => initiallinkdim for q in hamiltonian.lattice.q*(-initial_Lmax:initial_Lmax)]); right = U1Space(total_defect => 1))
        alg = MPSKit.DMRG2(; maxiter=maxiters, tol=energy_tol, trscheme = trunctol(; rtol = cutoff), verbosity=verbose ? 1 : 0)
        ψ, = MPSKit.find_groundstate(ψ₀, Huse, alg)
        states[1] = MPSKitState(hamiltonian, split(ψ), hamiltonian.defects)

        if nstates > 1
            isnothing(momentum) || throw(ArgumentError("Momentum-resolved excitations not supported for finite lattices"))
            alg2 = MPSKit.FiniteExcited(alg, weight)
            _, psis = MPSKit.excitations(Huse, alg2, (ψ,); num = nstates - 1)
            for n in 2:nstates
                states[n] = MPSKitState(hamiltonian, split(psis[n-1]), hamiltonian.defects)
            end
        end
    end

    return states
end

"""
`groundstate(hamiltonian)`
Returns the ground state of the Schwinger model Hamiltonian.

# Arguments
- `hamiltonian::SchwingerOperator`: Schwinger model Hamiltonian.
"""
function groundstate(hamiltonian::SchwingerOperator; kwargs...)
    return loweststates(hamiltonian, 1; kwargs...)[1]
end

"""
`energygap(hamiltonian)`
Returns the energy difference between the lowest two states of the Hamiltonian.

# Arguments
- `hamiltonian::SchwingerOperator`: Schwinger model Hamiltonian.
"""
function energygap(hamiltonian::SchwingerOperator; kwargs...)
    if isinf(hamiltonian.lattice.N)
        return energy(loweststates(hamiltonian, 2; kwargs...)[2]) # energy of first QP is the gap
    end
    return abs(-(map(energy, loweststates(hamiltonian, 2; kwargs...))...))
end

"""
`expectation(op, state)`

Return the expectation value of the operator `op` in `state`.

# Arguments
- `op::EDOperator``: operator.
- `state::EDState`: state.
"""
function expectation(op::EDOperator, state::EDState)
    validate_state_operator_compatibility(op, state)
    return dot(state.coeffs, op.matrix * state.coeffs)/real(dot(state, state))
end

"""
`expectation(op, state)`

Return the expectation value of the operator `op` in `state`.

# Arguments
- `op::ITensorOperator``: operator.
- `state::ITensorState`: state.
"""
function expectation(op::ITensorOperator, state::ITensorState)
    validate_state_operator_compatibility(op, state)
    return ITensors.inner(state.psi', op.mpo, state.psi)/real(dot(state, state))
end

"""
`expectation(op, state)`

Return the expectation value of the operator `op` in `state`.

# Arguments
- `op::MPSKitOperator``: operator.
- `state::MPSKitState`: state.
"""
function expectation(op::MPSKitOperator, state::MPSKitState)
    validate_state_operator_compatibility(op, state)
    return MPSKit.expectation_value(state.psi, op.lempo)
end

"""
`dot(psi1::EDState, psi2::EDState)`

Return the inner product of two Schwinger model states.

# Arguments
- `psi1::EDState`: state.
- `psi2::EDState`: state.
"""
function LinearAlgebra.dot(psi1::EDState, psi2::EDState)
    return dot(psi1.coeffs, psi2.coeffs)
end

"""
`dot(psi1::ITensorState, psi2::ITensorState)`

Return the inner product of two Schwinger model states.

# Arguments
- `psi1::ITensorState`: state.
- `psi2::ITensorState`: state.
"""
function LinearAlgebra.dot(psi1::ITensorState, psi2::ITensorState)
    return ITensors.inner(psi1.psi, psi2.psi)
end

"""
`dot(psi1::MPSKitState, psi2::MPSKitState)`

Return the inner product of two Schwinger model states.

# Arguments
- `psi1::MPSKitState`: state.
- `psi2::MPSKitState`: state.
"""
function LinearAlgebra.dot(psi1::MPSKitState, psi2::MPSKitState)
    return dot(psi1.psi, psi2.psi)
end

"""
`act(op, state)`

Apply the operator `op` to the state `state`.

# Arguments
- `op::ITensorOperator`: operator.
- `state::ITensorState`: state.
"""
function act(op::ITensorOperator, state::ITensorState)
    validate_state_operator_compatibility(op, state)
    return ITensorState(state.hamiltonian, apply(op.mpo, state.psi), state.defects)
end

function Base.:*(op::ITensorOperator, state::ITensorState)
    return act(op, state)
end

"""
`act(op, state)`

Apply the operator `op` to the state `state`.

# Arguments
- `op::EDOperator`: operator.
- `state::EDState`: state.
"""
function act(op::EDOperator, state::EDState)
    validate_state_operator_compatibility(op, state)
    return EDState(state.hamiltonian, op.matrix * state.coeffs, state.defects, op.out_charge)
end

function Base.:*(op::EDOperator, state::EDState)
    return act(op, state)
end

"""
`act(op, state)`

Apply the operator `op` to the state `state`.

# Arguments
- `op::MPSKitOperator`: operator.
- `state::MPSKitState`: state.
"""
function act(op::MPSKitOperator, state::MPSKitState)
    validate_state_operator_compatibility(op, state)
    return MPSKitState(state.hamiltonian, op.lempo * state.psi, state.defects)
end

function Base.:*(op::MPSKitOperator, state::MPSKitState)
    return act(op, state)
end

"""
`energy(state)`

Return the expectation value of the Hamiltonian.

# Arguments
- `state::SchwingerState`: Schwinger model state (ED, MPS, or MPSKit).
"""
function energy(state::Union{EDState,ITensorState,MPSKitState}; warn::Bool = true)
    if warn && isinf(lattice(state).N)
        @warn "Calculating energy for infinite lattice. Consider using energy_density instead."
    end
    return real(expectation(state.hamiltonian, state))
end

function energy(state::MPSKitQPState; warn::Bool = true)
    return real(MPSKitLEMPO.expectation_gap(state.psi, state.hamiltonian.lempo))
end

"""
`energy_density(state)`

Return the energy density (energy per unit length) of the state.
# Arguments
- `state::SchwingerState`: Schwinger model state (ED, MPS, or MPSKit).
"""
function energy_density(state::Union{EDState,ITensorState,MPSKitState})
    if isinf(lattice(state).N)
        # For infinite lattices, energy gives energy of two-site unit cell
        return energy(state, warn = false)/(2*lattice(state).a)
    end
    return energy(state)/lattice(state).L
end

"""
`energy_density(state, site)`

Return the energy density (energy per unit length) at `site`: the electric energy on link
`site`, the mass on site `site`, and the hopping (and hopping-mass) on the bond
`(site, site+1)`, all divided by the lattice spacing `a`. Summing over all sites and
multiplying by `a` gives the total energy (see [`energy_densities`](@ref)).

# Arguments
- `state::SchwingerState`: Schwinger model state.
- `site::Int`: site.
"""
# The mass operator on `site` (an on-site term).
_mass_op(state::EDState, site::Int) = EDMass(lattice(state), site; bare = false,
    L_max = state.hamiltonian.L_max, universe = state.hamiltonian.universe, charge = state.net_charge)
_mass_op(state::ITensorState, site::Int) = ITensorMass(lattice(state), site; bare = false,
    L_max = state.hamiltonian.L_max, universe = state.hamiltonian.universe)
_mass_op(state::MPSKitState, site::Int) = MPSKitMass(lattice(state), site; bare = false,
    universe = state.hamiltonian.universe)

# The energy on bond/link ℓ: the electric energy on link ℓ plus the hopping and
# hopping-mass on bond (ℓ, ℓ+1) — everything that lives *between* sites ℓ and ℓ+1.
function _bond_energy(state::EDState, ℓ::Int)
    lat = lattice(state); kw = (; L_max = state.hamiltonian.L_max,
                                universe = state.hamiltonian.universe, charge = state.net_charge)
    return _gauge_link_energy(state, ℓ) +
           real(expectation(EDHopping(lat, ℓ; bare = false, kw...), state)) +
           real(expectation(EDHoppingMass(lat, ℓ; bare = false, kw...), state))
end
function _bond_energy(state::ITensorState, ℓ::Int)
    lat = lattice(state); kw = (; L_max = state.hamiltonian.L_max, universe = state.hamiltonian.universe)
    return _gauge_link_energy(state, ℓ) +
           real(expectation(ITensorHopping(lat, ℓ; bare = false, kw...), state)) +
           real(expectation(ITensorHoppingMass(lat, ℓ; bare = false, kw...), state))
end
function _bond_energy(state::MPSKitState, ℓ::Int)
    lat = lattice(state); u = state.hamiltonian.universe; N = Int(lat.N)
    e = _gauge_link_energy(state, ℓ) + real(expectation(MPSKitHoppingMass(lat, ℓ; bare = false, universe = u), state))
    ℓ < N && (e += sum(real(expectation(op, state)) for op in MPSKitHopping(lat, ℓ; bare = false, universe = u)))
    return e
end

# expectation of the electric energy (a/2)(Lₗᵢₙₖ+θ)² on a single `link`
_gauge_link_energy(state::EDState, link::Int) = real(expectation(
    EDGaugeKinetic(lattice(state), link; bare = false, L_max = state.hamiltonian.L_max,
                   universe = state.hamiltonian.universe, charge = state.net_charge), state))
_gauge_link_energy(state::ITensorState, link::Int) = real(expectation(
    ITensorGaugeKinetic(lattice(state), link; bare = false, L_max = state.hamiltonian.L_max,
                        universe = state.hamiltonian.universe), state))
_gauge_link_energy(state::MPSKitState, link::Int) = real(expectation(
    MPSKitGaugeKinetic(lattice(state), link; universe = state.hamiltonian.universe), state))

# Average a per-bond quantity `b(ℓ)` (electric + hopping + hopping-mass, all of which live
# between sites ℓ and ℓ+1) over the two bonds neighboring `site`: an interior site gets
# ½·(left bond + right bond); a bond with a single neighboring site (a true open boundary,
# `boundary = true` for the rightmost bond of a finite lattice) is assigned in full, so the
# per-site energies still sum to the total. Mass is on-site and is added separately.
function _averaged_bond(b, site::Int, N::Int, boundary::Bool)
    e = site ≥ 2 ? 0.5 * b(site - 1) : 0.0
    e += site < N ? 0.5 * b(site) : (boundary ? b(site) : 0.0)
    return e
end

function energy_density(state::Union{EDState,ITensorState,MPSKitState}, site::Int)
    if state isa MPSKitState && _isfinitewindow(state)
        return _energy_densities_window(state)[site]   # wavepacket on an infinite background
    end
    isinf(lattice(state).N) && throw(ArgumentError("per-site energy_density requires a finite lattice"))
    N = Int(lattice(state).N)
    # energy per unit length: per-site energy / a (so Σ energy_densities · a = total energy)
    return (real(expectation(_mass_op(state, site), state)) +
            _averaged_bond(ℓ -> _bond_energy(state, ℓ), site, N, true)) / lattice(state).a
end

# The left/right "wings" of a wavepacket `WindowMPS` are the infinite vacuum. Expose them as
# ordinary (infinite) `MPSKitState`s so existing observables (e.g. `energy_density`) can be
# evaluated on them — used to supply the boundary bond that the finite window otherwise drops.
function _window_vacua(state::MPSKitState)
    ψ = state.psi
    ψ isa WindowMPS || throw(ArgumentError("_window_vacua requires a WindowMPS-backed state"))
    return MPSKitState(state.hamiltonian, ψ.left_gs), MPSKitState(state.hamiltonian, ψ.right_gs)
end

# Per-site energy of a wavepacket `WindowMPS` (infinite background), for the whole window.
# The gauge-integrated electric energy is non-local, but the per-link `(a/2)⟨(Lₙ+θ)²⟩` comes
# from the bond (flux) distribution via `link_expectation`, and the on-site mass / two-site
# hopping (+ hopping-mass) from local expectation values; the mass is on-site and the bond
# energy is averaged over the two bonds neighboring each site. The wavepacket is unnormalized,
# so the norm `dot(ψ,ψ)` is computed once and divided out (rather than per `expectation_value`
# call, which recomputes it — the dominant cost). Operators are cached by site/bond parity and
# contracted directly via `contract_mpo_expval1/2`. Currently supports F = 1.
# (`energy_density(state, site)` on a window indexes into this.)
function _energy_densities_window(state::MPSKitState)
    lat = lattice(state); ψ = state.psi
    lat.F == 1 || throw(ArgumentError("energy_density on a window currently supports F = 1"))
    W = length(ψ); u = state.hamiltonian.universe
    sp = get_mpskit_spaces(lat); P(n) = sp[mod1(n, length(sp))]
    q = lat.q; a = lat.a
    nrm2 = real(dot(ψ, ψ))

    massop(site) = begin
        mop = zeros(ComplexF64, P(site) ← P(site))
        block(mop, U1Irrep(0)) .= -0.5
        block(mop, U1Irrep(isodd(site) ? -q : q)) .= 0.5
        lat.mlat[mod1(site, length(lat.mlat))][1] * mop
    end
    hopcache = Dict{Bool,Any}()                          # hopping operator depends only on bond parity
    hopop(ℓ) = get!(hopcache, isodd(ℓ)) do
        raw = nothing
        for qs in (q, -q)
            openT  = ones(ComplexF64, U1Space(0 => 1) ⊗ P(ℓ)     ← P(ℓ)     ⊗ U1Space(qs => 1))
            closeT = ones(ComplexF64, U1Space(qs => 1) ⊗ P(ℓ + 1) ← P(ℓ + 1) ⊗ U1Space(0 => 1))
            @tensor t[-1 -2; -3 -4] := openT[1, -1; -3, 2] * closeT[2, -2; -4, 1]
            raw = raw === nothing ? t : raw + t
        end
        raw
    end

    masses = [real(MPSKit.contract_mpo_expval1(ψ.AC[site], massop(site))) / nrm2 for site in 1:W]
    bonds = map(1:W-1) do ℓ                              # each internal bond once
        θℓ = Float64(lat.θ2π[mod1(ℓ, length(lat.θ2π))]) + u
        b = link_expectation(ψ, ℓ, r::U1Irrep -> (a / 2) * (r.charge + θℓ)^2) / nrm2
        coeff = 1/(2a) + (-1)^(ℓ + 1) * lat.mprime[mod1(ℓ, length(lat.mprime))][1]
        hop = real(MPSKit.contract_mpo_expval2(ψ.AC[ℓ], ψ.AR[ℓ + 1], hopop(ℓ))) / nrm2
        b + coeff * hop
    end
    eds = [(masses[site] + _averaged_bond(ℓ -> bonds[ℓ], site, W, false)) / a for site in 1:W]   # per unit length

    # The two outermost sites miss the bond into the wing, which spikes their energy density.
    # For a `WindowMPS` the wings are the explicit infinite vacuum, so replace those sites with
    # the wing's vacuum energy density (≈ what a bulk vacuum site carries). This adds energy
    # living partly in the wings, so Σ(eds)·a no longer equals the window energy exactly — the
    # sum rule is therefore not asserted for window states. A bare `FiniteMPS` has no wing (a
    # genuine open boundary), so it is left as is.
    if ψ isa WindowMPS
        lv, rv = _window_vacua(state)
        eds[1]   = real(energy_density(lv))
        eds[end] = real(energy_density(rv))
    end
    return eds
end

"""
`energy_densities(state)`

Return the list of energy densities (energy per unit length) of `state` on sites 1 through
N; their sum times the lattice spacing `a` is the total energy of the state.

# Arguments
- `state::SchwingerState`: Schwinger model state.
"""
function energy_densities(state::SchwingerState)
    if state isa MPSKitState && _isfinitewindow(state)
        return _energy_densities_window(state)            # batched: norm/bonds computed once
    end
    N = isinf(lattice(state).N) ? 2 : Int(lattice(state).N)
    return energy_density.(Ref(state), 1:N)
end

"""
`energy(state)`

Return the expectation value of the Hamiltonian.

# Arguments
- `state::BasisState`: Schwinger model basis state.
"""
function energy(state::BasisState)
    N, F = state.lattice.N, state.lattice.F
    efs = electricfields(state) # includes θ/2π
    occs = occupations(state)

    electricenergy = (state.lattice.a/2) * sum(efs.^2)
    massenergy = sum((-1)^j * state.lattice.mlat[j][k] * occs[j,k] for j=1:N, k=1:F)

    return electricenergy + massenergy
end

"""
`occupation(state, site)`

Return the expectations of χ†χ operators of each flavor on a given site.

# Arguments
- `state::ITensorState`: Schwinger model state.
- `site::Int`: the lattice site.
"""
function occupation(state::ITensorState, site::Int)
    N, F = state.hamiltonian.lattice.N, state.hamiltonian.lattice.F
    if !(1 ≤ site ≤ N)
        throw(ArgumentError("Site must be between 1 and N"))
    end
    psi = state.psi
    return expect(psi, "Sz", sites=((site-1)*F+1):(site*F)) .+ 1/2
end

"""
`occupations(state)`

Return an NxF matrix of the expectations of χ†χ operators on each site.

# Arguments
- `state::ITensorState`: Schwinger model state.
"""
function occupations(state::ITensorState)
    N, F = Int(state.hamiltonian.lattice.N), state.hamiltonian.lattice.F
    psi = state.psi
    return transpose(reshape(expect(psi, "Sz", sites=1:N*F) .+ 1/2, (F,N)))
end

"""
`occupation(state, site)`

Return the expectations of χ†χ operators of each flavor on a given site.

# Arguments
- `state::MPSKitState`: Schwinger model state.
- `site::Int`: the lattice site.
"""
function occupation(state::MPSKitState, site::Int)
    N, F = state.hamiltonian.lattice.N, state.hamiltonian.lattice.F
    psi = state.psi
    if _isfinitewindow(state)
        N_sites = length(psi) ÷ F
        1 ≤ site ≤ N_sites || throw(ArgumentError("Site must be between 1 and $N_sites"))
    elseif isinf(N) && !(1 ≤ site ≤ 2)
        return occupation(state, site % 2 == 0 ? 2 : 1)
    elseif isfinite(N) && !(1 ≤ site ≤ Int(N))
        throw(ArgumentError("Site must be between 1 and N"))
    end

    lat = state.hamiltonian.lattice
    space = get_mpskit_spaces(lat)[isodd(site) ? 1 : F + 1]
    occop = zeros(space ← space)
    if isodd(site)
        block(occop, U1Irrep(0)) .= 1.0
    else
        block(occop, U1Irrep(lattice(state).q)) .= 1.0
    end

    occs = zeros(F)
    for k in 1:F
        idx = _matter_index(lat, state.defects, site, k)   # skip defect sites
        occs[k] = real(MPSKit.expectation_value(psi, idx => occop))
    end
    return occs
end

"""
`occupations(state)`

Return an NxF matrix of the expectations of χ†χ operators on each site.

# Arguments
- `state::MPSKitState`: Schwinger model state.
"""
function occupations(state::MPSKitState)
    lat = state.hamiltonian.lattice
    defects = state.defects
    N, F = lat.N, lat.F
    psi = state.psi
    N = _isfinitewindow(state) ? (length(psi) - length(defects)) ÷ F : (isinf(N) ? 2 : Int(N))
    occs = zeros(N, F)

    # matter physical spaces (independent of the inserted defect sites)
    matterspaces = get_mpskit_spaces(lat)
    occopodd  = zeros(matterspaces[1]   ← matterspaces[1])
    occopeven = zeros(matterspaces[F+1] ← matterspaces[F+1])
    block(occopodd,  U1Irrep(0)) .= 1.0
    block(occopeven, U1Irrep(lat.q)) .= 1.0

    norm2 = real(MPSKit.dot(psi, psi))

    for j in 1:N
        for k in 1:F
            idx = _matter_index(lat, defects, j, k)   # skip inserted defect sites
            occs[j, k] = real(MPSKit.contract_mpo_expval1(psi.AC[idx], (iseven(j) ? occopeven : occopodd), psi.AC[idx]))/norm2
        end
    end
    return occs
end

"""
`occupation(state, site)`

Return the expectations of χ†χ operators of each flavor on a given site.

# Arguments
- `state::BasisState`: Schwinger model basis state.
- `site::Int`: the lattice site.
"""
function occupation(state::BasisState, site::Int)
    if !(1 ≤ site ≤ state.lattice.N)
        throw(ArgumentError("Site must be between 1 and N"))
    end
    return state.occupations[site]
end

"""
`occupations(state)`

Return an NxF matrix of the expectations of χ†χ operators on each site.

# Arguments
- `state::BasisState`: Schwinger model basis state.
"""
function occupations(state::BasisState) 
    return state.occupations
end

"""
`occupation(state, site)`

Return the expectations of χ†χ operators of each flavor on a given site.

# Arguments
- `state::EDState`: Schwinger model basis state.
- `site::Int`: the lattice site.
"""
function occupation(state::EDState, site::Int)
    if !(1 ≤ site ≤ state.hamiltonian.lattice.N)
        throw(ArgumentError("Site must be between 1 and N"))
    end
    states = _edbasis(state)
    occs = zeros(F)
    for (coeff, state) in zip(state.coeffs, states)
        occs .+= real(abs2(coeff)) .* occupations(state)[site]
    end
    return occs/real(dot(state, state))
end

"""
`occupations(state)`

Return an NxF matrix of the expectations of χ†χ operators on each site.

# Arguments
- `state::EDState`: Schwinger model basis state.
"""
function occupations(state::EDState)
    N, F = Int(state.hamiltonian.lattice.N), state.hamiltonian.lattice.F
    states = _edbasis(state)
    occs = zeros(N,F)
    for (coeff, state) in zip(state.coeffs, states)
        occs .+= real(abs2(coeff)) .* occupations(state)
    end
    return occs/real(dot(state, state))
end

"""
`scalardensity(state, site)`

Return the scalar density at site `site`.

# Arguments
- `state::SchwingerState`: Schwinger model state.
- `site::Int`: site.
"""
function scalardensity(state::SchwingerState, site::Int)
    N = lattice(state).N
    if isfinite(N) && !(1 ≤ site ≤ Int(N))
        throw(ArgumentError("site must be between 1 and N"))
    end
    N = isinf(N) ? 2 : Int(N)
    before = site == 1 ? N : site - 1
    after = site == N ? 1 : site + 1
    occs = sum(occupations(state), dims=2)
    return 1/(lattice(state).a)*((-1)^(before) * occs[before]/4 + (-1)^site * occs[site]/2 + (-1)^(after) * occs[after]/4)
end

"""
`scalardensities(state)`

Return the list of scalar densities of `state` on sites 1 through N.

# Arguments
- `state::SchwingerState`: Schwinger model state.
- `site::Int`: site.
"""
function scalardensities(state::SchwingerState)
    N = isinf(lattice(state).N) ? 2 : Int(lattice(state).N)
    return scalardensity.(Ref(state), 1:N)
end

"""
`scalar(state)`

Return the expectation value of the scalar condensate, ⟨H_mass⟩/L.

# Arguments
- `state::SchwingerState`: Schwinger model state.
"""
function scalar(state::SchwingerState)
    return mean(scalardensities(state))
end

"""
`pseudoscalardensity(state, n)`

Return the pseudoscalar density at site `n`.

# Arguments
- `state::SchwingerState`: Schwinger model state.
- `n::Int`: site.
"""
function pseudoscalardensity(state::ITensorState, site::Int)
    N = Int(lattice(state).N)
    if isfinite(N) && !(1 ≤ site ≤ N)
        throw(ArgumentError("site must be between 1 and N"))
    end
    before = site == 1 ? N : site - 1
    function p(n)
        return real(expectation(ITensorHoppingMass(lattice(state), n; bare = true, L_max = state.hamiltonian.L_max, universe = state.hamiltonian.universe), state))
    end
    return 1/(lattice(state).a)*(p(site)/2 + p(before)/2)
end

function pseudoscalardensity(state::EDState, site::Int)
    N = Int(lattice(state).N)
    if isfinite(N) && !(1 ≤ site ≤ N)
        throw(ArgumentError("site must be between 1 and N"))
    end
    before = site == 1 ? N : site - 1
    function p(n)
        return real(expectation(EDHoppingMass(lattice(state), n; bare = true, L_max = state.hamiltonian.L_max, universe = state.hamiltonian.universe, charge = state.net_charge), state))
    end
    return 1/(lattice(state).a)*(p(site)/2 + p(before)/2)
end

function pseudoscalardensity(::BasisState, ::Int)
    return 0
end

function pseudoscalardensity(state::MPSKitState, site::Int)
    if isfinite(lattice(state).N) && !(1 ≤ site ≤ Int(lattice(state).N))
        throw(ArgumentError("site must be between 1 and N"))
    end
    N = isinf(lattice(state).N) ? 2 : Int(lattice(state).N)
    before = site == 1 ? N : site - 1
    function p(n)
        return real(expectation(MPSKitHoppingMass(lattice(state), n; bare = true, universe = state.hamiltonian.universe), state))
    end
    return 1/(lattice(state).a)*(p(site)/2 + p(before)/2)
end

"""
`pseudoscalardensities(state)`

Return the list of pseudoscalar densities of `state` on sites 1 through N.

# Arguments
- `state::SchwingerState`: Schwinger model state.
- `site::Int`: site.
"""
function pseudoscalardensities(state::SchwingerState)
    N = isinf(lattice(state).N) ? 2 : Int(lattice(state).N)
    return pseudoscalardensity.(Ref(state), 1:N)
end

"""
`pseudoscalar(state)`

Return the expectation value of the pseudoscalar condensate, ⟨H_hoppingmass⟩/L.

# Arguments
- `state::SchwingerState`: Schwinger model state.
"""
function pseudoscalar(state::SchwingerState)
    return mean(pseudoscalardensities(state))
end

"""
`charges(state)`

Return a list of the expectations of Q operators on each site.

# Arguments
- `state::SchwingerState`: Schwinger model state.
"""
function charges(state::SchwingerState)
    F    = lattice(state).F
    occs = occupations(state)
    N    = size(occs, 1)
    return (sum(occs, dims=2) + (F .* [-1/2 + (-1)^(n)/2 for n=1:N])) .* lattice(state).q
end

"""
`charge(state, site)`

Return the expectation value of the charge operator on site `site`.

# Arguments
- `state::SchwingerState`: Schwinger model state.
- `site::Int`: site.
"""
function charge(state::SchwingerState, site::Int)
    return charges(state)[site]
end

"""
`L₀(state)`

Return the expectation value of L₀.

# Arguments
- `state::SchwingerState`: Schwinger model state.
"""
function L₀(state::BasisState)
    return state.L₀
end

function L₀(state::ITensorState)
    N, F = lattice(state).N, lattice(state).F
    return (lattice(state).periodic ? lattice(state).q * expect(state.psi, "L0", sites=N*F + 1) : 0) + state.hamiltonian.universe
end

function L₀(state::EDState)
    states = _edbasis(state)
    return sum(abs2(coeff) * L₀(state) for (coeff, state) in zip(state.coeffs, states))
end

function L₀(state::MPSKitState)
    # MPSKitState cannot be periodic
    return state.hamiltonian.universe
end

"""
`electricfield(state, link)`

Return the expectation of (L + θ/2π) on the link `link`.

# Arguments
- `state::SchwingerState`: Schwinger model state.
- `link::Int`: link.
"""
function electricfield(state::SchwingerState, link::Int)
    if isfinite(lattice(state).N) && !(1 ≤ link ≤ lattice(state).N)
        throw(ArgumentError("link must be between 1 and N"))
    end
    if isinf(lattice(state).N)
        link = link % 2 == 0 ? 2 : 1
    end
    return sum(charges(state)[1:link]) .+ L₀(state) .+ lattice(state).θ2π
end

"""
`electricfields(state)`

Return a list of the expectations of (L + θ/2π) operators on each link.

# Arguments
- `state::SchwingerState`: Schwinger model state.
"""
function electricfields(state::SchwingerState)
    lat = lattice(state)
    Qs = charges(state)

    N = isinf(lattice(state).N) ? 2 : Int(lattice(state).N)
    efs = accumulate(+, Qs) .+ L₀(state) .+ Base.convert(Vector{Float64}, lat.θ2π[1:N])
    # static defect charges are not in the (unaltered) lattice θ2π nor in `charges`;
    # add their flux contribution so the reported field on the real links is physical.
    # Uniform across backends (ED/ITensors keep the original lattice; MPSKit hides its
    # extra site), since `defects(state)` is the same regardless of representation.
    ds = defects(state)
    if !isempty(ds)
        efs = efs .+ _defect_theta_shift(lat, ds)
    end
    return efs
end

# On an infinite lattice there is no boundary to anchor Gauss's law, so the generic
# accumulate-charge formula above cannot resolve the background field — it returns the same
# value for the two θ=π vacua (and is unreliable in general). Measure ⟨L+θ⟩ on each link
# directly from the bond instead, via `link_expectation` (as the energy-density code does).
# This works for any MPS over an infinite lattice — the uniform `InfiniteMPS` vacua as well as
# a finite window (`WindowMPS`/`FiniteMPS`) cut from one, e.g. an evolving wavepacket.
function electricfields(state::MPSKitState)
    lat = lattice(state)
    isinf(lat.N) || return invoke(electricfields, Tuple{SchwingerState}, state)
    u  = state.hamiltonian.universe
    ψ  = state.psi
    nrm2 = ψ isa MPSKit.InfiniteMPS ? 1.0 : real(dot(ψ, ψ))   # windows can drift from norm 1
    return [real(link_expectation(ψ, ℓ,
                r::U1Irrep -> Float64(r.charge) + Float64(lat.θ2π[mod1(ℓ, length(lat.θ2π))]) + u) / nrm2)
            for ℓ in 1:length(ψ)]
end

function electricfield(state::MPSKitState, link::Int)
    lat = lattice(state)
    isinf(lat.N) || return invoke(electricfield, Tuple{SchwingerState,Int}, state, link)
    return electricfields(state)[mod1(link, length(state.psi))]
end

function partialcode(state::BasisState, range::UnitRange{Int})
    N, F = Int(lattice(state).N), lattice(state).F
    if range.start > N
        return BitVector(), state.L₀
    else
        fermioncode = reshape(occupations(state)[range,:], (length(range)*F,))
        if range.start > 1
            return fermioncode, state.L₀
        else
            return fermioncode, 0
        end
    end
end

"""
`entanglement(state, bisection)`

Return the von Neumann entanglement entropy -tr(ρₐ log(ρₐ)), where a is the subsystem of sites 1..bisection

# Arguments
- `state::EDState`: Schwinger model state.
- `bisection::Int`: bisection index.
"""
function entanglement(state::EDState, bisection::Int)
    N = Int(lattice(state).N)
    basis = _edbasis(state)
    range1 = bisection ≤ N/2 ? (1:bisection) : (bisection+1:N)
    range2 = bisection ≤ N/2 ? (bisection+1:N) : (1:bisection)
    codes = partialcode.(basis, Ref(range1)) # don't broadcast range1
    index = Dict([code => i for (i, code) in enumerate(Set(codes))])

    grouped_states = Dict{Tuple{BitVector, Int}, Vector{Int}}()
    for (i, basis_state) in enumerate(basis)
        key = partialcode(basis_state, range2)
        if haskey(grouped_states, key)
            push!(grouped_states[key], i)
        else
            grouped_states[key] = [i]
        end
    end

    partialtrace = zeros(ComplexF64, length(index), length(index))
    for group in values(grouped_states)
        for i in group
            for j in group
                partialtrace[index[codes[i]], index[codes[j]]] += state.coeffs[i] * conj(state.coeffs[j])
            end
        end
    end

    factorization = LinearAlgebra.eigen(partialtrace)
    return -real(sum(p == 0 ? 0 : p * log(abs(p)) for p in factorization.values))
end

"""
`entanglement(state, bisection)`

Return the von Neumann entanglement entropy -tr(ρₐ log(ρₐ)), where a is the subsystem of sites 1..bisection

# Arguments
- `state::ITensorState`: Schwinger model state.
- `bisection::Int`: bisection index.
"""
function entanglement(state::ITensorState, bisection::Int)
    F = lattice(state).F
    psi = state.psi
    psiorth = ITensorMPS.orthogonalize(psi, bisection*F)
    _,S,_ = svd(psiorth[bisection*F], (linkinds(psiorth, bisection*F-1)..., siteinds(psiorth,bisection*F)...))
    return -sum(p * log(p) for p in diag(S) .^ 2)
end

"""
`entanglement(state, bisection)`

Return the von Neumann entanglement entropy -tr(ρₐ log(ρₐ)), where a is the subsystem of sites 1..bisection

# Arguments
- `state::MPSKitState`: Schwinger model state.
- `bisection::Int`: bisection index.
"""
function entanglement(state::MPSKitState, bisection::Int)
    F = lattice(state).F
    psi = state.psi
    # Get the entanglement entropy at the bond after site bisection*F
    # MPSKit provides entanglement_entropy function
    bond_idx = bisection * F
    return MPSKit.entropy(psi, bond_idx)
end

"""
`entanglements(state)`

Return a list of the von Neumann entanglement entropies for each bisection of the lattice.

# Arguments
- `state::SchwingerState`: Schwinger model state.
"""
function entanglements(state::SchwingerState)
    N = isinf(lattice(state).N) ? 2 : Int(lattice(state).N)
    return [entanglement(state, i) for i=1:(N - (lattice(state).periodic || isinf(lattice(state).N) ? 0 : 1))]
end