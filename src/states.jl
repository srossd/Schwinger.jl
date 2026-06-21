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
"""
function loweststates(hamiltonian::ITensorOperator, nstates::Int; 
    maxiters::Int = 500, initiallinkdim::Int = 4, maxlinkdim::Int = 600, energy_tol::Real = 1E-6, weight::Real = 100, outputlevel::Int = 0, minsweeps::Int = 5)

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
function loweststates(hamiltonian::EDOperator, nstates::Int; ncv::Int = max(20, 2*nstates + 1))
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
`lowest_states(hamiltonian, nstates)`

Returns the lowest few eigenstates of the Schwinger model Hamiltonian using MPSKit.

# Arguments
- `hamiltonian::MPSKitOperator`: Schwinger model Hamiltonian.
- `nstates::Int`:: number of states to determine.
"""
function loweststates(hamiltonian::MPSKitOperator, nstates::Int;
    maxiters::Int = 500, initiallinkdim::Int = 20, bonddim::Union{Nothing,Int} = nothing, initial_Lmax::Int = 3, energy_tol::Real = 1E-6, cutoff::Real = 1E-10, weight::Real = 100., verbose::Bool = false, momentum::Union{Real, Nothing} = nothing)

    initiallinkdim = something(bonddim, initiallinkdim)   # `bonddim` is an alias

    H = hamiltonian.lempo
    spaces = get_mpskit_spaces(hamiltonian.lattice)   # defect sites are added in the finite branch
    total_defect = sum(d.charge for d in hamiltonian.defects; init = 0)
    states = Vector{Union{MPSKitState,MPSKitQPState}}(undef, nstates)
    if isinf(hamiltonian.lattice.N)
        ψ₀ = MPSKit.InfiniteMPS(spaces, [U1Space([q => initiallinkdim for q in hamiltonian.lattice.q*(-initial_Lmax:initial_Lmax)]) for _ in 1:length(spaces)]) #TODO: add attenuation near theta = pi
        if abs(hamiltonian.lattice.θ2π[1] - 0.5) < 0.1
            ψ₀ = attenuateLinks(ψ₀, hamiltonian.lattice.θ2π[1] < 0.5 ? [U1Irrep(0), U1Irrep(0)] : [U1Irrep(-1), U1Irrep(-1)], 0.01)
        end
        alg = MPSKit.VUMPS(; maxiter=maxiters, tol=energy_tol, verbosity=verbose ? 1 : 0) #TODO: use VUMPSSVDCut once upstream bug fixed
        ψ, envs, _ = MPSKit.find_groundstate(ψ₀, H, alg)
        states[1] = MPSKitState(hamiltonian, ψ)
        
        if nstates > 1
            isnothing(momentum) && (momentum = 0.0)
            algqp = MPSKit.QuasiparticleAnsatz()
            _, psis = MPSKit.excitations(H, algqp, momentum, ψ, envs; num = nstates - 1)
            for n in 2:nstates
                states[n] = MPSKitQPState(hamiltonian, psis[n-1])
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