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
        if sum(occupations) != (lattice.N * lattice.F) ÷ 2
            throw(ArgumentError("occupations must have N*F/2 occupied sites"))
        end
        if q < 1
            throw(ArgumentError("q must be ≥ 1"))
        end
        new(lattice, occupations, L₀, q)
    end
end

function nstates(N::Int, F::Int; L_max::Int = 3)
    return (2*L_max+1)*binomial(N*F, N*F÷2)
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
@memoize function schwingerbasis(lattice::Lattice; L_max::Union{Nothing,Int} = nothing, universe::Int = 0)
    N, F = Int(lattice.N), lattice.F
    if !isnothing(L_max) && L_max < 0
        throw(ArgumentError("L_max must be non-negative"))
    end
    if !isnothing(L_max) && L_max > 0 && !lattice.periodic
        throw(ArgumentError("L_max must be 0 for open boundary conditions"))
    end

    L_max = isnothing(L_max) ? (lattice.periodic ? 3 : 0) : L_max

    states = Vector{BasisState}(undef, nstates(Int(N), F; L_max = L_max))
    stateidx = 1
    for L₀ in ((-L_max:L_max) .* lattice.q) .+ universe
        for comb in combinations(1:N*F, N*F÷2)
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

@memoize function positionindex(lattice::Lattice; L_max::Union{Nothing,Int} = nothing, universe::Int = 0)
    states = schwingerbasis(lattice; L_max = L_max, universe = universe)
    return Dict{Tuple{BitMatrix,Int},Int}((states[i].occupations, states[i].L₀) => i for i in eachindex(states))
end

"""
`EDState(hamiltonian, coeffs)`

A Schwinger model state represented as a linear combination of basis states.
"""
struct EDState <: SchwingerState
    hamiltonian::EDOperator
    coeffs::Vector{ComplexF64}

    function EDState(hamiltonian::EDOperator, coeffs::Vector{ComplexF64})
        new(hamiltonian, coeffs)
    end
end

function lattice(state::EDState)
    return state.hamiltonian.lattice
end

function Base.:*(scalar::Number, state::EDState)
    return EDState(state.hamiltonian, scalar * state.coeffs)
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

    function ITensorState(hamiltonian::ITensorOperator, psi::ITensorMPS.MPS)
        new(hamiltonian, psi)
    end
end

function lattice(state::ITensorState)
    return state.hamiltonian.lattice
end

function Base.:*(scalar::Number, state::ITensorState)
    return ITensorState(state.hamiltonian, scalar * state.psi)
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

    function MPSKitState(hamiltonian::MPSKitOperator, psi::MPSKit.AbstractMPS)
        isinf(hamiltonian.lattice) == (psi isa InfiniteMPS) || throw(ArgumentError("Hamiltonian and state must both be finite or infinite"))
        (psi isa InfiniteMPS) || (length(psi) == Int(hamiltonian.lattice.N) * hamiltonian.lattice.F) || throw(ArgumentError("State length does not match Hamiltonian lattice size"))
        new(hamiltonian, psi)
    end
end

function lattice(state::MPSKitState)
    return state.hamiltonian.lattice
end

function Base.:*(scalar::Number, state::MPSKitState)
    return MPSKitState(state.hamiltonian, scalar * state.psi)
end

function Base.:*(state::MPSKitState, scalar::Number)
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
    maxiters::Int = 500, bonddim::Int = 100, energy_tol::Real = 1E-8, weight::Real = 100., verbose::Bool = false)

    H = hamiltonian.lempo
    spaces = get_mpskit_spaces(hamiltonian.lattice)
    states = Vector{MPSKitState}(undef, nstates)
    if isinf(hamiltonian.lattice.N)
        nstates == 1 || throw(ArgumentError("Only ground state (nstates=1) supported for infinite lattices"))
        ψ₀ = MPSKit.InfiniteMPS(spaces, [U1Space([q => bonddim for q in hamiltonian.lattice.q*(-3:3)]) for _ in 1:length(spaces)]) #TODO: add attenuation near theta = pi
        alg = MPSKit.VUMPS(; maxiter=maxiters, tol=energy_tol, verbosity=verbose ? 1 : 0) #TODO: use VUMPSSVDCut once upstream bug fixed
        ψ, = MPSKit.find_groundstate(ψ₀, H, alg)
        states[1] = MPSKitState(hamiltonian, ψ)
        #TODO: QuasiparticleAnsatz for excited states in infinite case
    else
        ψ₀ = MPSKit.FiniteMPS(rand, ComplexF64, spaces, U1Space([q => bonddim for q in hamiltonian.lattice.q*(-4:4)]))
        alg = MPSKit.DMRG(; maxiter=maxiters, verbosity=verbose ? 1 : 0) #TODO: replace with DMRG2 once upstream bug fixed
        ψ, = MPSKit.find_groundstate(ψ₀, H, alg)
        states[1] = MPSKitState(hamiltonian, ψ)
        if nstates > 1
            alg2 = MPSKit.FiniteExcited(alg, weight)
            _, psis = MPSKit.excitations(H, alg2, (states[1].psi,); num = nstates - 1)
            for n in 2:nstates
                states[n] = MPSKitState(hamiltonian, psis[n-1])
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
    return ITensorState(state.hamiltonian, apply(op.mpo, state.psi))
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
    return EDState(state.hamiltonian, op.matrix * state.coeffs)
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
    return MPSKitState(state.hamiltonian, op.lempo * state.psi)
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
    if isinf(N) && !(1 ≤ site ≤ 2)
        return occupation(state, site % 2 == 0 ? 2 : 1)
    end
    if isfinite(N) && !(1 ≤ site ≤ Int(N))
        throw(ArgumentError("Site must be between 1 and N"))
    end
    psi = state.psi
    
    space = physicalspace(psi, (site - 1) * F + 1)
    occop = zeros(space ← space)
    if isodd(site)
        block(occop, U1Irrep(0)) .= 1.0
    else
        block(occop, U1Irrep(lattice(state).q)) .= 1.0
    end

    # MPSKit expectation values for Sz operator on each flavor at the site
    occs = zeros(F)
    for k in 1:F
        idx = (site - 1) * F + k
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
    N, F = state.hamiltonian.lattice.N, state.hamiltonian.lattice.F
    N = isinf(N) ? 2 : Int(N)
    psi = state.psi
    occs = zeros(N, F)
    
    oddspace = physicalspace(psi, 1)
    evenspace = physicalspace(psi, F + 1)
    occopodd = zeros(oddspace ← oddspace)
    occopeven = zeros(evenspace ← evenspace)
    block(occopodd, U1Irrep(0)) .= 1.0
    block(occopeven, U1Irrep(lattice(state).q)) .= 1.0

    for j in 1:N
        for k in 1:F
            idx = (j - 1) * F + k
            occs[j, k] = real(MPSKit.expectation_value(psi, idx => (iseven(j) ? occopeven : occopodd)))
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
    states = schwingerbasis(state.hamiltonian.lattice; L_max = state.hamiltonian.L_max)
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
    states = schwingerbasis(state.hamiltonian.lattice; L_max = state.hamiltonian.L_max)
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
        return real(expectation(EDHoppingMass(lattice(state), n; bare = true, L_max = state.hamiltonian.L_max, universe = state.hamiltonian.universe), state))
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
    N = isinf(lattice(state).N) ? 2 : Int(lattice(state).N)
    F = lattice(state).F
    return (sum(occupations(state), dims=2) + (F .* [-1/2 + (-1)^(n)/2 for n=1:N])) .* lattice(state).q
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
    states = schwingerbasis(lattice(state); L_max = state.hamiltonian.L_max, universe = state.hamiltonian.universe)
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
    return accumulate(+, Qs) .+ L₀(state) .+ Base.convert(Vector{Float64},  lat.θ2π[1:N])
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
    basis = schwingerbasis(lattice(state); L_max = state.hamiltonian.L_max)
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