abstract type SchwingerState{N,F} end

"""
`SchwingerBasisState{N,F}(occupations, L0)`

A Schwinger model basis state.
"""
struct SchwingerBasisState{N,F} <: SchwingerState{N,F}
    lattice::SchwingerLattice{N,F}
    occupations::BitMatrix
    L₀::Int
    q::Int

    function SchwingerBasisState(lattice::SchwingerLattice{N,F}, occupations::BitMatrix, L₀::Int, q::Int) where {N,F}
        if size(occupations) != (N,F)
            throw(ArgumentError("occupations must be an NxF BitMatrix"))
        end
        if q < 1
            throw(ArgumentError("q must be ≥ 1"))
        end
        new{N,F}(lattice, occupations, L₀, q)
    end
end

function nstates(N::Int, F::Int; L_max::Int = 3)
    return (2*L_max+1)*binomial(N*F, N*F÷2)
end

function lattice(state::SchwingerBasisState{N,F}) where {N,F}
    return state.lattice
end

"""
`basis(lattice; L_max, q, universe)`

Returns a basis of Schwinger model states.

# Arguments
- `lattice::SchwingerLattice`
- `L_max::Int`
- `q::Int`
- `universe::Int`
"""
@memoize function basis(lattice::SchwingerLattice{N,F}; L_max::Union{Nothing,Int} = nothing, q::Int = 1, universe::Int = 0) where {N,F}
    if !isnothing(L_max) && L_max < 0
        throw(ArgumentError("L_max must be non-negative"))
    end
    if !isnothing(L_max) && L_max > 0 && !lattice.periodic
        throw(ArgumentError("L_max must be 0 for open boundary conditions"))
    end
    if q < 1
        throw(ArgumentError("q must be ≥ 1"))
    end

    L_max = isnothing(L_max) ? (lattice.periodic ? 3 : 0) : L_max
    states = Vector{SchwingerBasisState{N,F}}(undef, nstates(N, F; L_max = L_max))
    stateidx = 1
    for L₀ in ((-L_max:L_max) .* q) .+ universe
        for comb in combinations(1:N*F, N*F÷2)
            occupations = BitMatrix(undef, N, F)
            for idx in comb
                occupations[idx] = true
            end
            states[stateidx] = SchwingerBasisState(lattice, occupations, L₀, q)
            stateidx += 1
        end
    end
    return states
end

@memoize function positionindex(lattice::SchwingerLattice{N,F}; L_max::Union{Nothing,Int} = nothing, q::Int = 1, universe::Int = 0) where {N,F}
    states = basis(lattice; L_max = L_max, q = lattice.q, universe = universe)
    return Dict{Tuple{BitMatrix,Int},Int}((states[i].occupations, states[i].L₀) => i for i in eachindex(states))
end

"""
`SchwingerEDState{N,F}(hamiltonian, coeffs)`

A Schwinger model state represented as a linear combination of basis states.
"""
struct SchwingerEDState{N,F} <: SchwingerState{N,F}
    hamiltonian::EDOperator{N,F}
    coeffs::Vector{ComplexF64}

    function SchwingerEDState(hamiltonian::EDOperator{N,F}, coeffs::Vector{ComplexF64}) where {N,F}
        new{N,F}(hamiltonian, coeffs)
    end
end

function lattice(state::SchwingerEDState{N,F}) where {N,F}
    return state.hamiltonian.lattice
end

"""
`SchwingerMPSState{N,F}(hamiltonian, psi)`

A Schwinger model MPS.
"""
struct SchwingerMPS{N,F} <: SchwingerState{N,F}
    hamiltonian::MPOOperator{N,F}
    psi::MPS

    function SchwingerMPS(hamiltonian::MPOOperator{N,F}, psi::MPS) where {N,F}
        new{N,F}(hamiltonian, psi)
    end
end

function lattice(state::SchwingerMPS{N,F}) where {N,F}
    return state.hamiltonian.lattice
end

"""
`lowest_states(hamiltonian, nstates)`

Returns the lowest few eigenstates of the Schwinger model Hamiltonian.

# Arguments
- `hamiltonian::MPOOperator`: Schwinger model Hamiltonian.
- `nstates::Int`:: number of states to determine.
"""
function loweststates(hamiltonian::MPOOperator{N,F}, nstates::Int; 
    maxiters::Int = 500, initiallinkdim::Int = 4, maxlinkdim::Int = 600, energy_tol::Real = 1E-6, weight::Real = 100, outputlevel::Int = 0, minsweeps::Int = 5) where {N,F}

    H = hamiltonian.mpo

    states = Vector{SchwingerMPS}(undef, nstates)
    for idx in 1:nstates
        state = [n == N * F + 1 ? hamiltonian.L_max + 1 : isodd(floor((n-1)/F)) ? "Up" : "Dn" for n=1:(N * F + (hamiltonian.lattice.periodic ? 1 : 0))]
        psi = random_mps(sites(hamiltonian.lattice; L_max = hamiltonian.L_max), state; linkdims = initiallinkdim)

        sweeps = Sweeps(maxiters)
        maxdim!(sweeps, maxlinkdim)

        if outputlevel >= 1
            println("Running DMRG for state $idx")
        end
        _, psi = ITensorMPS.dmrg(H, [state.psi for state in states[1:idx-1]], psi, sweeps; weight = weight, observer = DMRGObserver(;energy_tol = energy_tol,minsweeps=minsweeps), outputlevel = outputlevel)
        if outputlevel >= 1
            println()
        end

        states[idx] = SchwingerMPS(hamiltonian, psi)
    end

    return states
end

function loweststates(hamiltonian::EDOperator, nstates::Int; ncv::Int = max(20, 2*nstates + 1))
    _, vecs = if nstates + 2 < size(hamiltonian.matrix)[1]
        try
            Arpack.eigs(hamiltonian.matrix; nev=nstates, ncv = ncv, which=:SR)
        catch e
            if e isa Arpack.XYAUPD_Exception
                @warn "Arpack.eigs failed to converge, increasing ncv to $(2*ncv)"
                return loweststates(hamiltonian, nstates; ncv = 2*ncv)
            else
                rethrow(e)
            end
        end
    else
        eigs = eigen(Matrix(hamiltonian.matrix))
        (eigs.values[1:nstates], eigs.vectors[:,1:nstates])
    end
    return [SchwingerEDState(hamiltonian, vecs[:,n]) for n=1:nstates]
end

"""
`groundstate(hamiltonian)`
Returns the ground state of the Schwinger model Hamiltonian.

# Arguments
- `hamiltonian::SchwingerOperator`: Schwinger model Hamiltonian.
"""
function groundstate(hamiltonian::SchwingerOperator{N,F}; kwargs...) where {N,F}
    return loweststates(hamiltonian, 1; kwargs...)[1]
end

"""
`energygap(hamiltonian)`
Returns the energy difference between the lowest two states of the Hamiltonian.

# Arguments
- `hamiltonian::SchwingerOperator`: Schwinger model Hamiltonian.
"""
function energygap(hamiltonian::SchwingerOperator{N,F}; kwargs...) where {N,F}
    return abs(-(map(energy, loweststates(hamiltonian, 2; kwargs...))...))
end

"""
`expectation(op, state)`

Return the expectation value of the operator `op` in `state`.

# Arguments
- `op::EDOperator{N,F}``: operator.
- `state::SchwingerEDState`: state.
"""
function expectation(op::EDOperator{N,F}, state::SchwingerEDState{N,F}) where {N,F}
    return dot(state.coeffs, op.matrix * state.coeffs)
end

"""
`expectation(op, state)`

Return the expectation value of the operator `op` in `state`.

# Arguments
- `op::MPOOperator{N,F}``: operator.
- `state::SchwingerMPS`: state.
"""
function expectation(op::MPOOperator{N,F}, state::SchwingerMPS{N,F}) where {N,F}
    return inner(state.psi', op.mpo, state.psi)
end

"""
`energy(state)`

Return the expectation value of the Hamiltonian.

# Arguments
- `state::SchwingerEDState`: Schwinger model state.
"""
function energy(state::Union{SchwingerEDState{N,F},SchwingerMPS{N,F}}) where {N,F}
    return real(expectation(state.hamiltonian, state))
end

"""
`energy(state)`

Return the expectation value of the Hamiltonian.

# Arguments
- `state::SchwingerBasisState`: Schwinger model basis state.
"""
function energy(state::SchwingerBasisState{N,F}) where {N,F}
    efs = electricfields(state) # includes θ/2π
    occs = occupations(state)

    electricenergy = (state.lattice.a/2) * sum(efs.^2)
    massenergy = sum((-1)^(j-1) * state.lattice.mlat[j][k] * occs[j,k] for j=1:N, k=1:F)

    return electricenergy + massenergy
end

"""
`occupations(state)`

Return an NxF matrix of the expectations of χ†χ operators on each site.

# Arguments
- `state::SchwingerMPS{N,F}`: Schwinger model state.
"""
function occupations(state::SchwingerMPS{N,F}) where {N,F}
    psi = state.psi
    return transpose(reshape(expect(psi, "Sz", sites=1:N*F) + [1/2 for n=1:N*F], (F,N)))
end

"""
`occupations(state)`

Return an NxF matrix of the expectations of χ†χ operators on each site.

# Arguments
- `state::SchwingerBasisState`: Schwinger model basis state.
"""
function occupations(state::SchwingerBasisState{N,F}) where {N,F}
    return state.occupations
end

"""
`occupations(state)`

Return an NxF matrix of the expectations of χ†χ operators on each site.

# Arguments
- `state::SchwingerBasisState`: Schwinger model basis state.
"""
function occupations(state::SchwingerEDState{N,F}) where {N,F}
    states = basis(state.hamiltonian.lattice; L_max = state.hamiltonian.L_max)
    occs = zeros(N,F)
    for (coeff, state) in zip(state.coeffs, states)
        occs .+= abs2(coeff) .* occupations(state)
    end
    return occs
end

"""
`scalarvev(state)`

Return the VEV of the scalar condensate L⁻¹ ∑ (-1)ⁿ χ†ₙχₙ

# Arguments
- `state::SchwingerState`: Schwinger model state.
"""
function scalarvev(state::SchwingerState{N,F}) where {N,F}
    occs = occupations(state)[:,1]
    return (repeat([1,-1],N÷2)'occs)/lattice(state).L
end

"""
`pseudoscalarvev(state)`

Return the VEV of the pseudoscalar condensate L⁻¹ ∑ (-1)ⁿ χ†ₙχₙ

# Arguments
- `state::SchwingerState`: Schwinger model state.
"""
function pseudoscalarvev(state::SchwingerMPS{N,F}) where {N,F}
    return real(expectation(MPOHoppingMass(lattice(state); bare = true), state)/lattice(state).L)
end

function pseudoscalarvev(state::SchwingerEDState{N,F}) where {N,F}
    return real(expectation(EDHoppingMass(lattice(state); bare = true), state)/lattice(state).L)
end

function pseudoscalarvev(state::SchwingerBasisState{N,F}) where {N,F}
    return 0
end

"""
`charges(state)`

Return a list of the expectations of Q operators on each site and for each known eigenstate.

# Arguments
- `state::SchwingerState`: Schwinger model state.
"""
function charges(state::SchwingerState{N,F}) where {N,F}
    return (sum(occupations(state), dims=2) + (F .* [-1/2 + (-1)^(n - 1)/2 for n=1:N])) .* lattice(state).q
end

"""
`L₀(state)`

Return the expectation value of L₀.

# Arguments
- `state::SchwingerState`: Schwinger model state.
"""
function L₀(state::SchwingerBasisState{N,F}) where {N,F}
    return lattice(state).periodic ? state.L₀ : 0
end

function L₀(state::SchwingerMPS{N,F}) where {N,F}
    return lattice(state).periodic ? lattice(state).q * expect(state.psi, "L0", sites=N*F + 1) : 0
end

function L₀(state::SchwingerEDState{N,F}) where {N,F}
    states = basis(lattice(state); L_max = state.hamiltonian.L_max, q = lattice(state).q)
    return sum(abs2(coeff) * L₀(state) for (coeff, state) in zip(state.coeffs, states))
end

"""
`electricfields(state)`

Return a list of the expectations of (L + θ/2π) operators on each link.

# Arguments
- `state::SchwingerMPS`: Schwinger model state.
"""
function electricfields(state::SchwingerState{N,F}) where {N,F}
    lat = lattice(state)
    Qs = charges(state)

    return circshift(accumulate(+, Qs) .+ L₀(state),1) .+ lat.θ2π
end

"""
`entanglements(state)`

Return a list of the entanglement entropies for each bisection of the lattice.

# Arguments
- `state::SchwingerMPS`: Schwinger model state.
"""
function entanglements(state::SchwingerMPS{N,F}) where {N,F}
    data = Vector{Float64}(undef, N*F - (lattice.periodic ? 0 : 1))
    psi = state.psi

    for j=1:N*F - (lattice.periodic ? 0 : 1)
        orthogonalize!(psi, j)
        _,S,_ = svd(psi[j], (linkinds(psi, j-1)..., siteinds(psi,j)...))
        for n=1:dim(S, 1)
            p = S[n,n]^2
            data[j] -= p * log(p)
        end
    end

    return data
end