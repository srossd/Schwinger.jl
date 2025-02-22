abstract type SchwingerState{N,F} end

"""
`SchwingerBasisState{N,F}(occupations, L0)`

A Schwinger model basis state.
"""
struct SchwingerBasisState{N,F} <: SchwingerState{N,F}
    lattice::SchwingerLattice{N,F}
    occupations::BitMatrix
    L₀::Int

    function SchwingerBasisState(lattice::SchwingerLattice{N,F}, occupations::BitMatrix, L₀::Int) where {N,F}
        if size(occupations) != (N,F)
            throw(ArgumentError("occupations must be an NxF BitMatrix"))
        end
        new{N,F}(lattice, occupations, L₀)
    end
end

function nstates(N::Int, F::Int; L_max::Int = 3)
    return (2*L_max+1)*binomial(N*F, N*F÷2)
end

function lattice(state::SchwingerBasisState{N,F}) where {N,F}
    return state.lattice
end

"""
`basis(N, F; L_max)`

Returns a basis of Schwinger model states.

# Arguments
- `N::Int`: number of sites.
- `F::Int`: number of flavors.
- `L_max::Int`: maximum value of |L₀|.
"""

@memoize function basis(lattice::SchwingerLattice{N,F}; L_max::Union{Nothing,Int} = nothing) where {N,F}
    if !isnothing(L_max) && L_max < 0
        throw(ArgumentError("L_max must be non-negative"))
    end
    if !isnothing(L_max) && L_max > 0 && !lattice.periodic
        throw(ArgumentError("L_max must be 0 for open boundary conditions"))
    end

    L_max = isnothing(L_max) ? (lattice.periodic ? 3 : 0) : L_max
    states = Vector{SchwingerBasisState{N,F}}(undef, nstates(N, F; L_max = L_max))
    stateidx = 1
    for L₀ in -L_max:L_max
        for comb in combinations(1:N*F, N*F÷2)
            occupations = BitMatrix(undef, N, F)
            for idx in comb
                occupations[idx] = true
            end
            states[stateidx] = SchwingerBasisState(lattice, occupations, L₀)
            stateidx += 1
        end
    end
    return states
end

"""
`SchwingerEDState{N,F}(hamiltonian, coeffs)`

A Schwinger model state represented as a linear combination of basis states.
"""
struct SchwingerEDState{N,F} <: SchwingerState{N,F}
    hamiltonian::SchwingerHamiltonian{N,F}
    coeffs::Vector{ComplexF64}

    function SchwingerEDState(hamiltonian::SchwingerHamiltonian{N,F}, coeffs::Vector{ComplexF64}) where {N,F}
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
    hamiltonian::SchwingerHamiltonian{N,F}
    psi::MPS

    function SchwingerMPS(hamiltonian::SchwingerHamiltonian{N,F}, psi::MPS) where {N,F}
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
- `hamiltonian::MPOHamiltonian`: Schwinger model Hamiltonian.
- `nstates::Int`:: number of states to determine.
"""

function loweststates(hamiltonian::MPOHamiltonian{N,F}, nstates::Int; 
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

function loweststates(hamiltonian::EDHamiltonian, nstates::Int; shift=-100)
    vals, vecs = if nstates + 2 < size(hamiltonian.matrix)[1]
         Arpack.eigs(hamiltonian.matrix; nev=nstates, which=:LM, sigma = shift) # LM gives eigenvalues closest to sigma
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
- `hamiltonian::SchwingerHamiltonian`: Schwinger model Hamiltonian.
"""

function groundstate(hamiltonian::SchwingerHamiltonian{N,F}; kwargs...) where {N,F}
    return loweststates(hamiltonian, 1; kwargs...)[1]
end

"""
`energygap(hamiltonian)`
Returns the energy difference between the lowest two states of the Hamiltonian.

# Arguments
- `hamiltonian::SchwingerHamiltonian`: Schwinger model Hamiltonian.
"""

function energygap(hamiltonian::SchwingerHamiltonian{N,F}; kwargs...) where {N,F}
    return abs(-(map(energy, loweststates(hamiltonian, 2; kwargs...))...))
end

"""
`energy(state)`

Return the expectation value of the Hamiltonian.

# Arguments
- `state::SchwingerMPS`: Schwinger model state.
"""
function energy(state::SchwingerMPS{N,F}) where {N,F}
    H = state.hamiltonian.mpo
    psi = state.psi
    return inner(psi', H, psi)
end

"""
`energy(state)`

Return the expectation value of the Hamiltonian.

# Arguments
- `state::SchwingerEDState`: Schwinger model state.
"""
function energy(state::SchwingerEDState{N,F}) where {N,F}
    H = state.hamiltonian.matrix
    coeffs = state.coeffs
    return real(dot(coeffs, H*coeffs))
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
`charges(state)`

Return a list of the expectations of Q operators on each site and for each known eigenstate.

# Arguments
- `state::SchwingerState`: Schwinger model state.
"""
function charges(state::SchwingerState{N,F}) where {N,F}
    return sum(occupations(state), dims=2) + (F .* [-1/2 + (-1)^(n - 1)/2 for n=1:N])
end

function L₀(state::SchwingerBasisState{N,F}) where {N,F}
    return lattice(state).periodic ? state.L₀ : 0
end

function L₀(state::SchwingerMPS{N,F}) where {N,F}
    return lattice(state).periodic ? expect(psi, "L0", sites=N*F + 1) : 0
end

function L₀(state::SchwingerEDState{N,F}) where {N,F}
    states = basis(lattice(state); L_max = state.hamiltonian.L_max)
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