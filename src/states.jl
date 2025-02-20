"""
`SchwingerState{N,F}(lattice, psi, energy)`

A Schwinger model state.
"""

struct SchwingerState{N,F}
    lattice::SchwingerLattice{N,F}
    psi::MPS

    function SchwingerState(lattice::SchwingerLattice{N,F}, psi::MPS) where {N,F}
        new{N,F}(lattice, psi)
    end
end

"""
`lowest_states(lattice, nstates)`

Returns the lowest few eigenstates of the Schwinger model Hamiltonian on the given lattice.

# Arguments
- `lattice::SchwingerLattice`: Schwinger model lattice.
- `nstates::Int`:: number of states to determine.
"""

function loweststates(lattice::SchwingerLattice{N,F}, nstates::Int; 
    L_max::Int, maxiters::Int = 500, initiallinkdim::Int = 4, maxlinkdim::Int = 600, energy_tol::Real = 1E-6, weight::Real = 100, outputlevel::Int = 1, minsweeps::Int = 10) where {N,F}

    H = hamiltonian(lattice; L_max = L_max)

    states = Vector{SchwingerState}(undef, nstates)
    for idx in 1:nstates
        state = [n == N * F + 1 ? L_max + 1 : isodd(floor((n-1)/F)) ? "Up" : "Dn" for n=1:(N * F + (lattice.periodic ? 1 : 0))]
        psi = random_mps(sites(lattice), state; linkdims = initiallinkdim)

        sweeps = Sweeps(maxiters)
        maxdim!(sweeps, maxlinkdim)

        if outputlevel >= 1
            println("Running DMRG for state $idx")
        end
        _, psi = ITensorMPS.dmrg(H, [state.psi for state in states[1:idx-1]], psi, sweeps; weight = weight, observer = DMRGObserver(;energy_tol = energy_tol,minsweeps=minsweeps), outputlevel = outputlevel)
        if outputlevel >= 1
            println()
        end

        states[idx] = SchwingerState(lattice, psi)
    end

    return states
end

"""
`groundstate(lattice)`
Returns the ground state of the Schwinger model Hamiltonian on the given lattice.

# Arguments
- `lattice::SchwingerLattice`: Schwinger model lattice.
"""

function groundstate(lattice::SchwingerLattice{N,F}; kwargs...) where {N,F}
    return loweststates(lattice, 1; kwargs...)[1]
end

"""
`energy(state)`

Return the Hamiltonian eigenvalues of all known eigenstates.

# Arguments
- `state::SchwingerState`: Schwinger model state.
"""
function energy(state::SchwingerState{N,F}) where {N,F}
    H = hamiltonian(state.lattice)
    psi = state.psi
    return inner(psi', H, psi)
end

"""
`occupations(state)`

Return a list of the expectations of χ†χ operators on each site.

# Arguments
- `state::SchwingerState`: Schwinger model state.
"""
function occupations(state::SchwingerState{N,F}) where {N,F}
    psi = state.psi
    return expect(psi, "Sz", sites=1:N*F) + [1/2 for n=1:N*F]
end

"""
`charges(state)`

Return a list of the expectations of Q operators on each site and for each known eigenstate.

# Arguments
- `state::SchwingerState`: Schwinger model state.
"""
function charges(state::SchwingerState{N,F}) where {N,F}
    occs = occupations(state)

    return occs + [-1/2 + (-1)^((n - 1) ÷ F)/2 for n=1:N*F]
end

"""
`electricfields(state)`

Return a list of the expectations of (L + θ/2π) operators on each link.

# Arguments
- `state::SchwingerState`: Schwinger model state.
"""
function electricfields(state::SchwingerState{N,F}) where {N,F}
    psi = state.psi
    lattice = state.lattice
    Qs = charges(state)

    L0 = lattice.periodic ? expect(psi, "L0", sites=N*F + 1) : 0

    return accumulate(+, Qs) .+ L0 .+ lattice.θ2π
end

"""
`entanglements(state)`

Return a list of the entanglement entropies for each bisection of the lattice.

# Arguments
- `state::SchwingerState`: Schwinger model state.
"""
function entanglements(state::SchwingerState{N,F}) where {N,F}
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