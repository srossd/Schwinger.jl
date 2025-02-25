abstract type SchwingerOperator{N,F} end

struct EDOperator{N,F} <: SchwingerOperator{N,F}
    lattice::SchwingerLattice{N,F}
    matrix::SparseMatrixCSC{ComplexF64,Int64}
    L_max::Int
    universe::Int

    function EDOperator(lattice::SchwingerLattice{N,F}, matrix::SparseMatrixCSC{ComplexF64,Int64}, L_max::Int, universe::Int) where {N,F}
        return new{N,F}(lattice, matrix, L_max, universe)
    end
end

struct MPOOperator{N,F} <: SchwingerOperator{N,F}
    lattice::SchwingerLattice{N,F}
    mpo::MPO
    L_max::Int
    universe::Int

    function MPOOperator(lattice::SchwingerLattice{N,F}, mpo::MPO, L_max::Int, universe::Int) where {N,F}
        return new{N,F}(lattice, mpo, L_max, universe)
    end
end

abstract type SchwingerHamiltonian{N,F} <: SchwingerOperator{N,F} end

struct EDHamiltonian{N,F} <: SchwingerHamiltonian{N,F}
    lattice::SchwingerLattice{N,F}
    matrix::SparseMatrixCSC{ComplexF64,Int64}
    L_max::Int64
    universe::Int64

    function EDHamiltonian(lattice::SchwingerLattice{N,F}, matrix::SparseMatrixCSC{ComplexF64,Int64}, L_max::Int, universe::Int64) where {N,F}
        new{N,F}(lattice, matrix, L_max, universe)
    end
end

struct MPOHamiltonian{N,F} <: SchwingerHamiltonian{N,F}
    lattice::SchwingerLattice{N,F}
    mpo::MPO
    L_max::Int64
    universe::Int64

    function MPOHamiltonian(lattice::SchwingerLattice{N,F}, mpo::MPO, L_max::Int, universe::Int) where {N,F}
        new{N,F}(lattice, mpo, L_max, universe)
    end
end

"""
`wilsonloop(hamiltonian, conjugate = false)`

Returns the spatial Wilson loop operator for `lattice`.

# Arguments
- `lattice::SchwingerLattice`: lattice.
- `conjugate::Bool`: Conjugation of the Wilson loop.
"""
function MPOWilsonLoop(lattice::SchwingerLattice{N,F}, conjugate::Bool = false; L_max::Union{Nothing,Int} = nothing, universe::Int = 0) where {N,F}
    if !lattice.periodic
        throw(DomainError(lattice, "Lattice must be periodic."))
    end

    if !isnothing(L_max) && L_max < 0
        throw(ArgumentError("L_max must be non-negative"))
    end
    if !isnothing(L_max) && L_max > 0 && !lattice.periodic
        throw(ArgumentError("L_max must be 0 for open boundary conditions"))
    end

    universe = mod(universe, lattice.q)
    if universe < 0
        universe += q
    end

    L_max = isnothing(L_max) ? (lattice.periodic ? 3 : 0) : L_max

    holonomy = OpSum()
    holonomy += (conjugate ? "lower" : "raise"),N * F + 1
    
    mpo = MPO(holonomy, sites(hamiltonian))
    return MPOOperator(hamiltonian, mpo, L_max, universe)
end

"""
`wilsonloop(hamiltonian, conjugate = false)`

Returns the spatial Wilson loop operator for `lattice`.

# Arguments
- `lattice::SchwingerLattice`: lattice.
- `conjugate::Bool`: Conjugation of the Wilson loop.
"""
function EDWilsonLoop(lattice::SchwingerLattice{N,F}, conjugate::Bool = false; L_max::Union{Nothing,Int} = nothing, universe::Int = 0) where {N,F}
    if !lattice.periodic
        throw(DomainError(lattice, "Lattice must be periodic."))
    end

    if !isnothing(L_max) && L_max < 0
        throw(ArgumentError("L_max must be non-negative"))
    end
    if !isnothing(L_max) && L_max > 0 && !lattice.periodic
        throw(ArgumentError("L_max must be 0 for open boundary conditions"))
    end

    universe = mod(universe, lattice.q)
    if universe < 0
        universe += q
    end

    L_max = isnothing(L_max) ? (lattice.periodic ? 3 : 0) : L_max

    states = basis(lattice; L_max = L_max, q = lattice.q, universe = universe)
    positions = Dict{Tuple{BitMatrix,Int},Int}((states[i].occupations, states[i].L₀) => i for i in eachindex(states))

    I = Vector{Int}(undef, length(states))
    J = Vector{Int}(undef, length(states))
    V = Vector{ComplexF64}(undef, length(states))
    idx = 1
    for i in eachindex(states)
        key = (states[i].occupations, states[i].L₀ + (conjugate ? -1 : 1))
        if haskey(positions, key)
            I[idx] = i
            J[idx] = positions[key]
            V[idx] = 1
            idx += 1
        end
    end

    matrix = sparse(I[1:idx-1], J[1:idx-1], V[1:idx-1], length(states), length(states))
    return EDOperator(lattice, matrix, L_max, universe)
end