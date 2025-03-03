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

function Base.:+(op1::EDOperator{N,F}, op2::EDOperator{N,F}) where {N,F}
    if op1.lattice != op2.lattice
        return ArgumentError("Cannot add operators on different lattices")
    end
    if op1.L_max != op2.L_max
        return ArgumentError("Cannot add operators with different L_max cutoffs: $(op1.L_max) ≠ $(op2.L_max)")
    end
    if op1.universe != op2.universe
        return ArgumentError("Cannot add operators from different universes: $(op1.universe) ≠ $(op2.universe)")
    end

    return EDOperator(op1.lattice, op1.matrix + op2.matrix, op1.L_max, op1.universe)
end

function Base.:*(op1::EDOperator{N,F}, op2::EDOperator{N,F}) where {N,F}
    if op1.lattice != op2.lattice
        return ArgumentError("Cannot multiply operators on different lattices")
    end
    if op1.L_max != op2.L_max
        return ArgumentError("Cannot multiply operators with different L_max cutoffs: $(op1.L_max) ≠ $(op2.L_max)")
    end
    if op1.universe != op2.universe
        return ArgumentError("Cannot multiply operators from different universes: $(op1.universe) ≠ $(op2.universe)")
    end

    return EDOperator(op1.lattice, op1.matrix * op2.matrix, op1.L_max, op1.universe)
end

@inline function Base.:*(scalar::Number, op::EDOperator{N,F}) where {N,F}
    return EDOperator(op.lattice, scalar * op.matrix, op.L_max, op.universe)
end

@inline function Base.:*(op::SchwingerOperator{N,F}, scalar::Number) where {N,F}
    return scalar * op
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

function Base.:+(op1::MPOOperator{N,F}, op2::MPOOperator{N,F}) where {N,F}
    if op1.lattice != op2.lattice
        return ArgumentError("Cannot add operators on different lattices")
    end
    if op1.L_max != op2.L_max
        return ArgumentError("Cannot add operators with different L_max cutoffs: $(op1.L_max) ≠ $(op2.L_max)")
    end
    if op1.universe != op2.universe
        return ArgumentError("Cannot add operators from different universes: $(op1.universe) ≠ $(op2.universe)")
    end

    return MPOOperator{N,F}(op1.lattice, op1.mpo + op2.mpo, op1.L_max, op1.universe)
end

function Base.:*(op1::MPOOperator{N,F}, op2::MPOOperator{N,F}) where {N,F}
    if op1.lattice != op2.lattice
        return ArgumentError("Cannot multiply operators on different lattices")
    end
    if op1.L_max != op2.L_max
        return ArgumentError("Cannot multiply operators with different L_max cutoffs: $(op1.L_max) ≠ $(op2.L_max)")
    end
    if op1.universe != op2.universe
        return ArgumentError("Cannot multiply operators from different universes: $(op1.universe) ≠ $(op2.universe)")
    end

    return MPOOperator(op1.lattice, ITensors.apply(op1.mpo, op2.mpo), op1.L_max, op1.universe)
end

@inline function Base.:*(scalar::Number, op::MPOOperator{N,F}) where {N,F}
    return MPOOperator(op.lattice, scalar * op.mpo, op.L_max, op.universe)
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
    
    mpo = MPO(holonomy, sites(lattice; L_max=L_max))
    return MPOOperator(lattice, mpo, L_max, universe)
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