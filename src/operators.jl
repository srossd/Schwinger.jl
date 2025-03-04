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