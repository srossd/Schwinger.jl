abstract type SchwingerOperator end

"""
    validate_operator_compatibility(op1, op2, operation)

Validate that two operators are compatible for the given operation.
Throws ArgumentError if operators have different lattices, L_max, or universes.
"""
function validate_operator_compatibility(op1::SchwingerOperator, op2::SchwingerOperator, operation::String)
    if op1.lattice != op2.lattice
        throw(ArgumentError("Cannot $operation operators on different lattices"))
    end
    if op1.L_max != op2.L_max
        throw(ArgumentError("Cannot $operation operators with different L_max cutoffs: $(op1.L_max) ≠ $(op2.L_max)"))
    end
    if op1.universe != op2.universe
        throw(ArgumentError("Cannot $operation operators from different universes: $(op1.universe) ≠ $(op2.universe)"))
    end
end

# Generic scalar multiplication (right multiply)
@inline function Base.:*(op::SchwingerOperator, scalar::Number) 
    return scalar * op
end

# =============================================================================
# ED Backend
# =============================================================================

struct EDOperator <: SchwingerOperator
    lattice::Lattice
    matrix::SparseMatrixCSC{ComplexF64,Int64}
    L_max::Int
    universe::Int

    function EDOperator(lattice::Lattice, matrix::SparseMatrixCSC{ComplexF64,Int64}, L_max::Int, universe::Int)
        N, F = lattice.N, lattice.F
        return new(lattice, matrix, L_max, universe)
    end
end

function Base.:+(op1::EDOperator, op2::EDOperator) 
    validate_operator_compatibility(op1, op2, "add")
    return EDOperator(op1.lattice, op1.matrix + op2.matrix, op1.L_max, op1.universe)
end

function Base.:*(op1::EDOperator, op2::EDOperator) 
    validate_operator_compatibility(op1, op2, "multiply")
    return EDOperator(op1.lattice, op1.matrix * op2.matrix, op1.L_max, op1.universe)
end

@inline function Base.:*(scalar::Number, op::EDOperator) 
    return EDOperator(op.lattice, scalar * op.matrix, op.L_max, op.universe)
end

# =============================================================================
# ITensors Backend
# =============================================================================

struct ITensorOperator <: SchwingerOperator
    lattice::Lattice
    mpo::ITensorMPS.MPO
    L_max::Int
    universe::Int

    function ITensorOperator(lattice::Lattice, mpo::ITensorMPS.MPO, L_max::Int, universe::Int)
        N, F = lattice.N, lattice.F
        return new(lattice, mpo, L_max, universe)
    end
end

function Base.:+(op1::ITensorOperator, op2::ITensorOperator) 
    validate_operator_compatibility(op1, op2, "add")
    return ITensorOperator(op1.lattice, op1.mpo + op2.mpo, op1.L_max, op1.universe)
end

function Base.:*(op1::ITensorOperator, op2::ITensorOperator) 
    validate_operator_compatibility(op1, op2, "multiply")
    return ITensorOperator(op1.lattice, ITensors.apply(op1.mpo, op2.mpo), op1.L_max, op1.universe)
end

@inline function Base.:*(scalar::Number, op::ITensorOperator) 
    return ITensorOperator(op.lattice, scalar * op.mpo, op.L_max, op.universe)
end

# =============================================================================
# MPSKit Backend
# =============================================================================

struct MPSKitOperator <: SchwingerOperator
    lattice::Lattice
    lempo::MPSKit.AbstractMPO
    universe::Int

    function MPSKitOperator(lattice::Lattice, lempo::MPSKit.AbstractMPO, universe::Int = 0)
        return new(lattice, lempo, universe)
    end
end

function Base.:+(op1::MPSKitOperator, op2::MPSKitOperator)
    validate_operator_compatibility(op1, op2, "add")
    return MPSKitOperator(op1.lattice, op1.lempo + op2.lempo, op1.universe)
end

function Base.:*(op1::MPSKitOperator, op2::MPSKitOperator)
    validate_operator_compatibility(op1, op2, "multiply")
    return MPSKitOperator(op1.lattice, op1.lempo * op2.lempo, op1.universe)
end

@inline function Base.:*(scalar::Number, op::MPSKitOperator)
    return MPSKitOperator(op.lattice, scalar * op.lempo, op.universe)
end
