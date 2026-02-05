module Schwinger

# Lattice
export Lattice, LatticeSize

# Backend types
export Backend, EDBackend, ITensorsBackend, MPSKitBackend
export ED, ITensorsMPS, MPSKitMPS
export resolve_backend, set_default_backend, get_default_backend

# State types
export SchwingerState, BasisState, EDState, ITensorState, MPSKitState

# Operator types
export SchwingerOperator, EDOperator, ITensorOperator, MPSKitOperator

# Unified API (backend-agnostic)
export Hamiltonian, GaugeKinetic, Mass, Hopping, HoppingMass
export WilsonLoop, WilsonLine, AverageElectricField

# Backend-specific operators (for backward compatibility)
export EDHamiltonian, ITensorHamiltonian, MPSKitHamiltonian
export EDGaugeKinetic, EDMass, EDHopping, EDHoppingMass
export ITensorGaugeKinetic, ITensorMass, ITensorHopping, ITensorHoppingMass
export MPSKitGaugeKinetic, MPSKitMass, MPSKitHopping, MPSKitHoppingMass
export EDWilsonLoop, ITensorWilsonLoop, MPSKitWilsonLoop
export EDWilsonLine, ITensorWilsonLine, MPSKitWilsonLine
export EDAverageElectricField, ITensorAverageElectricField, MPSKitAverageElectricField
# State operations
export lattice
export loweststates, groundstate, energygap
export expectation, act
export evolve

# Observables
export energy, energy_density, occupation, occupations
export charge, charges
export L₀, electricfield, electricfields
export entanglement, entanglements
export scalar, scalardensity, scalardensities
export pseudoscalar, pseudoscalardensity, pseudoscalardensities

# ITensors and ITensorMPS
using ITensors, ITensorMPS
using Observers

# MPSKit and TensorKit
using MPSKit, MPSKitLEMPO
using TensorKit

# Numerical libraries
using Arpack
using KrylovKit

# Standard libraries
using Combinatorics
using Statistics
using LinearAlgebra
using SparseArrays
using PeriodicArrays
using DataFrames

# Utilities
using Memoize

# Include files in dependency order
include("./lattice.jl")
include("./backends.jl")
include("./operators.jl")
include("./utility.jl")
include("./states.jl")
include("./hamiltonian.jl")
include("./wilson.jl")
include("./average_electric_field.jl")
include("./timeevolution.jl")
include("./api.jl")
# include("./hdf5.jl")

# =============================================================================
# Change MPSKit.excitations definition to allow for LEMPOS
# =============================================================================

function __init__()
    @eval function MPSKit.excitations(
            H::FiniteLEMPOHamiltonian, alg::MPSKit.FiniteExcited,
            states::Tuple{T, Vararg{T}};
            init = MPSKit.FiniteMPS(
                [ copy(first(states).AC[i]) for i in 1:length(first(states)) ]
            ), num = 1
        ) where {T <: MPSKit.FiniteMPS}
        num == 0 && return (scalartype(T)[], T[])

        super_op = MPSKit.LinearCombination(
            tuple(H, MPSKit.ProjectionOperator.(states)...),
            tuple(1.0, broadcast(x -> alg.weight, states)...)
        )
        envs = MPSKit.environments(init, super_op)
        ne, _ = MPSKit.find_groundstate(init, super_op, alg.gsalg, envs)

        nstates = (states..., ne)
        ens, excis = MPSKit.excitations(H, alg, nstates; init = init, num = num - 1)

        push!(ens, MPSKit.expectation_value(ne, H))
        push!(excis, ne)

        return ens, excis
    end
    @eval function MPSKit.excitations(H, alg::MPSKit.FiniteExcited, ψ::MPSKit.FiniteMPS; kwargs...)
        return MPSKit.excitations(H, alg, (ψ,); kwargs...)
    end
end

# =============================================================================
# END HACK
# =============================================================================

"""
A Julia package for the Hamiltonian lattice Schwinger model.

Supports three computational backends:
- `:ED` - Exact diagonalization using sparse matrices
- `:ITensors` - ITensors.jl MPO/MPS representations
- `:MPSKit` - MPSKitLEMPO.jl (LE)MPO/MPS representations
"""
Schwinger

end
