module Schwinger

export SchwingerLattice
export SchwingerState, SchwingerBasisState, SchwingerEDState, SchwingerMPS
export SchwingerOperator, EDOperator, MPOOperator
export basis, EDHamiltonian, MPOHamiltonian
export expectation, EDWilsonLoop, MPOWilsonLoop, EDWilsonLine, MPOWilsonLine
export EDGaugeKinetic, EDMass, EDHopping, EDHoppingMass
export MPOGaugeKinetic, MPOMass, MPOHopping, MPOHoppingMass
export EDAverageElectricField, MPOAverageElectricField
export loweststates, groundstate, energygap
export energy, occupation, occupations, charge, charges, L₀, electricfield, electricfields, entanglement, entanglements
export scalar, scalardensity, scalardensities
export pseudoscalar, pseudoscalardensity, pseudoscalardensities
export expectation, act
export evolve

using ITensors, ITensorMPS
using Observers

using Arpack
using KrylovKit

using Combinatorics
using Statistics
using LinearAlgebra
using SparseArrays
using DataFrames
using HDF5

using Memoize
using Parameters

include("./lattice.jl")
include("./operators.jl")
include("./states.jl")
include("./hamiltonian.jl")
include("./wilson.jl")
include("./average_electric_field.jl")
include("./timeevolution.jl")
include("./utility.jl")
include("./hdf5.jl")

"""
A Julia package for the Hamiltonian lattice Schwinger model.
"""
Schwinger

end
