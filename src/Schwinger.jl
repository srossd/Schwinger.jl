module Schwinger

export SchwingerLattice
export basis, EDHamiltonian, MPOHamiltonian
export expectation, EDWilsonLoop, MPOWilsonLoop, EDWilsonLine, MPOWilsonLine
export EDAverageElectricField, MPOAverageElectricField
export loweststates, groundstate, energygap
export energy, occupations, charges, electricfields, entanglements, scalarvev, pseudoscalarvev
export act, evolve

using ITensors, ITensorMPS
using Combinatorics
using Statistics
using SparseArrays
using Arpack
using KrylovKit

using Memoize
using Parameters

include("./lattice.jl")
include("./operators.jl")
include("./states.jl")
include("./hamiltonian.jl")
include("./wilson.jl")
include("./average_electric_field.jl")
include("./utility.jl")

"""
A Julia package for the Hamiltonian lattice Schwinger model.
"""
Schwinger

end
