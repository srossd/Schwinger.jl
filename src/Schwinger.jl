module Schwinger

export SchwingerLattice
export basis, EDHamiltonian, MPOHamiltonian
export expectation, EDWilsonLoop, MPOWilsonLoop
export loweststates, groundstate, energygap
export energy, occupations, charges, electricfields, entanglements, scalarvev

using ITensors, ITensorMPS
using Combinatorics
using SparseArrays
using Arpack

using Memoize
using Parameters

include("./lattice.jl")
include("./operators.jl")
include("./hamiltonian.jl")
include("./states.jl")
include("./utility.jl")

"""
A Julia package for the Hamiltonian lattice Schwinger model.
"""
Schwinger

end
