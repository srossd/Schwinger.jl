var documenterSearchIndex = {"docs":
[{"location":"index.html#Schwinger.jl","page":"Index","title":"Schwinger.jl","text":"","category":"section"},{"location":"index.html","page":"Index","title":"Index","text":"Schwinger","category":"page"},{"location":"index.html#Schwinger","page":"Index","title":"Schwinger","text":"A Julia package for the Hamiltonian lattice Schwinger model.\n\n\n\n\n\n","category":"module"},{"location":"index.html#Module-Index","page":"Index","title":"Module Index","text":"","category":"section"},{"location":"index.html","page":"Index","title":"Index","text":"Modules = [Schwinger]\nOrder   = [:constant, :type, :function, :macro]","category":"page"},{"location":"index.html#Detailed-API","page":"Index","title":"Detailed API","text":"","category":"section"},{"location":"index.html","page":"Index","title":"Index","text":"Modules = [Schwinger]\nOrder   = [:constant, :type, :function, :macro]","category":"page"},{"location":"index.html#Schwinger.EDHamiltonian-Union{Tuple{SchwingerLattice{N, F}}, Tuple{F}, Tuple{N}} where {N, F}","page":"Index","title":"Schwinger.EDHamiltonian","text":"EDHamiltonian(lattice) Computes the Hamiltonian for the Schwinger model.\n\nArguments\n\nlattice::SchwingerLattice: Schwinger model lattice.\n\n\n\n\n\n","category":"method"},{"location":"index.html#Schwinger.MPOHamiltonian-Union{Tuple{SchwingerLattice{N, F}}, Tuple{F}, Tuple{N}} where {N, F}","page":"Index","title":"Schwinger.MPOHamiltonian","text":"MPOhamiltonian(lattice)\n\nComputes the MPO Hamiltonian for the Schwinger model.\n\nArguments\n\nlattice::SchwingerLattice: Schwinger model lattice.\n\n\n\n\n\n","category":"method"},{"location":"index.html#Schwinger.SchwingerBasisState","page":"Index","title":"Schwinger.SchwingerBasisState","text":"SchwingerBasisState{N,F}(occupations, L0)\n\nA Schwinger model basis state.\n\n\n\n\n\n","category":"type"},{"location":"index.html#Schwinger.SchwingerEDState","page":"Index","title":"Schwinger.SchwingerEDState","text":"SchwingerEDState{N,F}(hamiltonian, coeffs)\n\nA Schwinger model state represented as a linear combination of basis states.\n\n\n\n\n\n","category":"type"},{"location":"index.html#Schwinger.SchwingerLattice","page":"Index","title":"Schwinger.SchwingerLattice","text":"SchwingerLattice(N, F; periodic=false, L_max=0)\n\nConstructs a SchwingerLattice for the Schwinger model.\n\nArguments\n\nN::Int: Number of sites.\nF::Int: Number of flavors.\nperiodic::Bool=false: Whether the lattice is periodic.\nL_max::Int=0: Maximum absolute value for L_0 when periodic.\n\nReturns\n\nA SchwingerLattice object.\n\n\n\n\n\n","category":"type"},{"location":"index.html#Schwinger.charges-Union{Tuple{Schwinger.SchwingerState{N, F}}, Tuple{F}, Tuple{N}} where {N, F}","page":"Index","title":"Schwinger.charges","text":"charges(state)\n\nReturn a list of the expectations of Q operators on each site and for each known eigenstate.\n\nArguments\n\nstate::SchwingerState: Schwinger model state.\n\n\n\n\n\n","category":"method"},{"location":"index.html#Schwinger.electricfields-Union{Tuple{Schwinger.SchwingerState{N, F}}, Tuple{F}, Tuple{N}} where {N, F}","page":"Index","title":"Schwinger.electricfields","text":"electricfields(state)\n\nReturn a list of the expectations of (L + θ/2π) operators on each link.\n\nArguments\n\nstate::SchwingerMPS: Schwinger model state.\n\n\n\n\n\n","category":"method"},{"location":"index.html#Schwinger.energy-Union{Tuple{Schwinger.SchwingerBasisState{N, F}}, Tuple{F}, Tuple{N}} where {N, F}","page":"Index","title":"Schwinger.energy","text":"energy(state)\n\nReturn the expectation value of the Hamiltonian.\n\nArguments\n\nstate::SchwingerBasisState: Schwinger model basis state.\n\n\n\n\n\n","category":"method"},{"location":"index.html#Schwinger.energy-Union{Tuple{Schwinger.SchwingerEDState{N, F}}, Tuple{F}, Tuple{N}} where {N, F}","page":"Index","title":"Schwinger.energy","text":"energy(state)\n\nReturn the expectation value of the Hamiltonian.\n\nArguments\n\nstate::SchwingerEDState: Schwinger model state.\n\n\n\n\n\n","category":"method"},{"location":"index.html#Schwinger.energy-Union{Tuple{Schwinger.SchwingerMPS{N, F}}, Tuple{F}, Tuple{N}} where {N, F}","page":"Index","title":"Schwinger.energy","text":"energy(state)\n\nReturn the expectation value of the Hamiltonian.\n\nArguments\n\nstate::SchwingerMPS: Schwinger model state.\n\n\n\n\n\n","category":"method"},{"location":"index.html#Schwinger.entanglements-Union{Tuple{Schwinger.SchwingerMPS{N, F}}, Tuple{F}, Tuple{N}} where {N, F}","page":"Index","title":"Schwinger.entanglements","text":"entanglements(state)\n\nReturn a list of the entanglement entropies for each bisection of the lattice.\n\nArguments\n\nstate::SchwingerMPS: Schwinger model state.\n\n\n\n\n\n","category":"method"},{"location":"index.html#Schwinger.occupations-Union{Tuple{Schwinger.SchwingerBasisState{N, F}}, Tuple{F}, Tuple{N}} where {N, F}","page":"Index","title":"Schwinger.occupations","text":"occupations(state)\n\nReturn an NxF matrix of the expectations of χ†χ operators on each site.\n\nArguments\n\nstate::SchwingerBasisState: Schwinger model basis state.\n\n\n\n\n\n","category":"method"},{"location":"index.html#Schwinger.occupations-Union{Tuple{Schwinger.SchwingerEDState{N, F}}, Tuple{F}, Tuple{N}} where {N, F}","page":"Index","title":"Schwinger.occupations","text":"occupations(state)\n\nReturn an NxF matrix of the expectations of χ†χ operators on each site.\n\nArguments\n\nstate::SchwingerBasisState: Schwinger model basis state.\n\n\n\n\n\n","category":"method"},{"location":"index.html#Schwinger.occupations-Union{Tuple{Schwinger.SchwingerMPS{N, F}}, Tuple{F}, Tuple{N}} where {N, F}","page":"Index","title":"Schwinger.occupations","text":"occupations(state)\n\nReturn an NxF matrix of the expectations of χ†χ operators on each site.\n\nArguments\n\nstate::SchwingerMPS{N,F}: Schwinger model state.\n\n\n\n\n\n","category":"method"}]
}
