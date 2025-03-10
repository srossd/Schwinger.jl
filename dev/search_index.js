var documenterSearchIndex = {"docs":
[{"location":"man/timeevolution.html#Time-evolution","page":"Time evolution","title":"Time evolution","text":"","category":"section"},{"location":"man/timeevolution.html","page":"Time evolution","title":"Time evolution","text":"Schwinger.jl supports time-evolving states. In the exact diagonalization framework, this is accomplished using Krylov methods from KrylovKit.jl. With matrix product operators, this is accomplished using the time-dependent variational principle algorithm (TDVP), as implemented in ITensorMPS.jl.","category":"page"},{"location":"man/timeevolution.html","page":"Time evolution","title":"Time evolution","text":"The evolve function evolves a state forwards in time, and monitors any given observables. It returns a final state along with a DataFrame of the observables. For example, here is a simulation of flux unwinding.","category":"page"},{"location":"man/timeevolution.html","page":"Time evolution","title":"Time evolution","text":"using Schwinger, Plots\n\nlat = SchwingerLattice{10,1}(L = 2, periodic = true)\ngs = groundstate(EDHamiltonian(lat))\n\n_, df = evolve(EDWilsonLoop(lat) * gs, 15; \n    nsteps = 30, \n    observable = (ψ, t) -> sum(electricfields(ψ))/10\n)\n\nscatter(df.time, df.observable, xlabel = \"gt\", ylabel = \"Average electric field\", label = \"Schwinger.jl\")\nplot!(0:.1:15, [cos(t/√(π)) for t in 0:.1:15], label = \"Exact\")","category":"page"},{"location":"man/timeevolution.html","page":"Time evolution","title":"Time evolution","text":"evolve","category":"page"},{"location":"man/timeevolution.html#Schwinger.evolve","page":"Time evolution","title":"Schwinger.evolve","text":"evolve(state::SchwingerEDState{N,F}, t::Real, observable::Function)\n\nEvolve the state by a time t, monitoring an observable.\n\nArguments\n\nstate::SchwingerEDState{N,F}: The state to evolve.\nt::Real: The time to evolve by.\nobservable::Function: The observable (ψ, t) -> obs to monitor.\n\n\n\n\n\nevolve(state::SchwingerMPS{N,F}, t::Real)\n\nEvolve the state by a time t.\n\nArguments\n\nstate::SchwingerMPS{N,F}: The state to evolve.\nt::Real: The time to evolve by.\n\n\n\n\n\n","category":"function"},{"location":"man/states.html#States","page":"States","title":"States","text":"","category":"section"},{"location":"man/states.html","page":"States","title":"States","text":"Lattice states in Schwinger.jl are represented by the abstract type SchwingerState, with three descendants: ","category":"page"},{"location":"man/states.html","page":"States","title":"States","text":"SchwingerBasisState: a state specified by the eigenvalues of occupation operators chi^dagger_nalphachi_nalpha and L_0\nSchwingerEDState: a linear combination of SchwingerBasisStates\nSchwingerMPS: a matrix product state, stored as an MPS object using ITensorMPS.jl","category":"page"},{"location":"man/states.html","page":"States","title":"States","text":"Given a state, we can find the expectation values of the occupation operators and electric field operators:","category":"page"},{"location":"man/states.html","page":"States","title":"States","text":"using Schwinger\nlat = SchwingerLattice{6,1}(a = 10) # towards the lattice strong coupling limit ga -> infty\ngs = groundstate(EDHamiltonian(lat))\n\noccupations(gs), electricfields(gs)","category":"page"},{"location":"man/states.html","page":"States","title":"States","text":"We can also evaluate the entanglement entropies of each bisection of the lattice:","category":"page"},{"location":"man/states.html","page":"States","title":"States","text":"using Schwinger\nlat = SchwingerLattice{20,1}()\ngs = groundstate(MPOHamiltonian(lat))\nentanglements(gs)","category":"page"},{"location":"man/states.html","page":"States","title":"States","text":"Several other useful functions are detailed below.","category":"page"},{"location":"man/states.html","page":"States","title":"States","text":"SchwingerState\nSchwingerBasisState\nSchwingerEDState\nSchwingerMPS\noccupations\ncharges\nelectricfields\nentanglements\nenergy\nL₀\nscalarvev\npseudoscalarvev","category":"page"},{"location":"man/states.html#Schwinger.SchwingerState","page":"States","title":"Schwinger.SchwingerState","text":"SchwingerState{N,F}\n\nAbstract type for Schwinger model states.\n\n\n\n\n\n","category":"type"},{"location":"man/states.html#Schwinger.SchwingerBasisState","page":"States","title":"Schwinger.SchwingerBasisState","text":"SchwingerBasisState{N,F}(occupations, L0)\n\nA Schwinger model basis state.\n\n\n\n\n\n","category":"type"},{"location":"man/states.html#Schwinger.SchwingerEDState","page":"States","title":"Schwinger.SchwingerEDState","text":"SchwingerEDState{N,F}(hamiltonian, coeffs)\n\nA Schwinger model state represented as a linear combination of basis states.\n\n\n\n\n\n","category":"type"},{"location":"man/states.html#Schwinger.SchwingerMPS","page":"States","title":"Schwinger.SchwingerMPS","text":"SchwingerMPSState{N,F}(hamiltonian, psi)\n\nA Schwinger model MPS.\n\n\n\n\n\n","category":"type"},{"location":"man/states.html#Schwinger.occupations","page":"States","title":"Schwinger.occupations","text":"occupations(state)\n\nReturn an NxF matrix of the expectations of χ†χ operators on each site.\n\nArguments\n\nstate::SchwingerMPS{N,F}: Schwinger model state.\n\n\n\n\n\noccupations(state)\n\nReturn an NxF matrix of the expectations of χ†χ operators on each site.\n\nArguments\n\nstate::SchwingerBasisState: Schwinger model basis state.\n\n\n\n\n\noccupations(state)\n\nReturn an NxF matrix of the expectations of χ†χ operators on each site.\n\nArguments\n\nstate::SchwingerBasisState: Schwinger model basis state.\n\n\n\n\n\n","category":"function"},{"location":"man/states.html#Schwinger.charges","page":"States","title":"Schwinger.charges","text":"charges(state)\n\nReturn a list of the expectations of Q operators on each site and for each known eigenstate.\n\nArguments\n\nstate::SchwingerState: Schwinger model state.\n\n\n\n\n\n","category":"function"},{"location":"man/states.html#Schwinger.electricfields","page":"States","title":"Schwinger.electricfields","text":"electricfields(state)\n\nReturn a list of the expectations of (L + θ/2π) operators on each link.\n\nArguments\n\nstate::SchwingerMPS: Schwinger model state.\n\n\n\n\n\n","category":"function"},{"location":"man/states.html#Schwinger.entanglements","page":"States","title":"Schwinger.entanglements","text":"entanglements(state)\n\nReturn a list of the von Neumann entanglement entropies for each bisection of the lattice.\n\nArguments\n\nstate::SchwingerState: Schwinger model state.\n\n\n\n\n\n","category":"function"},{"location":"man/states.html#Schwinger.energy","page":"States","title":"Schwinger.energy","text":"energy(state)\n\nReturn the expectation value of the Hamiltonian.\n\nArguments\n\nstate::SchwingerEDState: Schwinger model state.\n\n\n\n\n\nenergy(state)\n\nReturn the expectation value of the Hamiltonian.\n\nArguments\n\nstate::SchwingerBasisState: Schwinger model basis state.\n\n\n\n\n\n","category":"function"},{"location":"man/states.html#Schwinger.L₀","page":"States","title":"Schwinger.L₀","text":"L₀(state)\n\nReturn the expectation value of L₀.\n\nArguments\n\nstate::SchwingerState: Schwinger model state.\n\n\n\n\n\n","category":"function"},{"location":"man/states.html#Schwinger.scalarvev","page":"States","title":"Schwinger.scalarvev","text":"scalarvev(state)\n\nReturn the VEV of the scalar condensate L⁻¹ ∑ (-1)ⁿ χ†ₙχₙ\n\nArguments\n\nstate::SchwingerState: Schwinger model state.\n\n\n\n\n\n","category":"function"},{"location":"man/states.html#Schwinger.pseudoscalarvev","page":"States","title":"Schwinger.pseudoscalarvev","text":"pseudoscalarvev(state)\n\nReturn the VEV of the pseudoscalar condensate L⁻¹ ∑ (-1)ⁿ χ†ₙχₙ\n\nArguments\n\nstate::SchwingerState: Schwinger model state.\n\n\n\n\n\n","category":"function"},{"location":"man/operators.html#Operators","page":"Operators","title":"Operators","text":"","category":"section"},{"location":"man/operators.html","page":"Operators","title":"Operators","text":"Operators in Schwinger.jl are represented by the abstract type SchwingerOperator, with two descendants:","category":"page"},{"location":"man/operators.html","page":"Operators","title":"Operators","text":"EDOperator: an operator represented as a matrix for its action on a basis of SchwingerBasisStates\nMPOOperator: a matrix product operator, stored as an MPO object using ITensorMPS.jl","category":"page"},{"location":"man/operators.html","page":"Operators","title":"Operators","text":"We can evaluate the expectation of an operator in a state:","category":"page"},{"location":"man/operators.html","page":"Operators","title":"Operators","text":"using Schwinger\nlat = SchwingerLattice{10,1}()\ngs = groundstate(EDHamiltonian(lat))\n\nsum(electricfields(gs))/10, expectation(EDAverageElectricField(lat), gs)","category":"page"},{"location":"man/operators.html","page":"Operators","title":"Operators","text":"This can also be carried out manually by acting on the state with the operator:","category":"page"},{"location":"man/operators.html","page":"Operators","title":"Operators","text":"using LinearAlgebra\ndot(gs, EDAverageElectricField(lat) * gs)","category":"page"},{"location":"man/operators.html","page":"Operators","title":"Operators","text":"expectation\nact","category":"page"},{"location":"man/operators.html#Schwinger.expectation","page":"Operators","title":"Schwinger.expectation","text":"expectation(op, state)\n\nReturn the expectation value of the operator op in state.\n\nArguments\n\nop::EDOperator{N,F}`: operator.\nstate::SchwingerEDState: state.\n\n\n\n\n\nexpectation(op, state)\n\nReturn the expectation value of the operator op in state.\n\nArguments\n\nop::MPOOperator{N,F}`: operator.\nstate::SchwingerMPS: state.\n\n\n\n\n\n","category":"function"},{"location":"man/operators.html#Schwinger.act","page":"Operators","title":"Schwinger.act","text":"act(op, state)\n\nApply the operator op to the state state.\n\nArguments\n\nop::MPOOperator{N,F}: operator.\nstate::SchwingerMPS{N,F}: state.\n\n\n\n\n\nact(op, state)\n\nApply the operator op to the state state.\n\nArguments\n\nop::EDOperator{N,F}: operator.\nstate::SchwingerEDState{N,F}: state.\n\n\n\n\n\n","category":"function"},{"location":"man/operators.html#Exact-diagonalization","page":"Operators","title":"Exact diagonalization","text":"","category":"section"},{"location":"man/operators.html","page":"Operators","title":"Operators","text":"EDWilsonLoop\nEDWilsonLine\nEDAverageElectricField","category":"page"},{"location":"man/operators.html#Schwinger.EDWilsonLoop","page":"Operators","title":"Schwinger.EDWilsonLoop","text":"wilsonloop(hamiltonian, conjugate = false)\n\nReturns the spatial Wilson loop operator for lattice.\n\nArguments\n\nlattice::SchwingerLattice: lattice.\nconjugate::Bool: Conjugation of the Wilson loop.\n\n\n\n\n\n","category":"function"},{"location":"man/operators.html#Schwinger.EDWilsonLine","page":"Operators","title":"Schwinger.EDWilsonLine","text":"EDWilsonLine(lattice, conjugate = false, flavor = 1, start = 1, finish = N)\n\nReturns the spatial Wilson line operator for lattice.\n\nArguments\n\nlattice::SchwingerLattice: lattice.\nconjugate::Bool: Conjugation of the Wilson line.\nflavor::Int: Flavor of the Wilson line.\nstart::Int: Starting site of the Wilson line.\nfinish::Int: Finishing site of the Wilson line.\n\n\n\n\n\n","category":"function"},{"location":"man/operators.html#Schwinger.EDAverageElectricField","page":"Operators","title":"Schwinger.EDAverageElectricField","text":"EDAverageElectricField(lattice; power = 1, L_max = nothing, universe = 0)\n\nConstruct an EDOperator that computes the average electric field (raised to some power).\n\nArguments\n\nlattice::SchwingerLattice{N,F}: The lattice to compute the average electric field on.\npower::Int=1: The power to raise the electric field to.\nL_max::Union{Nothing,Int}=nothing: The maximum absolute value of L₀.\nuniverse::Int=0: The universe to compute the average electric field in.\nsitelist::Union{Nothing,Vector{Int}}=nothing: List of sites to average over.\n\n\n\n\n\n","category":"function"},{"location":"man/operators.html#Matrix-product-operators","page":"Operators","title":"Matrix product operators","text":"","category":"section"},{"location":"man/operators.html","page":"Operators","title":"Operators","text":"MPOWilsonLoop\nMPOWilsonLine\nMPOAverageElectricField","category":"page"},{"location":"man/operators.html#Schwinger.MPOWilsonLoop","page":"Operators","title":"Schwinger.MPOWilsonLoop","text":"wilsonloop(hamiltonian, conjugate = false)\n\nReturns the spatial Wilson loop operator for lattice.\n\nArguments\n\nlattice::SchwingerLattice: lattice.\nconjugate::Bool: Conjugation of the Wilson loop.\n\n\n\n\n\n","category":"function"},{"location":"man/operators.html#Schwinger.MPOWilsonLine","page":"Operators","title":"Schwinger.MPOWilsonLine","text":"MPOWilsonLine(lattice, conjugate = false, flavor = 1, start = 1, finish = N)\n\nReturns the spatial Wilson line operator for lattice.\n\nArguments\n\nlattice::SchwingerLattice: lattice.\nconjugate::Bool: Conjugation of the Wilson line.\nflavor::Int: Flavor of the Wilson line.\nstart::Int: Starting site of the Wilson line.\nfinish::Int: Finishing site of the Wilson line.\n\n\n\n\n\n","category":"function"},{"location":"man/operators.html#Schwinger.MPOAverageElectricField","page":"Operators","title":"Schwinger.MPOAverageElectricField","text":"MPOAverageElectricField(lattice; power = 1, L_max = nothing, universe = 0)\n\nConstruct an MPOOperator that computes the average electric field (raised to some power).\n\nArguments\n\nlattice::SchwingerLattice{N,F}: The lattice to compute the average electric field on.\npower::Int=1: The power to raise the electric field to.\nL_max::Union{Nothing,Int}=nothing: The maximum absolute value of L₀.\nuniverse::Int=0: The universe to compute the average electric field in.\nsitelist::Union{Nothing,Vector{Int}}=nothing: List of sites to average over.\n\n\n\n\n\n","category":"function"},{"location":"index.html#Schwinger.jl","page":"Index","title":"Schwinger.jl","text":"","category":"section"},{"location":"index.html","page":"Index","title":"Index","text":"A package for working with the Hamiltonian lattice Schwinger model.","category":"page"},{"location":"index.html#Table-of-contents","page":"Index","title":"Table of contents","text":"","category":"section"},{"location":"index.html","page":"Index","title":"Index","text":"Pages = [\"home.md\", \"man/lattices.md\", \"man/hamiltonian.md\", \"man/states.md\", \"man/operators.md\", \"man/timeevolution.md\"]\nDepth = 4","category":"page"},{"location":"index.html#Installation","page":"Index","title":"Installation","text":"","category":"section"},{"location":"index.html","page":"Index","title":"Index","text":"To install Schwinger.jl, use the following command in the Julia REPL:","category":"page"},{"location":"index.html","page":"Index","title":"Index","text":"using Pkg\nPkg.add(\"https://github.com/srossd/Schwinger.jl\")","category":"page"},{"location":"index.html#Key-Features","page":"Index","title":"Key Features","text":"","category":"section"},{"location":"index.html","page":"Index","title":"Index","text":"Schwinger.jl provides functions for computing the following properties of the Hamiltonian lattice Schwinger model:","category":"page"},{"location":"index.html","page":"Index","title":"Index","text":"Energy levels: Ground state, energy gap, or higher excited states.\nCorrelators: Expectation values of various operators in the ground state (or other states).\nTime evolution: Action of time evolution on a given initial state.","category":"page"},{"location":"index.html","page":"Index","title":"Index","text":"All of these features are implemented using both exact diagonalization (ED) and matrix product operators/states (MPO).","category":"page"},{"location":"index.html#Usage-Example","page":"Index","title":"Usage Example","text":"","category":"section"},{"location":"index.html","page":"Index","title":"Index","text":"Here is a basic example of how to use Schwinger.jl to calculate the average electric field in the one-flavor Schwinger model at mg = 01 as a function of theta. Note that the mass shift (see here) is included by default.","category":"page"},{"location":"index.html","page":"Index","title":"Index","text":"using Schwinger\nusing Base.MathConstants\nusing Plots\n\nfunction avgE(θ2π, m, shift = true)\n    lat = shift ? \n        SchwingerLattice{10,1}(θ2π = θ2π, m = m, periodic = true) : \n        SchwingerLattice{10,1}(θ2π = θ2π, mlat = m, periodic = true)\n    gs = groundstate(EDHamiltonian(lat))\n    return real(expectation(EDAverageElectricField(lat), gs))\nend\n\nm = 0.05\nθ2πs = 0:0.025:1\navgEs_shift = map(x -> avgE(x, m), θ2πs)\navgEs_noshift = map(x -> avgE(x, m, false), θ2πs)\n\n# See eq (24) of https://arxiv.org/abs/2206.05308\nperturbative = [(exp(γ)/√(π))*m*sin(2π*θ2π) - \n    (8.9139*exp(2γ)/(4π))*(m^2)*sin(4π*θ2π) for θ2π in θ2πs]\n\nscatter(θ2πs, avgEs_noshift, \n    label=\"Without mass shift\", \n    xlabel=\"θ/2π\", \n    ylabel=\"Average Electric Field\", \n    color=\"lightblue\")\nscatter!(θ2πs, avgEs_shift, \n    label=\"With mass shift\", \n    xlabel=\"θ/2π\", \n    ylabel=\"Average Electric Field\", \n    legend=:topright, \n    color=:orange)\nplot!(θ2πs, perturbative, \n    label=\"Perturbative\", \n    color=:black)","category":"page"},{"location":"index.html","page":"Index","title":"Index","text":"Here is an example of calculating the expectation value of the square of the mean electric field in the two-flavor Schwinger model at theta = pi, for a four-site periodic lattice, giving a very rough look at the phase diagram.","category":"page"},{"location":"index.html","page":"Index","title":"Index","text":"using Schwinger\nusing Plots\n\nfunction avgE2(m1, m2)\n    lat = SchwingerLattice{4,2}(θ2π = 0.5, m = (m1, m2), periodic = true)\n    gs = groundstate(EDHamiltonian(lat))\n    return real(expectation(EDAverageElectricField(lat; power=2), gs))\nend\n\nms = -1:0.05:1\nheatmap(ms, ms, avgE2, xlabel = \"m₁/g\", ylabel = \"m₂/g\", c = :berlin)","category":"page"},{"location":"man/hamiltonian.html#Hamiltonian","page":"Hamiltonian","title":"Hamiltonian","text":"","category":"section"},{"location":"man/hamiltonian.html","page":"Hamiltonian","title":"Hamiltonian","text":"The Hamiltonian for the lattice Schwinger model is","category":"page"},{"location":"man/hamiltonian.html","page":"Hamiltonian","title":"Hamiltonian","text":"beginsplit\nH = frac(qg)^2 a2sum_n=1^N left(L_n + fractheta2piright)^2 - fraci2asum_n=1^Nsum_alpha=1^F left(chi^dagger_nalpha chi_n+1alpha - chi^dagger_n+1alpha chi_nalpharight) \n+underbraceleft(m - frac(qg)^2 F a8right)_m_textlatsum_n=1^N sum_alpha=1^F (-1)^n chi^dagger_nalpha chi_nalpha + m sum_n=1^Nsum_alpha=1^F (-1)^nleft(chi^dagger_n-1alphachi_nalpha + chi^dagger_n+1alphachi_nalpharight)\nendsplit","category":"page"},{"location":"man/hamiltonian.html","page":"Hamiltonian","title":"Hamiltonian","text":"This is supplemented by the Gauss law","category":"page"},{"location":"man/hamiltonian.html","page":"Hamiltonian","title":"Hamiltonian","text":"L_n = L_n-1 + Q_n qquad Q_n equiv qleft(sum_alpha=1^F chi^dagger_nalpha chi_nalpha - begincases F  ntext odd  0  ntext even endcasesright)","category":"page"},{"location":"man/hamiltonian.html","page":"Hamiltonian","title":"Hamiltonian","text":"In Schwinger.jl this be constructed using two strategies, exact diagonalization (ED) or matrix product operators (MPO).","category":"page"},{"location":"man/hamiltonian.html#Exact-diagonalization","page":"Hamiltonian","title":"Exact diagonalization","text":"","category":"section"},{"location":"man/hamiltonian.html","page":"Hamiltonian","title":"Hamiltonian","text":"When using exact diagonalization, Schwinger.jl constructs a basis of states that diagonalize the operators chi^dagger_nalpha chi_nalpha and L_0. It then builds a sparse matrix for the Hamiltonian acting on this basis.","category":"page"},{"location":"man/hamiltonian.html","page":"Hamiltonian","title":"Hamiltonian","text":"using Schwinger\nlat = SchwingerLattice{12,1}();\nham = EDHamiltonian(lat);\n\nham.matrix","category":"page"},{"location":"man/hamiltonian.html","page":"Hamiltonian","title":"Hamiltonian","text":"When q  1, the universe (i.e., the allowed values of L_n modulo q) can be specified; the default value is 0. When the lattice is periodic, the maximum absolute value of L_0 can be set using L_max; the default value is 3.","category":"page"},{"location":"man/hamiltonian.html","page":"Hamiltonian","title":"Hamiltonian","text":"Using Arpack.jl, Schwinger.jl can find the lowest eigenstates of a Hamiltonian.","category":"page"},{"location":"man/hamiltonian.html","page":"Hamiltonian","title":"Hamiltonian","text":"using Schwinger\nlat = SchwingerLattice{10,1}(periodic = true);\nham = EDHamiltonian(lat);\n\nmap(energy, loweststates(ham, 5))","category":"page"},{"location":"man/hamiltonian.html","page":"Hamiltonian","title":"Hamiltonian","text":"EDHamiltonian\nEDGaugeKinetic\nEDHopping\nEDMass\nEDHoppingMass","category":"page"},{"location":"man/hamiltonian.html#Schwinger.EDHamiltonian","page":"Hamiltonian","title":"Schwinger.EDHamiltonian","text":"EDHamiltonian(lattice) Computes the Hamiltonian for the Schwinger model.\n\nArguments\n\nlattice::SchwingerLattice: Schwinger model lattice.\n\n\n\n\n\n","category":"function"},{"location":"man/hamiltonian.html#Schwinger.EDGaugeKinetic","page":"Hamiltonian","title":"Schwinger.EDGaugeKinetic","text":"EDGaugeKinetic(lattice) Computes the gauge kinetic operator ∑(Lₙ+θ/2π)² for the Schwinger model.\n\nArguments\n\nlattice::SchwingerLattice: Schwinger model lattice.\n\n\n\n\n\n","category":"function"},{"location":"man/hamiltonian.html#Schwinger.EDHopping","page":"Hamiltonian","title":"Schwinger.EDHopping","text":"EDHopping(lattice) Computes the hopping term -i ∑(χ†ₙ χₙ₊₁ - χ†ₙ₊₁ χₙ) for the Schwinger model.\n\nArguments\n\nlattice::SchwingerLattice: Schwinger model lattice.\n\n\n\n\n\n","category":"function"},{"location":"man/hamiltonian.html#Schwinger.EDMass","page":"Hamiltonian","title":"Schwinger.EDMass","text":"EDMass(lattice) Computes the mass operator ∑ (-1)ⁿ χ†ₙχₙ for the Schwinger model.\n\nArguments\n\nlattice::SchwingerLattice: Schwinger model lattice.\n\n\n\n\n\n","category":"function"},{"location":"man/hamiltonian.html#Schwinger.EDHoppingMass","page":"Hamiltonian","title":"Schwinger.EDHoppingMass","text":"EDHoppingMass(lattice) Computes the hopping-type mass term i/2 ∑(-1)^n (χ†ₙ₊₁ χₙ + χ†ₙ₋₁ χₙ) for the Schwinger model.\n\nArguments\n\nlattice::SchwingerLattice: Schwinger model lattice.\n\n\n\n\n\n","category":"function"},{"location":"man/hamiltonian.html#Matrix-product-operator","page":"Hamiltonian","title":"Matrix product operator","text":"","category":"section"},{"location":"man/hamiltonian.html","page":"Hamiltonian","title":"Hamiltonian","text":"Especially for larger lattice sizes where exact diagonalization is infeasible, Schwinger.jl can instead construct a matrix product operator representation of the Hamiltonian. It uses ITensors.jl and ITensorMPS.jl as the backend for all calculations with this form of the Hamiltonian.","category":"page"},{"location":"man/hamiltonian.html","page":"Hamiltonian","title":"Hamiltonian","text":"using Schwinger\nlat = SchwingerLattice{10,1}(periodic = true);\n\n[energygap(EDHamiltonian(lat)), energygap(MPOHamiltonian(lat))]","category":"page"},{"location":"man/hamiltonian.html","page":"Hamiltonian","title":"Hamiltonian","text":"MPOHamiltonian\nMPOGaugeKinetic\nMPOHopping\nMPOMass\nMPOHoppingMass","category":"page"},{"location":"man/hamiltonian.html#Schwinger.MPOHamiltonian","page":"Hamiltonian","title":"Schwinger.MPOHamiltonian","text":"MPOHamiltonian(lattice)\n\nComputes the MPO Hamiltonian for the Schwinger model.\n\nArguments\n\nlattice::SchwingerLattice: Schwinger model lattice.\n\n\n\n\n\n","category":"function"},{"location":"man/hamiltonian.html#Schwinger.MPOGaugeKinetic","page":"Hamiltonian","title":"Schwinger.MPOGaugeKinetic","text":"MPOGaugeKinetic(lattice)\n\nComputes the MPO gauge kinetic operator for the Schwinger model.\n\nArguments\n\nlattice::SchwingerLattice: Schwinger model lattice.\n\n\n\n\n\n","category":"function"},{"location":"man/hamiltonian.html#Schwinger.MPOHopping","page":"Hamiltonian","title":"Schwinger.MPOHopping","text":"MPOHopping(lattice)\n\nComputes the MPO hopping operator for the Schwinger model.\n\nArguments\n\nlattice::SchwingerLattice: Schwinger model lattice.\n\n\n\n\n\n","category":"function"},{"location":"man/hamiltonian.html#Schwinger.MPOMass","page":"Hamiltonian","title":"Schwinger.MPOMass","text":"MPOMass(lattice)\n\nComputes the MPO mass operator for the Schwinger model.\n\nArguments\n\nlattice::SchwingerLattice: Schwinger model lattice.\n\n\n\n\n\n","category":"function"},{"location":"man/hamiltonian.html#Schwinger.MPOHoppingMass","page":"Hamiltonian","title":"Schwinger.MPOHoppingMass","text":"MPOHoppingMass(lattice)\n\nComputes the MPO hopping-mass operator for the Schwinger model.\n\nArguments\n\nlattice::SchwingerLattice: Schwinger model lattice.\n\n\n\n\n\n","category":"function"},{"location":"man/lattices.html#Lattices","page":"Lattices","title":"Lattices","text":"","category":"section"},{"location":"man/lattices.html","page":"Lattices","title":"Lattices","text":"Any calculation in Schwinger.jl will start with a SchwingerLattice. The number of sites N and number of flavors F are type parameters. Other parameters can be specified as keyword arguments:","category":"page"},{"location":"man/lattices.html","page":"Lattices","title":"Lattices","text":"periodic: whether the lattice is periodic\nq: the integer charge of the fermions\nθ2π: the theta-angle (divided by 2pi)\na: the lattice spacing (in coupling units)\nm: the physical mass (in coupling units); the mass shift is applied automatically\nmlat: the mass parameter in the Hamiltonian\nmprime: the coefficient of the hopping-type mass term","category":"page"},{"location":"man/lattices.html","page":"Lattices","title":"Lattices","text":"The sites of the lattice are indexed from 1 to N. The electric field operators are laid out as in the diagram below.","category":"page"},{"location":"man/lattices.html","page":"Lattices","title":"Lattices","text":"(Image: A Schwinger model lattice)","category":"page"},{"location":"man/lattices.html","page":"Lattices","title":"Lattices","text":"Here alpha = 1ldotsF is a flavor index.","category":"page"},{"location":"man/lattices.html","page":"Lattices","title":"Lattices","text":"For details of how these parameters enter into the Hamiltonian, see here.","category":"page"},{"location":"man/lattices.html","page":"Lattices","title":"Lattices","text":"SchwingerLattice","category":"page"},{"location":"man/lattices.html#Schwinger.SchwingerLattice","page":"Lattices","title":"Schwinger.SchwingerLattice","text":"SchwingerLattice{N,F}(;kwargs...)\n\nConstructs a SchwingerLattice for the Schwinger model.\n\nArguments\n\nperiodic::Bool=false: Whether the lattice is periodic.\nq::Int=1: Charge.\nL::Union{Nothing,Real}=nothing: Length of the lattice.\na::Union{Nothing,Real}=nothing: Lattice spacing.\nm::Union{Real,NTuple{N,Real},NTuple{F,Real},NTuple{N,NTuple{F,Real}}=0.: Mass parameter.\nmlat::Union{Nothing,Real,NTuple{N,Real},NTuple{F,Real},NTuple{N,NTuple{F,Real}}=nothing: Local mass parameter.\nmprime::Union{Real,NTuple{N,Real},NTuple{F,Real},NTuple{N,NTuple{F,Real}}=0.: Prime mass parameter.\nθ2π::Union{Real,NTuple{N,Real}}=0.: Theta angle.\n\nReturns\n\nA SchwingerLattice object.\n\n\n\n\n\n","category":"type"}]
}
