# States

Lattice states in `Schwinger.jl` are represented by the abstract type `SchwingerState`, with three descendants: 
- `SchwingerBasisState`: a state specified by the eigenvalues of occupation operators $\chi^\dagger_{n,\alpha}\chi_{n,\alpha}$ and $L_0$
- `SchwingerEDState`: a linear combination of `SchwingerBasisState`s
- `SchwingerMPS`: a matrix product state, stored as an `MPS` object using [`ITensorMPS.jl`](https://github.com/ITensor/ITensorMPS.jl)

Given a state, we can find the expectation values of the occupation operators and electric field operators:
```@example sc
using Schwinger
lat = SchwingerLattice{6,1}(a = 10) # towards the lattice strong coupling limit ga -> infty
gs = groundstate(EDHamiltonian(lat))

occupations(gs), electricfields(gs)
```

With a `SchwingerMPS` state, we can evaluate the entanglement entropies of each bisection of the lattice:
```@example entanglement
using Schwinger
lat = SchwingerLattice{20,1}()
gs = groundstate(MPOHamiltonian(lat))
entanglements(gs)
```

Several other useful functions are detailed below.

```@docs
SchwingerState
SchwingerBasisState
SchwingerEDState
SchwingerMPS
occupations
charges
electricfields
entanglements
energy
L₀
scalarvev
pseudoscalarvev
```