# States

Lattice states in `Schwinger.jl` are represented by the abstract type `SchwingerState`, with four concrete types: 
- `BasisState`: a state specified by the eigenvalues of occupation operators $\chi^\dagger_{n,\alpha}\chi_{n,\alpha}$ and $L_0$
- `EDState`: a linear combination of `BasisState`s
- `ITensorState`: a matrix product state using [`ITensorMPS.jl`](https://github.com/ITensor/ITensorMPS.jl)
- `MPSKitState`: a matrix product state using [`MPSKit.jl`](https://github.com/maartenvd/MPSKit.jl)

Given a state, we can find the expectation values of the occupation operators and electric field operators:
```@example sc
using Schwinger
lat = Lattice(6; F = 1, a = 10) # towards the lattice strong coupling limit ga -> infty
gs = groundstate(Hamiltonian(lat; backend=:ED))

occupations(gs), electricfields(gs)
```

We can also evaluate the entanglement entropies of each bisection of the lattice:
```@example entanglement
using Schwinger
lat = Lattice(20; F = 1, )
gs = groundstate(Hamiltonian(lat; backend=:ITensors))
entanglements(gs)
```

When using an infinite lattice with `MPSKit`, we can compute these quantities in the thermodynamic limit.
```@example infentanglement
using Schwinger
lat = Lattice(Inf; F = 1, )
gs = groundstate(Hamiltonian(lat; backend=:MPSKit))
entanglements(gs)
```
The two numbers here correspond to the two possible bisections of the two-site unit cell (corresponding to an odd and an site from the middle of a large finite lattice).

Several other useful functions are detailed below.

```@docs
SchwingerState
BasisState
EDState
ITensorState
MPSKitState
occupation
occupations
charge
charges
electricfield
electricfields
entanglement
entanglements
energy
L₀
scalar
scalardensity
scalardensities
pseudoscalar
pseudoscalardensity
pseudoscalardensities
```