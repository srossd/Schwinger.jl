# Operators

Operators in `Schwinger.jl` are represented by the abstract type `SchwingerOperator`, with two descendants:
- `EDOperator`: an operator represented as a matrix for its action on a basis of `SchwingerBasisState`s
- `MPOOperator`: a matrix product operator, stored as an `MPO` object using [`ITensorMPS.jl`](https://github.com/ITensor/ITensorMPS.jl)

We can evaluate the expectation of an operator in a state:
```@example avgE
using Schwinger
lat = SchwingerLattice{10,1}()
gs = groundstate(EDHamiltonian(lat))

sum(electricfields(gs))/10, expectation(EDAverageElectricField(lat), gs)
```

This can also be carried out manually by acting on the state with the operator:
```@example avgE
using LinearAlgebra
dot(gs, EDAverageElectricField(lat) * gs)
```

```@docs
expectation
act
```

## Exact diagonalization

```@docs
EDWilsonLoop
EDWilsonLine
EDAverageElectricField
```

## Matrix product operators

```@docs
MPOWilsonLoop
MPOWilsonLine
MPOAverageElectricField
```
