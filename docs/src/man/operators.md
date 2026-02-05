# Operators

Operators in `Schwinger.jl` are represented by the abstract type `SchwingerOperator`, with three concrete types:
- `EDOperator`: an operator represented as a matrix for its action on a basis of `BasisState`s
- `ITensorOperator`: a matrix product operator using [`ITensorMPS.jl`](https://github.com/ITensor/ITensorMPS.jl)
- `MPSKitOperator`: a matrix product operator using [`MPSKit.jl`](https://github.com/maartenvd/MPSKit.jl)

## Unified API

Like the Hamiltonian, other operators can be constructed using the unified API with the `backend` keyword:

```julia
W = WilsonLoop(lattice; backend=:ED)
W = WilsonLoop(lattice; backend=:ITensors)
W = WilsonLoop(lattice; backend=:MPSKit)
```

We can evaluate the expectation of an operator in a state:
```@example avgE
using Schwinger
lat = Lattice(10; F = 1, )
gs = groundstate(Hamiltonian(lat; backend=:ED))

sum(electricfields(gs))/10, expectation(AverageElectricField(lat; backend=:ED), gs)
```

This can also be carried out manually by acting on the state with the operator:
```@example avgE
using LinearAlgebra
dot(gs, AverageElectricField(lat; backend=:ED) * gs)
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

## ITensors operators

```@docs
ITensorWilsonLoop
ITensorWilsonLine
ITensorAverageElectricField
```

## MPSKit operators

```@docs
MPSKitWilsonLoop
MPSKitWilsonLine
MPSKitAverageElectricField
```
