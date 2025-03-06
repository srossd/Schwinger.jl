# Time evolution

`Schwinger.jl` supports time-evolving states. In the exact diagonalization framework, this is accomplished using Krylov methods from [`KrylovKit.jl`](https://github.com/Jutho/KrylovKit.jl). With matrix product operators, this is accomplished using the time-dependent variational principle algorithm (TDVP), as implemented in [`ITensorMPS.jl`](https://github.com/ITensor/ITensorMPS.jl).

The `evolve` function evolves a state forwards in time, and monitors any given observables. It returns a final state along with a `DataFrame` of the observables. For example, here is a simulation of flux unwinding.
```@example time
using Schwinger, Plots

lat = SchwingerLattice{10,1}(L = 2, periodic = true)
gs = groundstate(EDHamiltonian(lat))

_, df = evolve(EDWilsonLoop(lat) * gs, 10; nsteps = 30, observable = (ψ, t) -> sum(electricfields(ψ))/10)

scatter!(df.time, df.observable, xlabel = "gt", ylabel = "Average electric field")
```

```@docs
evolve
```