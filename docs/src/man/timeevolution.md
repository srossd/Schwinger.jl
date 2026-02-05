# Time evolution

`Schwinger.jl` supports time-evolving states using all three backends:
- **ED**: Uses Krylov methods from [`KrylovKit.jl`](https://github.com/Jutho/KrylovKit.jl)
- **ITensors**: Uses the time-dependent variational principle (TDVP) from [`ITensorMPS.jl`](https://github.com/ITensor/ITensorMPS.jl)
- **MPSKit**: Uses TDVP from [`MPSKit.jl`](https://github.com/maartenvd/MPSKit.jl)

The `evolve` function evolves a state forwards in time, and monitors any given observables. It returns a final state along with a `DataFrame` of the observables. For example, here is a simulation of flux unwinding.
```@example time
using Schwinger, Plots

lat = Lattice(10; F = 1, L = 2, periodic = true)
gs = groundstate(Hamiltonian(lat; backend=:ED))

_, df = evolve(WilsonLoop(lat; backend=:ED) * gs, 15; 
    nsteps = 30, 
    observable = (ψ, t) -> sum(electricfields(ψ))/10
)

scatter(df.time, df.observable, xlabel = "gt", ylabel = "Average electric field", label = "Schwinger.jl")
plot!(0:.1:15, [cos(t/√(π)) for t in 0:.1:15], label = "Exact")
```

```@docs
evolve
```