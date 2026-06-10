# Defects

A `DefectCharge` is a static (external) charge inserted into the lattice. It enters
only through Gauss's law: a charge `q` placed on the left of link `n` shifts the
electric flux on every link `≥ n` by `+q`. A pair of opposite charges therefore
sources a flux string between them, against which the dynamical matter rearranges.

Defects are passed to `Hamiltonian` with the `defects` keyword and work with every
backend:
- **ED** / **ITensors** realise them as a `θ2π` step (no extra site);
- **MPSKit** inserts a genuine extra lattice site carrying the charge, inert under
  the matter Hamiltonian but contributing to the gauge accumulation. The inserted
  site is invisible to `Schwinger.jl` observables. (This is overkill for an abelian theory, but it generalizes straightforwardly to non-abelian charges.)

## Example: static charges on a 32-site lattice

Here we place charges `±1` at sites 8 and 24 of a 32-site lattice with `q = 4` and
`m = 0`, compute the ground state with the ITensors backend, and plot the scalar and pseudoscalar densities. The chiral condensate rearranges across the charged region (see [here](https://arxiv.org/pdf/2210.04237)).

```@example defects
using Schwinger, Plots

lat  = Lattice(32; F = 1, q = 4, m = 0)
defs = [DefectCharge(8, 1), DefectCharge(24, -1)]
gs   = groundstate(Hamiltonian(lat; backend = :ITensors, defects = defs))

S = real.(scalardensities(gs))
P = real.(pseudoscalardensities(gs))

p = plot(1:length(S), S, label = "scalar",
          xlabel = "Site", ylabel = "Density", title = "(Pseudo)scalar density")
plot!(p, 1:length(P), P, label = "pseudoscalar")
vline!(p, [8, 24], ls = :dash, color = :gray, label = "")

p
```

## Manipulating defects

A state remembers the defects it was built with, and the static charges can be moved
or added/removed after the fact. Relocating a charge leaves the matter wavefunction
untouched (so `occupations` are preserved while `electricfields` shift); inserting or
removing one changes the total charge and is the way to begin or end a Wilson line.

```@docs
DefectCharge
defects
translate_defect
insert_defect
remove_defect
```
