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

On an infinite lattice we can also probe the single-particle spectrum using the quasiparticle ansatz. For instance, at $m = 0$ the energy of the lightest quasiparticle will converge to the Schwinger boson mass $1/\sqrt{\pi} \approx 0.564$ as we take $a\to 0$:
```@example infqp
using Schwinger
lat = Lattice(Inf; F = 1, a = 0.5, m = 0)
H = Hamiltonian(lat; backend = :MPSKit)
gs, qp = loweststates(H, 2)
real(energy(qp))   # quasiparticle mass, close to 1/√π
```

The quasiparticle ansatz is a momentum eigenstate. Passing a `momentum` (in units of the coupling $g$, like the rest of the code) to `loweststates` builds the excitation at nonzero momentum — a *moving* quasiparticle (only available on an infinite lattice). We can then use `wavepacket` to build a finite-width, spatially localized packet. Here we build a packet from a quasiparticle of momentum $p = 1.0\,g$, center it on the middle 16 sites of a 64-site window, and plot its energy density relative to the vacuum (dropping the two window-boundary sites, which carry edge artifacts):
```@example infqp
using Plots, LaTeXStrings
qp_moving = loweststates(H, 2; momentum = 1.0)[2]
wp = wavepacket(qp_moving, 64; support = 25:40)
ed = real.(energy_densities(wp))
ed = ed .- ed[5]   # subtract the (uniform) far-field vacuum value
x = (1:64) .* lat.a   # physical position of each site
plot(x[2:63], ed[2:63]; xlabel = L"x", ylabel = "Energy density above vacuum",
     legend = false, title = "Quasiparticle wavepacket (m = 0, ag = 0.5, p = 1.0g)")
```

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