## Installation

To install Schwinger.jl, use the following command in the Julia REPL:

```julia
using Pkg
Pkg.add("https://github.com/srossd/Schwinger.jl")
```

## Key Features

Schwinger.jl provides functions for computing the following properties of the Hamiltonian lattice Schwinger model:

- **Energy levels**: Ground state, energy gap, or higher excited states.
- **Correlators**: Expectation values of various operators in the ground state (or other states).
- **Time evolution**: Action of time evolution on a given initial state.

All of these features are implemented using both exact diagonalization (ED) and matrix product operators/states (MPO).

## Usage Example

Here is a basic example of how to use Schwinger.jl to calculate the average electric field in the one-flavor Schwinger model at $m/g = 0.1$ as a function of $\theta$. Note that the mass shift (see [here](https://arxiv.org/abs/2206.05308)) is included by default.

```@example avgE
using Schwinger
using Base.MathConstants
using Plots

function avgE(θ2π, m, shift = true)
    lat = shift ? SchwingerLattice{10,1}(θ2π = θ2π, m = m, periodic = true) : SchwingerLattice{10,1}(θ2π = θ2π, mlat = m, periodic = true)
    gs = groundstate(EDHamiltonian(lat))
    return real(expectation(EDAverageElectricField(lat), gs))
end

m = 0.05
θ2πs = 0:0.025:1
avgEs_shift = map(x -> avgE(x, m), θ2πs)
avgEs_noshift = map(x -> avgE(x, m, false), θ2πs)

# See eq (24) of https://arxiv.org/abs/2206.05308
perturbative = [(exp(γ)/√(π))*m*sin(2π*θ2π) - (8.9139*exp(2γ)/(4π))*(m^2)*sin(4π*θ2π) for θ2π in θ2πs]

scatter(θ2πs, avgEs_shift, label="With mass shift", xlabel="θ/2π", ylabel="Average Electric Field", legend=:topright, color=:orange)
scatter!(θ2πs, avgEs_noshift, label="Without mass shift", xlabel="θ/2π", ylabel="Average Electric Field", color="lightblue")
plot!(θ2πs, perturbative, label="Perturbative", color=:black)
```

Here is an example of calculating the expectation value of the square of the mean electric field in the two-flavor Schwinger model at $\theta = \pi$, for a four-site periodic lattice, giving a very rough look at the [phase diagram](https://arxiv.org/abs/2305.04437).

```@example twoflavor
using Schwinger
using Plots

function avgE2(m1, m2)
    lat = SchwingerLattice{4,2}(θ2π = 0.5, m = (m1, m2), periodic = true)
    gs = groundstate(EDHamiltonian(lat))
    return real(expectation(EDAverageElectricField(lat; power=2), gs))
end

ms = 0:0.05:1.5
m1s = repeat(reshape(ms, 1, :), length(ms), 1)
m2s = repeat(ms, 1, length(ms))
vals = map(avgE2, m1s, m2s)

contour(ms, ms, vals, fill = true, lw = 0, xlabel = "m₁/g", ylabel = "m₂/g", clabel = "⟨E²⟩")
```