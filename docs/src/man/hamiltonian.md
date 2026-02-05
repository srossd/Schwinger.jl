# Hamiltonian

The Hamiltonian for the lattice Schwinger model is
```math
\begin{split}
H = &\frac{(qg)^2 a}{2}\sum_{n=1}^N \left(L_n + \frac{\theta}{2\pi}\right)^2 - \frac{i}{2a}\sum_{n=1}^N\sum_{\alpha=1}^F \left(\chi^\dagger_{n,\alpha} \chi_{n+1,\alpha} - \chi^\dagger_{n+1,\alpha} \chi_{n,\alpha}\right) \\
&+\underbrace{\left(m - \frac{(qg)^2 F a}{8}\right)}_{m_\text{lat}}\sum_{n=1}^N \sum_{\alpha=1}^F (-1)^n \chi^\dagger_{n,\alpha} \chi_{n,\alpha} + m' \sum_{n=1}^N\sum_{\alpha=1}^F (-1)^{n+1}\left(\chi^\dagger_{n-1,\alpha}\chi_{n,\alpha} + \chi^\dagger_{n+1,\alpha}\chi_{n,\alpha}\right)
\end{split}
```
This is supplemented by the Gauss law
```math
L_n = L_{n-1} + Q_n, \qquad Q_n \equiv q\left(\sum_{\alpha=1}^F \chi^\dagger_{n,\alpha} \chi_{n,\alpha} - \begin{cases} F & n\text{ odd} \\ 0 & n\text{ even} \end{cases}\right).
```

In `Schwinger.jl`, the Hamiltonian can be constructed using three computational backends:
- **Exact Diagonalization (ED)**: Uses sparse matrices for small systems
- **ITensors**: Uses MPO/MPS via ITensors.jl and ITensorMPS.jl
- **MPSKit**: Uses MPO/MPS via MPSKit.jl

## Unified API

The unified API allows you to construct operators using any backend by specifying the `backend` keyword argument:

```julia
# Using symbols
H = Hamiltonian(lattice; backend=:ED)
H = Hamiltonian(lattice; backend=:ITensors)  # Default
H = Hamiltonian(lattice; backend=:MPSKit)

# Using backend types directly
H = Hamiltonian(lattice; backend=EDBackend())
H = Hamiltonian(lattice; backend=ITensorsBackend())
H = Hamiltonian(lattice; backend=MPSKitBackend())
```

You can also set a default backend programmatically:

```julia
set_default_backend(:MPSKit)
H = Hamiltonian(lattice)  # Now uses MPSKit by default
```

## Exact diagonalization

When using exact diagonalization, `Schwinger.jl` constructs a basis of states that diagonalize the operators $\chi^\dagger_{n,\alpha} \chi_{n,\alpha}$ and $L_0$. It then builds a sparse matrix for the Hamiltonian acting on this basis.

```@example ed
using Schwinger
lat = Lattice(12; F = 1);
ham = Hamiltonian(lat; backend=:ED);

ham.matrix
```

You can also use the backend-specific constructor `EDHamiltonian` for backward compatibility.

When $q > 1$, the `universe` (i.e., the allowed values of $L_n$ modulo $q$) can be specified; the default value is 0. When the lattice is periodic, the maximum absolute value of $L_0$ can be set using `L_max`; the default value is 3.

Using [`Arpack.jl`](https://github.com/JuliaLinearAlgebra/Arpack.jl), `Schwinger.jl` can find the lowest eigenstates of a Hamiltonian.
```@example eigs
using Schwinger
lat = Lattice(10; F = 1, periodic = true);
ham = Hamiltonian(lat; backend=:ED);

map(energy, loweststates(ham, 5))
```

```@docs
EDHamiltonian
EDGaugeKinetic
EDHopping
EDMass
EDHoppingMass
```

## Matrix product operators

For larger lattice sizes where exact diagonalization is infeasible, `Schwinger.jl` can construct a matrix product operator representation of the Hamiltonian using two backends:

- **ITensors**: Uses [`ITensors.jl`](https://github.com/ITensor/ITensors.jl) and [`ITensorMPS.jl`](https://github.com/ITensor/ITensorMPS.jl)
- **MPSKit**: Uses [`MPSKit.jl`](https://github.com/maartenvd/MPSKit.jl)

```@example mpo
using Schwinger
lat = Lattice(10; F = 1);

# Compare energy gaps across backends
[
    energygap(Hamiltonian(lat; backend=:ED)),
    energygap(Hamiltonian(lat; backend=:ITensors)),
    energygap(Hamiltonian(lat; backend=:MPSKit))
]
```

### ITensors Backend

```@docs
ITensorHamiltonian
ITensorGaugeKinetic
ITensorHopping
ITensorMass
ITensorHoppingMass
```

### MPSKit Backend

```@docs
MPSKitHamiltonian
MPSKitGaugeKinetic
MPSKitHopping
MPSKitMass
MPSKitHoppingMass
```