# Hamiltonian

The Hamiltonian for the lattice Schwinger model is
```math
\begin{split}
H = &\frac{(qg)^2 a}{2}\sum_{n=1}^N \left(L_n + \frac{\theta}{2\pi}\right)^2 - \frac{i}{2a}\sum_{n=1}^N\sum_{\alpha=1}^F \left(\chi^\dagger_{n,\alpha} \chi_{n+1,\alpha} - \chi^\dagger_{n+1,\alpha} \chi_{n,\alpha}\right) \\
&+\underbrace{\left(m - \frac{(qg)^2 F a}{8}\right)}_{m_\text{lat}}\sum_{n=1}^N \sum_{\alpha=1}^F (-1)^n \chi^\dagger_{n,\alpha} \chi_{n,\alpha} + m' \sum_{n=1}^N\sum_{\alpha=1}^F (-1)^j\left(\chi^\dagger_{n-1,\alpha}\chi_{n,\alpha} + \chi^\dagger_{n+1,\alpha}\chi_{n,\alpha}\right)
\end{split}
```
This is supplemented by the Gauss law
```math
L_n = L_{n-1} + Q_n, \qquad Q_n \equiv q\left(\sum_{\alpha=1}^F \chi^\dagger_{n,\alpha} \chi_{n,\alpha} - \begin{cases} F & n\text{ odd} \\ 0 & n\text{ even} \end{cases}\right).
```

In `Schwinger.jl` this be constructed using two strategies, exact diagonalization (`ED`) or matrix product operators (`MPO`).

## Exact diagonalization

When using exact diagonalization, `Schwinger.jl` constructs a basis of states that diagonalize the operators $\chi^\dagger_{n,\alpha} \chi_{n,\alpha}$ and $L_0$. It then builds a sparse matrix for the Hamiltonian acting on this basis.

```@example ed
using Schwinger
lat = SchwingerLattice{12,1}();
ham = EDHamiltonian(lat);

ham.matrix
```

When $q > 1$, the `universe` (i.e., the allowed values of $L_n$ modulo $q$) can be specified; the default value is 0. When the lattice is periodic, the maximum absolute value of $L_0$ can be set using `L_max`; the default value is 3.

Using [`Arpack.jl`](https://github.com/JuliaLinearAlgebra/Arpack.jl), `Schwinger.jl` can find the lowest eigenstates of a Hamiltonian.
```@example eigs
using Schwinger
lat = SchwingerLattice{10,1}(periodic = true);
ham = EDHamiltonian(lat);

map(energy, loweststates(ham, 5))
```

```@docs
EDHamiltonian
EDGaugeKinetic
EDHopping
EDMass
EDHoppingMass
```

## Matrix product operator

Especially for larger lattice sizes where exact diagonalization is infeasible, `Schwinger.jl` can instead construct a matrix product operator representation of the Hamiltonian. It uses [`ITensors.jl`](https://github.com/ITensor/ITensors.jl) and [`ITensorMPS.jl`](https://github.com/ITensor/ITensorMPS.jl) as the backend for all calculations with this form of the Hamiltonian.

```@example mpo
using Schwinger
lat = SchwingerLattice{10,1}(periodic = true);

[energygap(EDHamiltonian(lat)), energygap(MPOHamiltonian(lat))]
```

```@docs
MPOHamiltonian
MPOGaugeKinetic
MPOHopping
MPOMass
MPOHoppingMass
```