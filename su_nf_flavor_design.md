# Design: optional SU(Nf) flavor symmetry for Nf equal-mass flavors

> **STATUS: IMPLEMENTED.** Enable with `Lattice(N; F=Nf, m=m0, flavor_sym=true)`; the MPSKit
> backend then builds the `U(1)×SU(F)` Hamiltonian (`src/flavor_symmetry.jl`), and
> `loweststates(Hamiltonian(lat; backend=:MPSKit), n)` returns its states. Validated in
> `test/flavor_symmetry.jl` (ground state == Schwinger.jl ED to <1e-6, SU(3) builds, argument
> checks). Requires equal flavor masses, `F ≥ 2`, non-periodic, `mprime = 0`, no defects.
> The prototype/validation scripts in `examples/su_nf_flavor_*.jl` remain as standalone references.

---


## Motivation

Schwinger.jl currently represents `F` flavors as `F` **separate U(1) sites** per staggered site: the
MPSKit unit cell is `2F` sites, each a 2-dim occupation space, and the hopping MPO has bond dimension
`2 + 2F` (see `matrices_hopping`, `get_mpskit_spaces` in `src/hamiltonian.jl`). The only symmetry kept
is the U(1) gauge charge; the flavors are distinguishable, so nothing enforces the flavor symmetry that
is exact when the masses are equal.

When all `Nf` flavors have **equal mass**, the theory has a global **SU(Nf) flavor symmetry**. The
reference `two_flavor_equal_setup.jl` shows, for `Nf=2`, how to gauge that symmetry into the MPS by
using the product category `U(1) × SU(2)`: all flavors live on **one** physical site whose Hilbert
space is organized into flavor multiplets. This document generalizes that to arbitrary `Nf` and
proposes how to fold it into Schwinger.jl as an opt-in.

The payoff is large and is **already demonstrated** by the validated prototype
(`examples/su_nf_flavor_prototype.jl`):

| quantity                    | current (U(1), F flavors) | SU(Nf)-symmetric        |
|-----------------------------|---------------------------|-------------------------|
| MPSKit unit cell            | `2F` sites                | **2 sites** (any Nf)    |
| on-site Hilbert space       | 2 (per mode)              | `2^Nf`, split into `Nf+1` flavor sectors `Λᵏ` |
| hopping MPO bond dim        | `2 + 2F`                  | **4** (any Nf)          |
| symmetry                    | U(1)                      | U(1) × SU(Nf)           |
| excitations labeled by      | charge only               | charge **and flavor irrep** |

Because the SU(Nf) multiplet structure is made *exact* (a quantum number) rather than *learned* by the
bond dimension, the MPS is dramatically more efficient, and — the real physics win — excited states can
be targeted directly by their flavor irrep (flavor-singlet vs. adjoint mesons, etc.).

---

## 1. The physics being encoded

On one staggered site there are `Nf` fermion modes. Their Fock space is `Λ*(ℂ^Nf) = ⊕ₖ Λᵏ(ℂ^Nf)`,
where `ℂ^Nf` is the SU(Nf) **fundamental**. The `k`-fermion sector is the `k`-th antisymmetric power
`Λᵏ` (fermions ⇒ antisymmetric), which is an **irreducible** SU(Nf) representation of dimension `C(Nf,k)`:

- `Λ⁰` = singlet (0 fermions), `Λ¹` = fundamental (1 fermion), `Λ²` = antisym 2-tensor, …, `Λ^Nf` = singlet.
- `Σₖ dim Λᵏ = 2^Nf` = full Fock space. (Validated: `Nf=2 → [1,2,1]`, `Nf=3 → [1,3,3,1]`, `Nf=4 → [1,4,6,4,1]`.)

The U(1) gauge charge of the `k`-fermion sector is `q·k` (even site) or `q·(k−Nf)` (odd/staggered site,
particle–hole shifted). So the physical site is the graded space

```
Pₑ = ⊕ₖ ( charge = q·k,      Λᵏ )        (even)
P_o = ⊕ₖ ( charge = q·(k−Nf), Λᵏ )        (odd)
```

matching the reference's `Ve = ((0,sing),(1,fund),(2,sing))`, `Vo = ((−2,sing),(−1,fund),(0,sing))` at
`Nf=2, q=1`.

### The one nontrivial ingredient: the creation operator

`χ†_f` (flavor `f`) transforms as the fundamental, so as a tensor operator it carries a `VS = (charge q, fund)`
auxiliary leg: `χ† : VS ⊗ Pₖ ← Pₖ`, mapping `Λᵏ → Λᵏ⁺¹`. In TensorKit the block structure already
carries the Clebsch–Gordan coefficients, so `χ†` is fixed by one reduced matrix element `rₖ` per
transition `Λᵏ → Λᵏ⁺¹`.

**Key result (derived, and validated numerically for Nf=2,3,4):**

```
rₖ = √(k+1),   independent of Nf.
```

Derivation: require canonical fermions. Tracing the anticommutator over the `k`-fermion sector,
`Σ_f χ†_f χ_f = k` gives `rₖ₋₁² · dim Λᵏ = k · dim Λᵏ ⇒ rₖ₋₁ = √k`; the other anticommutator
`Σ_f χ_f χ†_f = Nf − k` gives `rₖ² · dim Λᵏ⁺¹ = (Nf−k)·dim Λᵏ`, and since
`dim Λᵏ⁺¹/dim Λᵏ = (Nf−k)/(k+1)`, again `rₖ = √(k+1)`. Both channels agree for every `Nf`.

**Why this is clean:** because the physical space `P` contains *only* the antisymmetric irreps, the
operator `χ† : VS ⊗ P ← P` automatically projects onto the fermionic channel — the symmetric part of
`fund ⊗ Λᵏ` simply has no block in `P`. No manual antisymmetrization is needed. (Validated: with
`rₖ=√(k+1)`, the number operator `χ†χ` has eigenvalues `0,1,…,Nf` exactly, for Nf=2,3,4.)

Everything else is SU(Nf)-blind and identical to the reference / current code:
- **Jordan–Wigner string** `F = (−1)^{#ferm}` (diagonal in `k`).
- **Staggered mass** `= (−1)ⁿ(#ferm − Nf/2)`: `Mₑ(k)=k−Nf/2`, `M_o(k)=Nf/2−k`.
- **DKPZ shift** `mlat = m − Nf·q²·a/8` (⇒ `m − a/4` at Nf=2, matching the reference).
- **Electric energy** `(a/2)(L+θ/2π)²` — a LEMPO `link_fct` of the U(1) charge only (SU(Nf)-blind).
- **Hopping MPO**: bond dim 4 (levels: start / χ†-forward / χ-backward / done), the flavor index bundled
  in the `VS=fund` bond ⇒ **Nf-independent structure**.

---

## 2. Validation status

Prototype: `examples/su_nf_flavor_prototype.jl` (infinite builder). Finite cross-check:
`examples/su_nf_flavor_finite_check.jl`.

- ✅ SU(Nf) spaces build for Nf=2,3,4; `dim P = 2^Nf`, split into `Nf+1` sectors `[1,2,1] / [1,3,3,1] / [1,4,6,4,1]`.
- ✅ `rₖ=√(k+1)` gives a **canonical** fermion number operator (`χ†χ = 0,1,…,Nf`) for Nf=2,3,4.
- ✅ The full LEMPO assembles for Nf=2,3,4 with **hopping bond 4, cell 2 sites**.
- ✅ **Matches Schwinger.jl's existing code, machine precision** (`examples/su_nf_flavor_schwinger_check.jl`).
  On a finite N=4-site, F=2 chain (`a=1, m=0.5, θ=0`), the SU(2)-symmetric prototype gives
  `E = −2.0518252582`, agreeing with **Schwinger.jl's exact diagonalization to 6.7e-15** and with its
  MPSKit U(1) backend to 4.7e-10 (DMRG tolerance). This is the strong check: the two use *entirely
  different representations* — the prototype puts both flavors on one site as SU(2) multiplets,
  Schwinger.jl uses two separate U(1) sites — so it confirms the conventions agree (DKPZ shift,
  staggering, electric term, boundaries), not just internal consistency. (Open boundaries fix the
  electric field by Gauss law ⇒ `L_max=0` ⇒ ED is *exact*, no cutoff.)
- ✅ **Encoding cross-check** (`examples/su_nf_flavor_finite_check.jl`): on a finite 8-site chain, the
  same construction with `SUNIrrep{2}` vs `SU2Irrep` flavor irreps agrees to `|ΔE| = 4.3e-14`,
  isolating the SU(N) encoding itself.
- ✅ Nf=3 ground state converges (new; no closed-form reference).

A finite chain is used deliberately for all energy checks: infinite VUMPS is flaky here (LAPACK SVD
failures at small bond dim, slow at large), whereas finite DMRG / ED on a few unit cells is robust and,
for open boundaries, exact.

Two TensorKit gotchas the prototype resolves (worth carrying into the package):
1. **Use the concrete sector type**: `U1Irrep ⊠ typeof(lam(Nf,1))` (e.g. `SUNIrrep{Nf,·}`, or `SU3Irrep`),
   *not* the UnionAll `U1Irrep ⊠ SUNIrrep{Nf}` — the latter makes `@tensor conj(...)` fail with a
   `SUNIrrep{Nf,1}` → `SUNIrrep{Nf}` convert error.
2. **Real (Float64) tensors**, as in the reference: the operators are real in this basis, and a
   ComplexF64 `Elt` union breaks the scalar `W[1,1]=1.0` MPO assignment.

---

## 3. Proposed Schwinger.jl integration

### 3.1 API

Add an opt-in on `Lattice`:

```julia
Lattice(N; F = Nf, m = m0, flavor_sym = true, ...)   # SU(Nf)-symmetric MPSKit backend
```

`flavor_sym::Bool = false` (default keeps today's behavior exactly). When `true`:

- **Preconditions** (throw a clear error otherwise):
  - `F ≥ 2` and **all flavor masses equal** (`all(mlat[s] .== mlat[s][1])` per site; the symmetry is
    exact only at equal mass — this is the whole premise).
  - MPSKit backend only. ED and ITensors use the occupation basis and cannot carry SU(Nf) blocks; error
    (or fall back to the U(1) path with a warning) for those backends.
  - `q` and `θ2π` are unrestricted (electric term is SU(Nf)-blind).

### 3.2 Where the code goes

The SU(Nf) construction changes the *site content* (1 site/staggered-site vs `F`), so it is a parallel
MPSKit code path, dispatched on `lattice.flavor_sym`:

| current function            | flavor-sym variant                              |
|-----------------------------|-------------------------------------------------|
| `get_mpskit_spaces`         | build `Pₑ,P_o` from `Λᵏ` irreps (§1)            |
| `MPSKitGaugeKinetic`        | LEMPO with `link_fct` on U(1) charge, one per site |
| `matrices_hopping`/`Mass`   | the bond-4 hopping + mass MPO from §1 operators |

Concretely, port `examples/su_nf_flavor_prototype.jl`'s `build_su_nf` behind a
`get_mpskit_spaces(lattice)` / `MPSKitHamiltonian(lattice)` branch:

```julia
function get_mpskit_spaces(lattice::Lattice)
    lattice.flavor_sym && return _mpskit_spaces_flavorsym(lattice)   # §1 Pₑ,P_o, unit cell 2
    ...existing U(1) code...
end
```

`_flavorsym` helpers need only `(Nf=F, a, mlat, q, θ2π, N)`; the symmetry group is
`U1Irrep ⊠ typeof(lam(Nf,1))`. For `Nf=2` you may special-case `SU2Irrep` (built into TensorKit, no
`SUNRepresentations` dependency); for `Nf ≥ 3` use `SUNRepresentations.SUNIrrep{Nf}`.

### 3.3 States, excitations, observables

- **Excitations by flavor irrep** (IMPLEMENTED): `loweststates(H, n; momentum, flavor_irrep=R)` targets
  the excitations at flavor sector `(0, R)` (infinite lattice / `QuasiparticleAnsatz`). Helpers
  `flavor_singlet(lat)` / `flavor_adjoint(lat)` give the two meson channels (fund ⊗ fund̄ = singlet ⊕
  adjoint). Wired via `sector = _flavorsym_excitation_sector(lat, R)` in `src/states.jl`.
- **Observables** (IMPLEMENTED): `occupation(s)` collapse the flavor array to **one number per staggered
  site** — the common per-flavor value `⟨n_f⟩` — and `charge(s)`, `electricfield(s)`, `scalardensit*`,
  `scalar`, `entanglement(s)` follow, each matching the F-site ED computation to DMRG tolerance
  (`test/flavor_symmetry.jl`). Per-flavor-resolved quantities are, by construction, not accessible
  (flavor is summed); resolve them by SU(F) irrep instead. `energy_density` (per-*site* decomposition)
  is not yet flavor-aware — it errors and points to `energy(state)` for the total.
- **Dependency**: `SUNRepresentations` added to the project.

### 3.4 Interaction with the GFF program

The gravitational form factor machinery (`qp_matrix_element`, `examples/gff`) transfers directly: `T⁰⁰`
is built as a `FiniteLEMPOHamiltonian` from the same operators, `Σⱼ T⁰⁰ⱼ = H` still gives `G(0)=1`, and
the forward oracle still validates. The new capability is that the meson whose `G(t)` you measure can be
selected by flavor irrep — e.g. the GFF of the flavor-singlet vs. flavor-adjoint meson.

---

## 4. Summary

The generalization is clean because a single universal rule, `rₖ = √(k+1)`, turns the reference's SU(2)
construction into an arbitrary-`Nf` one, and TensorKit's antisymmetric-irrep blocks handle the rest. The
recommended path: (i) port the validated `build_su_nf` behind a `flavor_sym` branch of the MPSKit
backend, guarded by an equal-mass check; (ii) special-case `SU2Irrep` for `Nf=2`, `SUNIrrep{Nf}` above;
(iii) add flavor-irrep-resolved observables. Cost per staggered site rises as `2^Nf` in raw Hilbert-space
dimension but is offset by the exact multiplet blocking and the fixed bond-dim-4 hopping, and it unlocks
flavor-resolved spectroscopy and form factors that the U(1) representation cannot express.
