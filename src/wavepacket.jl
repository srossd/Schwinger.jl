# Quasiparticle wavepacket construction on an infinite (WindowMPS) background:
# the reflection-symmetric Bloch gauge, the domain-wall pseudo-inverse, and the
# automaton that builds a (multi-)wavepacket finite window. Included after states.jl.

"""
    reflection_symmetric_gauge(qp::MPSKit.QP; tol=1e-12, maxiter=1000)

Re-gauge a quasiparticle's B tensors into the reflection-symmetric gauge of Van Damme
et al. ([arXiv:2012.07243](https://arxiv.org/abs/2012.07243), eq. C17), returning a
vector of the re-gauged (full) B tensors, one per unit-cell site.

The Bloch ansatz `|Φ_p(B)⟩` is invariant under `B[i] → B[i] + AL[i]·Y[i+1] -
e^{-ip}·Y[i]·AR[i]`. Left/right gauge fixing breaks this freedom asymmetrically; the
reflection-symmetric choice instead fixes it by minimising

    Σ_i ‖Σ_s B[i]^s ⊗ (AL[i]^s)*‖² + ‖Σ_s B[i]^s ⊗ (AR[i]^s)*‖²

which treats the two spatial directions on equal footing. This gives localized states
(and wavepackets) that are more symmetric and lie more accurately within the
single-quasiparticle subspace. The minimisation is a linear least-squares problem,
solved here via the normal equations.
"""
function reflection_symmetric_gauge(qp::MPSKit.QP; tol::Real = 1e-12, maxiter::Int = 1000)
    L  = length(qp)
    p  = qp.momentum
    AL = [qp.left_gs.AL[i]  for i in 1:L]
    AR = [qp.right_gs.AR[i] for i in 1:L]
    B0 = [qp[i] for i in 1:L]

    # gauge transformation ΔB[i] = AL[i]·Y[i+1] - e^{-ip}·Y[i]·AR[i] (leaves |Φ_p⟩ invariant)
    dB(Y) = map(1:L) do i
        ip1 = mod1(i + 1, L)
        @plansor t[-1 -2; -3 -4] := AL[i][-1 -2; 1] * Y[ip1][1; -3 -4]
        @plansor t[-1 -2; -3 -4] -= exp(-im * p) * Y[i][-1; -3 1] * AR[i][1 -2; -4]
        t
    end
    # its adjoint (B-space → Y-space)
    dBadj(D) = map(1:L) do j
        jm1 = mod1(j - 1, L)
        @plansor za[-1; -2 -3] := conj(AL[jm1][1 2; -1]) * D[jm1][1 2; -2 -3]
        @plansor zb[-1; -2 -3] := D[j][-1 2; -2 1] * conj(AR[j][-3 2; 1])
        za - exp(im * p) * zb
    end
    # 𝒦_A(B)[i]: Hermitian PSD operator with ⟨B, 𝒦_A B⟩ = ‖Σ_s A^s ⊗ (B^s)*‖²
    Kop(B, A) = map(1:L) do i
        @tensor G[-1; -2] := A[i][1 -1; 2] * conj(A[i][1 -2; 2])
        @tensor KB[-1 -2; -3 -4] := G[-2; 1] * B[i][-1 1; -3 -4]
        KB
    end

    # minimise f(Y) = ⟨B', (𝒦_AL+𝒦_AR) B'⟩ with B' = B0 + dB(Y): normal equations H Y = rhs
    H(Y)  = dBadj([Kop(dB(Y), AL)[i] + Kop(dB(Y), AR)[i] for i in 1:L])
    rhs   = dBadj([-(Kop(B0, AL)[i] + Kop(B0, AR)[i]) for i in 1:L])
    Y0    = [zero(y) for y in dBadj(B0)]
    Ysol, info = linsolve(H, rhs, Y0; maxiter = maxiter, tol = tol)
    info.converged == 0 &&
        @warn "reflection_symmetric_gauge: linsolve did not converge (normres = $(info.normres))"
    ΔB = dB(Ysol)
    return [B0[i] + ΔB[i] for i in 1:L]
end

# Truncated pseudo-inverse of a ground-state center matrix, used as an AR→AL domain wall.
# The Schmidt spectrum decays steeply (e.g. 0.91, 0.41, 0.04, …, 1e-13) over a bond space
# far larger than the effective rank, so a plain `inv` is catastrophically ill-conditioned
# (‖C⁻¹‖ ~ 1e12) and amplifies the ground state's finite convergence error (~1e-8) into an
# O(1) defect at the wall. We invert only singular values above `rtol`·σmax and zero the
# rest, keeping ‖C⁻¹‖ ≲ 1/rtol·σmax; the dropped directions carry weight ~rtol², so the
# wall is accurate to that order while staying well-conditioned.
function _trunc_pinv(C::AbstractTensorMap; rtol::Real = 1e-2)
    # Threshold relative to the GLOBAL largest singular value across all charge blocks.
    # (A per-block maximum would set a tiny threshold in sub-dominant sectors and invert
    # their near-zero Schmidt values, blowing the state norm up to ~1e8.)
    smax = 0.0
    for (_, b) in blocks(C)
        sv = LinearAlgebra.svdvals(Matrix(b))
        isempty(sv) || (smax = max(smax, maximum(sv)))
    end
    cutoff = rtol * smax
    Ci = zero(C')
    for (c, b) in blocks(C)
        F = LinearAlgebra.svd(Matrix(b))
        Sinv = [s > cutoff ? inv(s) : zero(s) for s in F.S]
        block(Ci, c) .= F.Vt' * LinearAlgebra.Diagonal(Sinv) * F.U'
    end
    return Ci
end

"""
    wavepacket(state::MPSKitQPState, W::Int; support=1:W, weights=nothing)

Build a wavepacket `MPSKitState` from an infinite quasiparticle state.

The total window has `W` sites. The B tensor is summed only over the sites in
`support` (a `UnitRange`, default `1:W`), so sites outside `support` are pure
ground state. This lets you build a spatially localized wavepacket inside a
larger window — e.g. `W=100, support=10:20`.

`weights` gives one complex amplitude per site in `support`. If `nothing`, the
default is `exp(im·p·k)` for each site k ∈ support (the QP momentum phase),
giving a position-space truncation of a momentum eigenstate.

`gauge` selects how the residual Bloch gauge freedom on the B tensors is fixed:
`:symmetric` (default) uses the reflection-symmetric gauge of arXiv:2012.07243
(see [`reflection_symmetric_gauge`](@ref)), giving more symmetric, better-localized
wavepackets, with every QP centred identically on a mixed-canonical background;
`:left` uses MPSKit's left gauge.

For multiple QPs the vacuum returns AR→AL between them through a ground-state
center-matrix domain wall, built with a truncated pseudo-inverse: `wall_rtol`
(default `1e-2`, relative to the largest Schmidt value) sets the truncation. A
poorly-converged ground state has poorly-determined small Schmidt values that the
inverse amplifies into a defect at the wall; raise `wall_rtol` (e.g. `0.1`) for a
cruder ground state, or converge it better.

By default (`window = false`) the result is an `MPSKitState` wrapping a plain
`FiniteMPS` over the `W` window sites. With `window = true` it instead wraps a
`WindowMPS` whose left/right environments are the QP's `left_gs`/`right_gs`; the two
give identical local observables, but the `WindowMPS` carries the explicit infinite
environments (use it for real-time evolution that should keep the infinite background
fixed). The construction follows MPSKit's `convert(FiniteMPS, v::QP{FiniteMPS})`,
adapted for the infinite case.

!!! note
    Assumes a non-topological excitation (trivial auxiliary sector).
"""
function wavepacket(state::MPSKitQPState, W::Int;
                      support::UnitRange{Int}=1:W,
                      sigma::Union{Nothing,Real}=nothing,
                      center::Union{Nothing,Real}=nothing,
                      weights::Union{Nothing,AbstractVector}=nothing,
                      gauge::Symbol=:symmetric,
                      window::Bool=false,
                      wall_rtol::Real=1e-2)
    return wavepacket([state], W;
                      supports=[support], sigmas=[sigma],
                      centers=[center], weights=[weights], gauge=gauge, window=window,
                      wall_rtol=wall_rtol)
end

function wavepacket(states::AbstractVector{<:MPSKitQPState}, W::Int;
                      supports::AbstractVector=fill(1:W, length(states)),
                      sigmas::AbstractVector=fill(nothing, length(states)),
                      centers::AbstractVector=fill(nothing, length(states)),
                      weights::AbstractVector=fill(nothing, length(states)),
                      gauge::Symbol=:symmetric,
                      window::Bool=false,
                      wall_rtol::Real=1e-2)
    n = length(states)
    n > 0 || throw(ArgumentError("states must be non-empty"))
    length(supports) == n && length(sigmas) == n &&
        length(centers) == n && length(weights) == n ||
        throw(ArgumentError("supports, sigmas, centers, and weights must all have length $(n)"))

    # All wavepackets must share the same background
    left_gs  = states[1].psi.left_gs
    right_gs = states[1].psi.right_gs
    for i in 2:n
        states[i].psi.left_gs === left_gs && states[i].psi.right_gs === right_gs ||
            throw(ArgumentError("all states must share the same background MPS"))
    end

    Uc = length(left_gs)
    W % Uc == 0 || throw(ArgumentError("Window size W must be a multiple of the unit cell length $(Uc)"))

    # Validate supports: must start on odd sites and end on even sites (unit cell alignment)
    for (i, sup) in enumerate(supports)
        first(sup) >= 1 && last(sup) <= W ||
            throw(ArgumentError("supports[$i] must be contained in 1:$W"))
        isodd(first(sup)) ||
            throw(ArgumentError("supports[$i] must begin on an odd site (got $(first(sup)))"))
        iseven(last(sup)) ||
            throw(ArgumentError("supports[$i] must end on an even site (got $(last(sup)))"))
    end

    # Sort supports by starting site, reordering the associated arrays consistently
    perm     = sortperm(collect(supports), by=first)
    supports = collect(supports)[perm]
    states   = collect(states)[perm]
    sigmas   = collect(sigmas)[perm]
    centers  = collect(centers)[perm]
    weights  = collect(weights)[perm]

    # Validate that sorted supports are disjoint with at least 2 sites between them
    for i in 1:n-1
        last(supports[i]) + 2 < first(supports[i+1]) ||
            throw(ArgumentError("supports[$i] and supports[$(i+1)] must have at least 2 sites between them (not yet implemented overlap handling)"))
    end

    T = storagetype(MPSKit.site_type(left_gs))

    # Build weight dict k => w for QP ℓ
    function make_ws(ℓ)
        qp  = states[ℓ].psi
        p   = qp.momentum
        sup = supports[ℓ]
        support_ws = if !isnothing(weights[ℓ])
            collect(ComplexF64, weights[ℓ])
        else
            x0 = isnothing(centers[ℓ]) ? (first(sup) + last(sup)) / 2.0 : Float64(centers[ℓ])
            σ  = isnothing(sigmas[ℓ])  ? length(sup) / 4.0               : Float64(sigmas[ℓ])
            [exp(im * p * k) * exp(-(k - x0)^2 / (2σ^2)) for k in sup]
        end
        length(support_ws) == length(sup) ||
            throw(ArgumentError("length(weights[$ℓ]) must equal length(supports[$ℓ])=$(length(sup))"))
        return Dict(k => support_ws[idx] for (idx, k) in enumerate(sup))
    end

    # Fuse the auxiliary (excitation) leg of a raw B tensor into the left virtual leg
    function fuse_B(t, utl)
        frontmap = isomorphism(T, fuse(utl * MPSKit._firstspace(t)),
                               utl * MPSKit._firstspace(t))
        @plansor tt[-1 -2; -3] := t[1 -2; 2 -3] * frontmap[-1; 2 1]
        return tt
    end

    gauge in (:symmetric, :left) ||
        throw(ArgumentError("gauge must be :symmetric or :left (got $gauge)"))

    # Precompute fused, weighted B tensors for each QP at each site in its support, in the
    # chosen Bloch gauge. With `:left` every QP uses MPSKit's left gauge. With `:symmetric`
    # every QP uses the reflection-symmetric gauge of arXiv:2012.07243: each QP is made to
    # straddle its own AL→AR transition (see the mixed-canonical background below), so the
    # symmetric gauge centres it. Every particle in the window is treated identically.
    all_Bs = map(1:n) do ℓ
        qp  = states[ℓ].psi
        utl = MPSKit.auxiliaryspace(qp)
        ws  = make_ws(ℓ)
        Bs  = gauge === :symmetric ? reflection_symmetric_gauge(qp) : [qp[i] for i in 1:Uc]
        Dict(k => fuse_B(ws[k] * Bs[mod1(k, Uc)], utl) for k in supports[ℓ])
    end

    # AL/AR tensors for every window site
    ALs = [left_gs.AL[mod1(k, Uc)]  for k in 1:W]
    ARs = [right_gs.AR[mod1(k, Uc)] for k in 1:W]

    # Map each window site to its QP support index (0 = not in any support)
    site_to_support = zeros(Int, W)
    for ℓ in 1:n
        for k in supports[ℓ]
            site_to_support[k] = ℓ
        end
    end
    first_support_start = first(supports[1])     # supports are sorted by start
    last_support_end    = last(supports[n])

    # Mixed-canonical background. Every QP straddles AL→AR, so each leaves the vacuum in
    # right-canonical (AR) gauge. The window must reach the next QP (and any QP except the
    # last needs AL on its left) in left-canonical (AL) gauge, so each inter-QP gap carries
    # one AR→AL domain wall: the ground-state center matrix inverse C⁻¹ at a mid-gap bond
    # (ALₐ…ALᵦ ALᵦ₊₁… = Cₐ₋₁ ARₐ…ARᵦ Cᵦ⁻¹ ALᵦ₊₁…), absorbed into the AL just past the bond.
    # Left of the first QP is AL (left environment); right of the last QP is AR (right env).
    wall_bond = Dict{Int,Int}()                  # gap a (between support a, a+1) => bond b
    for a in 1:n-1
        ja, ia1 = last(supports[a]), first(supports[a+1])
        wall_bond[a] = ja + (ia1 - ja - 1) ÷ 2   # mid-gap bond (C⁻¹ between site b and b+1)
    end

    function vacuum_tensor(k)
        k < first_support_start && return ALs[k]
        k > last_support_end    && return ARs[k]
        a = findfirst(g -> last(supports[g]) < k < first(supports[g+1]), 1:n-1)
        b = wall_bond[a]
        k <= b && return ARs[k]                  # AR tail trailing QP a
        if k == b + 1                            # AR→AL domain wall, fused into AL[b+1]
            Cinv = _trunc_pinv(left_gs.C[mod1(b, Uc)]; rtol = wall_rtol)
            @plansor t[-1 -2; -3] := Cinv[-1; 1] * ALs[k][1 -2; -3]
            return t
        end
        return ALs[k]                            # AL approach to QP a+1
    end

    # Build the automaton block tensor at site k within support ℓ.
    #
    # The automaton has two local states per support:
    #   state 0 (iso1): "QP not yet inserted" → propagate with AL
    #   state 1 (iso2): "QP inserted"         → propagate with AR
    # so every QP straddles an AL→AR transition (centred by the symmetric gauge). The
    # vacuum returns AR→AL between consecutive QPs via a domain wall (see `vacuum_tensor`).
    #
    # Tensor shapes (VL = left_virtualspace, VR = right_virtualspace of AL[k]):
    #   left boundary  (k == first(sup)):  VL        ← Vphys ⊗ (VR ⊕ VR)
    #   interior sites (i < k < j):        (VL ⊕ VL) ← Vphys ⊗ (VR ⊕ VR)
    #   right boundary (k == last(sup)):   (VL ⊕ VL) ← Vphys ⊗ VR
    #
    # Isometries iso1, iso2 (both VV ← V) select the two halves of the doubled bond.
    # Following the MPSKit plansor embedding convention:
    #   left  embed: @plansor M[-1 -2; -3] += iso[-1; 1] * T[1 -2; -3]
    #   right embed: @plansor M[-1 -2; -3] += T[-1 -2; 1] * conj(iso[-3; 1])
    function build_automaton_tensor(k, ℓ)
        AL      = ALs[k]
        B       = all_Bs[ℓ][k]
        A_after = ARs[k]   # every QP straddles AL→AR; AR→AL reset happens in the gap

        sup          = supports[ℓ]
        is_left_bdy  = (k == first(sup))
        is_right_bdy = (k == last(sup))

        VL    = MPSKit.left_virtualspace(AL)
        VR    = MPSKit.right_virtualspace(AL)
        VL_dbl = TensorKit.:⊕(VL, VL)
        VR_dbl = TensorKit.:⊕(VR, VR)

        if is_left_bdy
            # Left boundary: 1D left bond → 2D right bond
            # Row vector [AL | B] — state 0 goes to AL, state 1 goes to B
            iso1_R = isometry(T, VR_dbl, VR)
            iso2_R = TensorKit.left_null(iso1_R)
            @plansor M[-1 -2; -3] := AL[-1 -2; 1] * conj(iso1_R[-3; 1])
            @plansor M[-1 -2; -3] += B[-1 -2; 1] * conj(iso2_R[-3; 1])

        elseif is_right_bdy
            # Right boundary: 2D left bond → 1D right bond
            # Column vector [B; A_after] — from state 0 via B, from state 1 via A_after
            iso1_L = isometry(T, VL_dbl, VL)
            iso2_L = TensorKit.left_null(iso1_L)
            @plansor M[-1 -2; -3] := iso1_L[-1; 1] * B[1 -2; -3]
            @plansor M[-1 -2; -3] += iso2_L[-1; 1] * A_after[1 -2; -3]

        else
            # Interior: full 2×2 block
            # [AL       B    ]
            # [0        A_after]
            iso1_L = isometry(T, VL_dbl, VL)
            iso2_L = TensorKit.left_null(iso1_L)
            iso1_R = isometry(T, VR_dbl, VR)
            iso2_R = TensorKit.left_null(iso1_R)
            @plansor M[-1 -2; -3] := iso1_L[-1; 1] * AL[1 -2; 2] * conj(iso1_R[-3; 2])
            @plansor M[-1 -2; -3] += iso1_L[-1; 1] * B[1 -2; 2] * conj(iso2_R[-3; 2])
            @plansor M[-1 -2; -3] += iso2_L[-1; 1] * A_after[1 -2; 2] * conj(iso2_R[-3; 2])
        end

        return M
    end

    # Assemble all window tensors
    tensors = map(1:W) do k
        ℓ = site_to_support[k]
        ℓ == 0 ? vacuum_tensor(k) : build_automaton_tensor(k, ℓ)
    end

    winfms = FiniteMPS(tensors; normalize=false)
    psi = window ? WindowMPS(left_gs, winfms, right_gs) : winfms
    return MPSKitState(states[1].hamiltonian, psi)
end
