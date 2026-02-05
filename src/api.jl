"""
Unified API for Schwinger.jl

This module provides backend-agnostic functions that dispatch to the appropriate
backend implementation based on the `backend` keyword argument.
"""

# =============================================================================
# Hamiltonian
# =============================================================================

"""
    Hamiltonian(lattice; backend=nothing, L_max=nothing, universe=0)

Construct a Schwinger model Hamiltonian using the specified backend.

# Arguments
- `lattice::Lattice`: The lattice configuration.
- `backend`: Backend to use. Can be a Symbol (`:ED`, `:ITensors`, `:MPSKit`), a Backend instance, or `nothing` to use the default backend.
- `L_max::Union{Nothing,Int}`: Maximum angular momentum cutoff.
- `universe::Int`: Charge sector.

# Returns
- `EDOperator` for `:ED` backend
- `ITensorOperator` for `:ITensors` backend
- `MPSKitOperator` for `:MPSKit` backend

# Examples
```julia
# Use default backend (ITensors unless changed via set_default_backend)
H = Hamiltonian(lattice)

# Specify backend explicitly
H = Hamiltonian(lattice; backend=:ED)
H = Hamiltonian(lattice; backend=:MPSKit)

# Change default backend
set_default_backend(:MPSKit)
H = Hamiltonian(lattice)  # Now uses MPSKit
```
"""
function Hamiltonian(lattice::Lattice; backend::Union{Symbol,Backend,Nothing}=nothing, kwargs...)
    backend = isnothing(backend) ? get_default_backend() : resolve_backend(backend)
    return Hamiltonian(lattice, backend; kwargs...)
end

Hamiltonian(lattice::Lattice, ::EDBackend; kwargs...) = EDHamiltonian(lattice; kwargs...)
Hamiltonian(lattice::Lattice, ::ITensorsBackend; kwargs...) = ITensorHamiltonian(lattice; kwargs...)
Hamiltonian(lattice::Lattice, ::MPSKitBackend; kwargs...) = MPSKitHamiltonian(lattice; kwargs...)

# =============================================================================
# GaugeKinetic
# =============================================================================

"""
    GaugeKinetic(lattice; backend=nothing, L_max=nothing, universe=0, bare=true)

Construct the gauge kinetic operator ∑(Lₙ+θ/2π)² using the specified backend.

# Arguments
- `backend`: Backend to use (`:ED`, `:ITensors`, `:MPSKit`), or `nothing` for default.
"""
function GaugeKinetic(lattice::Lattice; backend::Union{Symbol,Backend,Nothing}=nothing, kwargs...)
    backend = isnothing(backend) ? get_default_backend() : resolve_backend(backend)
    return GaugeKinetic(lattice, backend; kwargs...)
end

GaugeKinetic(lattice::Lattice, ::EDBackend; kwargs...) = EDGaugeKinetic(lattice; kwargs...)
GaugeKinetic(lattice::Lattice, ::ITensorsBackend; kwargs...) = ITensorGaugeKinetic(lattice; kwargs...)
GaugeKinetic(lattice::Lattice, ::MPSKitBackend; kwargs...) = MPSKitGaugeKinetic(lattice; kwargs...)

# =============================================================================
# Mass
# =============================================================================

"""
    Mass(lattice; backend=nothing, L_max=nothing, universe=0, bare=true)

Construct the mass operator ∑ (-1)ⁿ χ†ₙχₙ using the specified backend.

# Arguments
- `backend`: Backend to use (`:ED`, `:ITensors`, `:MPSKit`), or `nothing` for default.
"""
function Mass(lattice::Lattice; backend::Union{Symbol,Backend,Nothing}=nothing, kwargs...)
    backend = isnothing(backend) ? get_default_backend() : resolve_backend(backend)
    return Mass(lattice, backend; kwargs...)
end

Mass(lattice::Lattice, ::EDBackend; kwargs...) = EDMass(lattice; kwargs...)
Mass(lattice::Lattice, ::ITensorsBackend; kwargs...) = ITensorMass(lattice; kwargs...)
Mass(lattice::Lattice, ::MPSKitBackend; kwargs...) = MPSKitMass(lattice; kwargs...)

# =============================================================================
# Hopping
# =============================================================================

"""
    Hopping(lattice; backend=nothing, L_max=nothing, universe=0, bare=true)

Construct the hopping term -i ∑(χ†ₙ χₙ₊₁ - χ†ₙ₊₁ χₙ) using the specified backend.

# Arguments
- `backend`: Backend to use (`:ED`, `:ITensors`, `:MPSKit`), or `nothing` for default.
"""
function Hopping(lattice::Lattice; backend::Union{Symbol,Backend,Nothing}=nothing, kwargs...)
    backend = isnothing(backend) ? get_default_backend() : resolve_backend(backend)
    return Hopping(lattice, backend; kwargs...)
end

Hopping(lattice::Lattice, ::EDBackend; kwargs...) = EDHopping(lattice; kwargs...)
Hopping(lattice::Lattice, ::ITensorsBackend; kwargs...) = ITensorHopping(lattice; kwargs...)
Hopping(lattice::Lattice, ::MPSKitBackend; kwargs...) = MPSKitHopping(lattice; kwargs...)

# =============================================================================
# HoppingMass
# =============================================================================

"""
    HoppingMass(lattice; backend=nothing, L_max=nothing, universe=0, bare=true)

Construct the hopping-type mass term i/2 ∑(-1)^n (χ†ₙ₊₁ χₙ + χ†ₙ₋₁ χₙ) using the specified backend.

# Arguments
- `backend`: Backend to use (`:ED`, `:ITensors`, `:MPSKit`), or `nothing` for default.
"""
function HoppingMass(lattice::Lattice; backend::Union{Symbol,Backend,Nothing}=nothing, kwargs...)
    backend = isnothing(backend) ? get_default_backend() : resolve_backend(backend)
    return HoppingMass(lattice, backend; kwargs...)
end

HoppingMass(lattice::Lattice, ::EDBackend; kwargs...) = EDHoppingMass(lattice; kwargs...)
HoppingMass(lattice::Lattice, ::ITensorsBackend; kwargs...) = ITensorHoppingMass(lattice; kwargs...)
HoppingMass(lattice::Lattice, ::MPSKitBackend; kwargs...) = MPSKitHoppingMass(lattice; kwargs...)

# Site-specific version
function HoppingMass(lattice::Lattice, site::Int; backend::Union{Symbol,Backend,Nothing}=nothing, kwargs...)
    backend = isnothing(backend) ? get_default_backend() : resolve_backend(backend)
    return HoppingMass(lattice, site, backend; kwargs...)
end

HoppingMass(lattice::Lattice, site::Int, ::EDBackend; kwargs...) = EDHoppingMass(lattice, site; kwargs...)
HoppingMass(lattice::Lattice, site::Int, ::ITensorsBackend; kwargs...) = ITensorHoppingMass(lattice, site; kwargs...)
HoppingMass(lattice::Lattice, site::Int, ::MPSKitBackend; kwargs...) = MPSKitHoppingMass(lattice, site; kwargs...)

# =============================================================================
# WilsonLoop
# =============================================================================

"""
    WilsonLoop(lattice, conjugate=false; backend=nothing, L_max=nothing, universe=0)

Construct the spatial Wilson loop operator using the specified backend.

# Arguments
- `backend`: Backend to use (`:ED`, `:ITensors`, `:MPSKit`), or `nothing` for default.
"""
function WilsonLoop(lattice::Lattice, conjugate::Bool=false; backend::Union{Symbol,Backend,Nothing}=nothing, kwargs...)
    backend = isnothing(backend) ? get_default_backend() : resolve_backend(backend)
    return WilsonLoop(lattice, conjugate, backend; kwargs...)
end

WilsonLoop(lattice::Lattice, conjugate::Bool, ::EDBackend; kwargs...) = EDWilsonLoop(lattice, conjugate; kwargs...)
WilsonLoop(lattice::Lattice, conjugate::Bool, ::ITensorsBackend; kwargs...) = ITensorWilsonLoop(lattice, conjugate; kwargs...)
WilsonLoop(lattice::Lattice, conjugate::Bool, ::MPSKitBackend; kwargs...) = MPSKitWilsonLoop(lattice, conjugate; kwargs...)

# =============================================================================
# WilsonLine
# =============================================================================

"""
    WilsonLine(lattice, conjugate=false, flavor=1, start=nothing, finish=nothing;
               backend=nothing, L_max=nothing, universe=0)

Construct the spatial Wilson line operator using the specified backend.

# Arguments
- `backend`: Backend to use (`:ED`, `:ITensors`, `:MPSKit`), or `nothing` for default.
"""
function WilsonLine(lattice::Lattice, conjugate::Bool=false, flavor::Int=1,
                    start::Int=1, finish::Int=Int(lattice.N);
                    backend::Union{Symbol,Backend,Nothing}=nothing, kwargs...)
    backend = isnothing(backend) ? get_default_backend() : resolve_backend(backend)
    return WilsonLine(lattice, conjugate, flavor, start, finish, backend; kwargs...)
end

WilsonLine(lattice::Lattice, conjugate::Bool, flavor::Int, start, finish, ::EDBackend; kwargs...) =
    EDWilsonLine(lattice, conjugate, flavor, start, finish; kwargs...)
WilsonLine(lattice::Lattice, conjugate::Bool, flavor::Int, start, finish, ::ITensorsBackend; kwargs...) =
    ITensorWilsonLine(lattice, conjugate, flavor, start, finish; kwargs...)
WilsonLine(lattice::Lattice, conjugate::Bool, flavor::Int, start, finish, ::MPSKitBackend; kwargs...) =
    MPSKitWilsonLine(lattice, conjugate, flavor, start, finish; kwargs...)

# =============================================================================
# AverageElectricField
# =============================================================================

"""
    AverageElectricField(lattice; backend=nothing, power=1, L_max=nothing, universe=0, sitelist=nothing)

Construct an operator that computes the average electric field (raised to some power).

# Arguments
- `backend`: Backend to use (`:ED`, `:ITensors`, `:MPSKit`), or `nothing` for default.
"""
function AverageElectricField(lattice::Lattice; backend::Union{Symbol,Backend,Nothing}=nothing, kwargs...)
    backend = isnothing(backend) ? get_default_backend() : resolve_backend(backend)
    return AverageElectricField(lattice, backend; kwargs...)
end

AverageElectricField(lattice::Lattice, ::EDBackend; kwargs...) = EDAverageElectricField(lattice; kwargs...)
AverageElectricField(lattice::Lattice, ::ITensorsBackend; kwargs...) = ITensorAverageElectricField(lattice; kwargs...)
AverageElectricField(lattice::Lattice, ::MPSKitBackend; kwargs...) = MPSKitAverageElectricField(lattice; kwargs...)
