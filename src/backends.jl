"""
Backend types for Schwinger.jl

This module defines the backend abstraction layer that allows users to select
between different computational backends (ED, ITensors, MPSKit) using a unified API.
"""

"""
    Backend

Abstract type representing a computational backend for Schwinger model calculations.
"""
abstract type Backend end

"""
    EDBackend <: Backend

Exact diagonalization backend using sparse matrices.
"""
struct EDBackend <: Backend end

"""
    ITensorsBackend <: Backend

ITensors.jl backend using MPO/MPS representations.
"""
struct ITensorsBackend <: Backend end

"""
    MPSKitBackend <: Backend

MPSKit.jl backend using MPO/MPS representations.
"""
struct MPSKitBackend <: Backend end

# Singleton instances for convenient usage
const ED = EDBackend()
const ITensorsMPS = ITensorsBackend()
const MPSKitMPS = MPSKitBackend()

"""
    resolve_backend(b::Symbol) -> Backend
    resolve_backend(b::Backend) -> Backend

Convert a backend specifier (Symbol or Backend instance) to a Backend type.

# Supported symbols
- `:ED` - Exact diagonalization
- `:ITensors` - ITensors.jl
- `:MPSKit` - MPSKit.jl
"""
function resolve_backend(b::Symbol)
    b == :ED && return ED
    b == :ITensors && return ITensorsMPS
    b == :MPSKit && return MPSKitMPS
    throw(ArgumentError("Unknown backend: $b. Use :ED, :ITensors, or :MPSKit"))
end

resolve_backend(b::Backend) = b

# =============================================================================
# Default Backend Management
# =============================================================================

# Module-level variable to store the default backend
const _DEFAULT_BACKEND = Ref{Backend}(ITensorsBackend())

"""
    set_default_backend(backend::Union{Symbol,Backend})

Set the default backend for Schwinger.jl calculations.

# Arguments
- `backend`: Backend to use as default. Can be a Symbol (`:ED`, `:ITensors`, `:MPSKit`) or a Backend instance.

# Examples
```julia
set_default_backend(MPSKitBackend())
set_default_backend(:MPSKit)
```
"""
function set_default_backend(backend::Union{Symbol,Backend})
    _DEFAULT_BACKEND[] = resolve_backend(backend)
    return nothing
end

"""
    get_default_backend() -> Backend

Get the current default backend for Schwinger.jl calculations.
"""
get_default_backend() = _DEFAULT_BACKEND[]
