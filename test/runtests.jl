using Schwinger
using Test

@testset "Time-independent tests" begin
    @testset "ED vs MPO" begin
        include("ed_mpo.jl")
    end

    @testset "Ground state" begin
        include("groundstate.jl")
    end

    @testset "Mass gap" begin
        include("massgap.jl")
    end

    @testset "Defect charges" begin
        include("defects.jl")
    end

    @testset "Fermion field" begin
        include("fermionfield.jl")
    end

    @testset "Wilson lines" begin
        include("wilson_line.jl")
    end

    @testset "Energy densities" begin
        include("energy_density.jl")
    end

    @testset "Currents" begin
        include("currents.jl")
    end

    @testset "Pseudoscalar density conventions" begin
        include("pseudoscalar_convention.jl")
    end

    @testset "Flavor symmetry" begin
        include("flavor_symmetry.jl")
    end
end

@testset "Time evolution tests" begin
    @testset "ED vs MPO" begin
        include("ed_mpo_time.jl")
    end

    @testset "Wilson loop correlator" begin
        include("wilson_correlator.jl")
    end

    @testset "Evolve checkpoint hook" begin
        include("evolve_checkpoint.jl")
    end

    @testset "Energy density convention" begin
        include("energy_convention.jl")
    end
end

@testset "Infinite lattice tests" begin
    @testset "Ground state" begin
        include("infinite_groundstate.jl")
    end

    @testset "Reflection-symmetric QP gauge" begin
        include("wavepacket_gauge.jl")
    end

    @testset "EMT cell operator" begin
        include("emt.jl")
    end

    @testset "Adaptive window growth" begin
        include("grow_window.jl")
    end

    @testset "Momentum density" begin
        include("momentum_density.jl")
    end
end

@testset "Cached ground state reuse" begin
    include("cached_groundstate.jl")
end

@testset "Quick-win helpers" begin
    include("quickwins.jl")
end

@testset "loweststates accessors" begin
    include("loweststates_accessors.jl")
end

@testset "Detector pad option" begin
    include("detector_pad.jl")
end

@testset "apply_local on window" begin
    include("apply_local.jl")
end

@testset "Two-point correlator and chargeprofile" begin
    include("correlator_chargeprofile.jl")
end

@testset "MPSKit energycurrents (local path)" begin
    include("energy_local.jl")
end

@testset "Window workflow: gs -> WindowMPS -> act -> evolve+grow" begin
    include("window_workflow.jl")
end