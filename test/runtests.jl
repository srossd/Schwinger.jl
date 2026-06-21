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

    @testset "Zero-length Wilson line" begin
        include("wilson_line.jl")
    end

    @testset "Energy densities" begin
        include("energy_density.jl")
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
end

@testset "Infinite lattice tests" begin
    @testset "Ground state" begin
        include("infinite_groundstate.jl")
    end

    @testset "Reflection-symmetric QP gauge" begin
        include("wavepacket_gauge.jl")
    end
end