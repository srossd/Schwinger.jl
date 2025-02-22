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
end