using Schwinger
using Test

@testset "Time-independent tests" begin
    @testset "Ground state tests" begin
        include("ground_state_tests.jl")
    end

    # @testset "Energy gap tests" begin
    #     include("energy_gap_tests.jl")
    # end
end