using Base.MathConstants
using Statistics
using Integrals
using ProgressMeter

@testset "Average Electric Field" begin
    lattice = SchwingerLattice(10, 1, L_max = 3, θ2π = 0, a = 1, m = 0)
    @showprogress desc = "Average Electric Field: m = 0" for θ2π=0:.05:1
        setθ2π!(lattice, θ2π)
        findstates!(lattice, 1; outputlevel = 0)
        
        efs = electricfields(lattice)[1]
        avgE = mean(efs)

        @test avgE ≈ 0 atol=0.01
    end

    m = 0.1
    setmass!(lattice, m)
    @showprogress desc = "Average Electric Field: m/g = 0.1" for i=0:.05:1
        setθ2π!(lattice, θ2π)
        findstates!(lattice, 1; outputlevel = 0)
        
        efs = electricfields(lattice)[1]
        avgE = mean(efs)

        # See eq (25) of https://arxiv.org/abs/2206.05308
        exact = (exp(γ)/√(π))*m*sin(2*π*θ2π) - (8.9139*exp(2*γ)/(4*π))*(m^2)*sin(4*π*θ2π)

        @test avgE ≈ exact atol=0.01
    end
end

@testset "Chiral Condensate" begin
    lattice = SchwingerLattice(10, 1, L_max = 3, θ2π = 0, a = 1, m = 0)
    f(x, L) = 1/(1 - exp(L*cosh(x)/√(π)))

    @showprogress desc = "Chiral Condensate: m = 0" for i=0:20
        θ = 2*π*i/20
        setθ2π!(lattice, θ/(2*π))
        findstates!(lattice, 1; outputlevel = 0)
        
        efs = electricfields(lattice)[1]
        avgE = mean(efs)

        # See eq (25) of https://arxiv.org/abs/2206.05308
        exact = (exp(γ)/√(π))*m*sin(θ) - (8.9139*exp(2*γ)/(4*π))*(m^2)*sin(2*θ)

        @test avgE ≈ exact rtol=0.03
    end
end