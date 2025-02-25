using Base.MathConstants
using Statistics
using Integrals
using ProgressMeter

@testset "Average Electric Field" begin
    @showprogress desc = "Average Electric Field: m = 0" for θ2π=0:.05:1
        lattice = SchwingerLattice{10,1}(periodic=true, θ2π=θ2π)
        ground = groundstate(EDHamiltonian(lattice))
        
        efs = electricfields(ground)
        avgE = mean(efs)

        @test avgE ≈ 0 atol=0.01
    end

    m = 0.1
    @showprogress desc = "Average Electric Field: m/g = 0.1" for θ2π=0:.05:1
        lattice=SchwingerLattice{10,1}(periodic=true, θ2π=θ2π, m = m)
        ground = groundstate(EDHamiltonian(lattice))
        
        efs = electricfields(ground)
        avgE = mean(efs)

        # See eq (24) of https://arxiv.org/abs/2206.05308
        exact = (exp(γ)/√(π))*m*sin(2π*θ2π) - (8.9139*exp(2γ)/(4π))*(m^2)*sin(4π*θ2π)

        @test avgE ≈ exact atol=0.01
    end
end

@testset "Scalar VEV" begin
    f(x, L) = 1/(1 - exp(L*cosh(x)/√(π)))

    N=10
    @showprogress desc = "Chiral Condensate (m = 0): L dependence" for L=1:.5:10
        lattice = SchwingerLattice{N,1}(periodic=true, L=L)
        ground = groundstate(EDHamiltonian(lattice))
        
        condensate = scalarvev(ground)

        # See eq (7) of https://arxiv.org/abs/2206.05308
        exact = -exp(γ)/(2π^(3/2))*exp(2*solve(IntegralProblem(f, (0, Inf), L), QuadGKJL()).u)

        @test condensate ≈ exact rtol=0.03
    end

    L=8
    lattice = SchwingerLattice{N,1}(periodic=true, L=L)
    ground = groundstate(EDHamiltonian(lattice))
    condensate = scalarvev(ground)
    @showprogress desc = "Chiral Condensate (m = 0): θ dependence" for θ2π=0:.04:1
        lattice = SchwingerLattice{N,1}(periodic=true, L=L, θ2π=θ2π)
        ground = groundstate(EDHamiltonian(lattice))
        
        condensate_ratio = scalarvev(ground)/condensate

        exact = cos(2π*θ2π)

        @test condensate_ratio ≈ exact rtol=0.01
    end
end

@testset "Wilson loop VEV" begin
    @showprogress desc = "Wilson loop VEV" for L=1:10:1
        lattice = SchwingerLattice{10,1}(periodic=true, L=L)
        ground = groundstate(EDHamiltonian(lattice))
        
        vev = expectation(EDWilsonLoop(lattice), ground)

        exact = exp(-√(π)/4*L)

        @test vev ≈ exact rtol=0.05
    end
end