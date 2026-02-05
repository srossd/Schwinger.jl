using Base.MathConstants
using Statistics
using Integrals
using ProgressMeter

@testset "Average Electric Field" begin
    @showprogress desc = "Average Electric Field: m = 0" for θ2π=0:.05:1
        lattice = Lattice(10; periodic=true, θ2π=θ2π)
        ground = groundstate(Hamiltonian(lattice, EDBackend()))
        
        efs = electricfields(ground)
        avgE = mean(efs)

        @test avgE ≈ 0 atol=0.01
    end

    m = 0.1
    @showprogress desc = "Average Electric Field: m/g = 0.1" for θ2π=0:.05:1
        lattice=Lattice(10; periodic=true, θ2π=θ2π, m = m)
        ground = groundstate(Hamiltonian(lattice, EDBackend()))
        
        efs = electricfields(ground)
        avgE = mean(efs)

        # See eq (24) of https://arxiv.org/abs/2206.05308
        exact = (exp(γ)/√(π))*m*sin(2π*θ2π) - (8.9139*exp(2γ)/(4π))*(m^2)*sin(4π*θ2π)

        @test avgE ≈ exact atol=0.01
    end
end

@testset "Scalar VEV" begin
    f(x, L) = 1/(1 - exp(L*cosh(x)/√(π)))

    @showprogress desc = "Chiral Condensate (m = 0): L dependence" for L=1:.5:10
        lattice = Lattice(10; periodic=true, L=L)
        ground = groundstate(EDHamiltonian(lattice))
        
        condensate = scalar(ground)

        # See eq (7) of https://arxiv.org/abs/2206.05308
        exact = -exp(γ)/(2π^(3/2))*exp(2*solve(IntegralProblem(f, (0, Inf), L), QuadGKJL()).u)

        @test condensate ≈ exact rtol=0.03
    end

    L=8
    lattice = Lattice(10; periodic=true, L=L)
    ground = groundstate(EDHamiltonian(lattice))
    condensate = scalar(ground)
    @showprogress desc = "Chiral Condensate (m = 0): θ dependence" for θ2π=0:.04:1
        lattice = Lattice(10; periodic=true, L=L, θ2π=θ2π)
        ground = groundstate(EDHamiltonian(lattice))
        
        condensate_ratio = scalar(ground)/condensate

        exact = cos(2π*θ2π)

        @test condensate_ratio ≈ exact rtol=0.01
    end
end

@testset "Wilson loop VEV" begin
    @showprogress desc = "Wilson loop VEV" for L=1:10:1
        lattice = Lattice(10; F = 1, periodic=true, L=L)
        ground = groundstate(EDHamiltonian(lattice))
        
        vev = expectation(EDWilsonLoop(lattice), ground)

        exact = exp(-√(π)/4*L)

        @test vev ≈ exact rtol=0.05
    end
end

@testset "Higher charge" begin
    @showprogress desc = "Universes" for p=0:3
        lattice = Lattice(10; F = 1, periodic=true, q = 4)
        gs = groundstate(EDHamiltonian(lattice; universe = p))

        arg = atan(pseudoscalar(gs), scalar(gs))

        @test mod2pi(arg + 1) ≈ mod2pi(π*(1+p/2) + 1) atol=0.01
    end
end