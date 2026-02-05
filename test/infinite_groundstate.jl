using Test
using Schwinger
using ProgressMeter

@testset "Infinite vs Finite Lattice Groundstate Comparison" begin
    # Test parameters
    test_params = [
        (F = 1, a = 0.5, m = 1.0, θ2π = 0.33),
        (F = 1, a = 0.7, m = 0.5, θ2π = -0.3),
        (F = 1, a = 0.5, m = 1.0, θ2π = 0.1),
        (F = 1, a = 0.9, m = 1.5, θ2π = 0.25),
    ]
    
    # Tolerance for comparison
    tol = 1e-2
    bonddim = 50
    N_finite = 20
    middle_idxs = [9, 10]
    
    @showprogress for params in test_params
        # Create infinite and large finite lattices
        lat_inf = Lattice(Inf; params...)
        lat_finite = Lattice(N_finite; params...)
        
        # Compute groundstates
        mpskit_gs_inf = groundstate(Hamiltonian(lat_inf, MPSKitBackend()); bonddim=bonddim)
        mpskit_gs_finite = groundstate(Hamiltonian(lat_finite, MPSKitBackend()); bonddim=bonddim)
        
        # Occupations
        occ_inf = occupations(mpskit_gs_inf)
        occ_finite = occupations(mpskit_gs_finite)
        @test isapprox(occ_inf, occ_finite[middle_idxs]; rtol=tol)
        
        # Scalar densities
        scalar_inf = scalardensities(mpskit_gs_inf)
        scalar_finite = scalardensities(mpskit_gs_finite)
        @test isapprox(scalar_inf, scalar_finite[middle_idxs]; rtol=tol)
        
        # Pseudoscalar densities
        ps_inf = pseudoscalardensities(mpskit_gs_inf)
        ps_finite = pseudoscalardensities(mpskit_gs_finite)
        @test isapprox(ps_inf, ps_finite[middle_idxs]; rtol=tol)
        
        # Charges
        charges_inf = charges(mpskit_gs_inf)
        charges_finite = charges(mpskit_gs_finite)
        @test isapprox(charges_inf, charges_finite[middle_idxs]; rtol=tol)
        
        # Entanglements
        ent_inf = entanglements(mpskit_gs_inf)
        ent_finite = entanglements(mpskit_gs_finite)
        @test isapprox(ent_inf, ent_finite[middle_idxs]; rtol=tol)
    end
end