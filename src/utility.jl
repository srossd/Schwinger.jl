function ITensors.op(::OpName"L0", ::SiteType"Boson", d::Int)
    mat = zeros(d, d)
    for k in 1:d
      mat[k, k] = (k - 1) - (d ÷ 2)
    end
    return mat
end
function ITensors.op(::OpName"raise", ::SiteType"Boson", d::Int)
    mat = zeros(d, d)
    for k in 1:d-1
      mat[k + 1, k] = 1
    end
    return mat
end
function ITensors.op(::OpName"lower", ::SiteType"Boson", d::Int)
    mat = zeros(d, d)
    for k in 1:d-1
      mat[k, k + 1] = 1
    end
    return mat
end

tuplejoin(x) = x
tuplejoin(x, y) = (x..., y...)
tuplejoin(x, y, z...) = tuplejoin(tuplejoin(x, y), z...)

function wigner_string(N::Int, F::Int, flavor::Int)
    return tuplejoin((("Sz",(j-1)*F + flavor) for j in 1:N)...)
end

@memoize function sites(lattice::SchwingerLattice{N,F}; L_max::Int = 3) where {N,F}
    total_sites = lattice.periodic ? N * F + 1 : N * F
    sites = Vector{Index}(undef, total_sites)

    # Create each site in the SchwingerLattice
    for i in 1:total_sites
        if lattice.periodic && i == total_sites
            sites[i] = Index([QN() => 2L_max + 1]; tags="Boson,Site,n=$i")
        else
            sites[i] = siteind("S=1/2", i; conserve_qns = true)
        end
    end

    return sites
end