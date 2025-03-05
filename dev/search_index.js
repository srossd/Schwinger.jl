var documenterSearchIndex = {"docs":
[{"location":"index.html#Installation","page":"Index","title":"Installation","text":"","category":"section"},{"location":"index.html","page":"Index","title":"Index","text":"To install Schwinger.jl, use the following command in the Julia REPL:","category":"page"},{"location":"index.html","page":"Index","title":"Index","text":"using Pkg\nPkg.add(\"https://github.com/srossd/Schwinger.jl\")","category":"page"},{"location":"index.html#Key-Features","page":"Index","title":"Key Features","text":"","category":"section"},{"location":"index.html","page":"Index","title":"Index","text":"Schwinger.jl provides functions for computing the following properties of the Hamiltonian lattice Schwinger model:","category":"page"},{"location":"index.html","page":"Index","title":"Index","text":"Energy levels: Ground state, energy gap, or higher excited states.\nCorrelators: Expectation values of various operators in the ground state (or other states).\nTime evolution: Action of time evolution on a given initial state.","category":"page"},{"location":"index.html","page":"Index","title":"Index","text":"All of these features are implemented using both exact diagonalization (ED) and matrix product operators/states (MPO).","category":"page"},{"location":"index.html#Usage-Example","page":"Index","title":"Usage Example","text":"","category":"section"},{"location":"index.html","page":"Index","title":"Index","text":"Here is a basic example of how to use Schwinger.jl to calculate the average electric field in the one-flavor Schwinger model at mg = 01 as a function of theta. Note that the mass shift (see here) is included by default.","category":"page"},{"location":"index.html","page":"Index","title":"Index","text":"using Schwinger\nusing Base.MathConstants\nusing Plots\n\nfunction avgE(θ2π, m, shift = true)\n    lat = shift ? SchwingerLattice{10,1}(θ2π = θ2π, m = m, periodic = true) : SchwingerLattice{10,1}(θ2π = θ2π, mlat = m, periodic = true)\n    gs = groundstate(EDHamiltonian(lat))\n    return real(expectation(EDAverageElectricField(lat), gs))\nend\n\nm = 0.05\nθ2πs = 0:0.025:1\navgEs_shift = map(x -> avgE(x, m), θ2πs)\navgEs_noshift = map(x -> avgE(x, m, false), θ2πs)\n\n# See eq (24) of https://arxiv.org/abs/2206.05308\nperturbative = [(exp(γ)/√(π))*m*sin(2π*θ2π) - (8.9139*exp(2γ)/(4π))*(m^2)*sin(4π*θ2π) for θ2π in θ2πs]\n\nscatter(θ2πs, avgEs_shift, label=\"With mass shift\", xlabel=\"θ/2π\", ylabel=\"Average Electric Field\", legend=:topright, color=:orange)\nscatter!(θ2πs, avgEs_noshift, label=\"Without mass shift\", xlabel=\"θ/2π\", ylabel=\"Average Electric Field\", color=\"lightblue\")\nplot!(θ2πs, perturbative, label=\"Perturbative\", color=:black)","category":"page"},{"location":"index.html","page":"Index","title":"Index","text":"Here is an example of calculating the expectation value of the square of the mean electric field in the two-flavor Schwinger model at theta = pi, for a four-site periodic lattice, giving a very rough look at the phase diagram.","category":"page"},{"location":"index.html","page":"Index","title":"Index","text":"using Schwinger\nusing Plots\n\nfunction avgE2(m1, m2)\n    lat = SchwingerLattice{4,2}(θ2π = 0.5, m = (m1, m2), periodic = true)\n    gs = groundstate(EDHamiltonian(lat))\n    return real(expectation(EDAverageElectricField(lat; power=2), gs))\nend\n\nms = 0:0.05:1.5\nm1s = repeat(reshape(ms, 1, :), length(ms), 1)\nm2s = repeat(ms, 1, length(ms))\nvals = map(avgE2, m1s, m2s)\n\ncontour(ms, ms, vals, fill = true, lw = 0, xlabel = \"m₁/g\", ylabel = \"m₂/g\", clabel = \"⟨E²⟩\")","category":"page"}]
}
