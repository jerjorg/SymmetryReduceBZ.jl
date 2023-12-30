using Aqua, Test
using SymmetryReduceBZ

@testset "method ambiguity" Aqua.test_ambiguities(SymmetryReduceBZ)
Aqua.test_all(SymmetryReduceBZ, ambiguities=false)

@testset    "lattices"  include("lattices.jl")
@testset    "symmetry"  include("symmetry.jl")
@testset    "utilities" include("utilities.jl")
@testset    "plotting"  include("plotting.jl")
