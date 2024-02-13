using Aqua, Test
using SymmetryReduceBZ

@testset "method ambiguity" Aqua.test_ambiguities(SymmetryReduceBZ)
Aqua.test_all(SymmetryReduceBZ, ambiguities=false)

@testset    "lattices"  begin include("lattices.jl") end
@testset    "symmetry"  begin include("symmetry.jl") end
@testset    "utilities" begin include("utilities.jl") end
@testset    "plotting"  begin include("plotting.jl") end
