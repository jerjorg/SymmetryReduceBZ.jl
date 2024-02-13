using Aqua, Test, SymmetryReduceBZ

@testset "method ambiguity" begin Aqua.test_ambiguities(SymmetryReduceBZ) end
Aqua.test_all(SymmetryReduceBZ, ambiguities=false)
