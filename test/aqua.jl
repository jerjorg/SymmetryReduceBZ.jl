using Aqua, Test, SymmetryReduceBZ

@testset "method ambiguity" Aqua.test_ambiguities(SymmetryReduceBZ)
Aqua.test_all(SymmetryReduceBZ, ambiguities=false)
