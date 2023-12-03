using Aqua
using SymmetryReduceBZ

Aqua.test_all(SymmetryReduceBZ, ambiguities=(broken=true,))

include("lattices.jl")
include("symmetry.jl")
include("utilities.jl")
include("plotting.jl")
