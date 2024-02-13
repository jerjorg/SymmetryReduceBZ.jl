using Test

@testset    "aqua"      begin include("aqua.jl") end
@testset    "lattices"  begin include("lattices.jl") end
@testset    "symmetry"  begin include("symmetry.jl") end
@testset    "utilities" begin include("utilities.jl") end
@testset    "plotting"  begin include("plotting.jl") end
