module SymmetryReduceBZ

include("Lattices.jl")
include("Utilities.jl")
include("Symmetry.jl")
include("Plotting.jl")

import .Plotting: plot_convexhulls
import .Symmetry: calc_bz, calc_ibz

export calc_bz, calc_ibz, plot_convexhulls

if !isdefined(Base, :get_extension)
    include("../ext/SymmetryReduceBZPyPlotExt.jl")
    include("../ext/SymmetryReduceBZQHullExt.jl")
end

end #module
