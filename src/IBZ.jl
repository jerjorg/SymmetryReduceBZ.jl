module IBZ

include("Lattices.jl")
include("Plotting.jl")
include("Symmetry.jl")
include("Utilities.jl")

import .Plotting: plot_convexhulls
import .Symmetry: calc_bz, calc_ibz

export calc_bz, calc_ibz, plot_convexhulls

# using .Lattices, .Plotting, .Symmetry, .Utilities
#
# # lattices
# export get_recip_latvecs, minkowski_reduce, check_reduced
# # 2D Bravais lattices
# export genlat_SQR, genlat_HXG, genlat_REC, genlat_RECI, genlat_OBL
# # 3D Bravais lattices
# export genlat_CUB, genlat_FCC, genlat_BCC, genlat_TET, genlat_BCT, genlat_ORC,
#     genlat_ORCF, genlat_ORCI, genlat_ORCC, genlat_HEX, genlat_RHL, genlat_MCL,
#     genlat_MCLC, genlat_TRI
#
# # plotting
# export plot_convexhulls
#
# # symmetry
# export mapto_unitcell, calc_pointgroup, calc_spacegroup, calc_bz, calc_ibz

end #module
