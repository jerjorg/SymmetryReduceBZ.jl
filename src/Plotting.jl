module Plotting

export plot_convexhulls

@doc """
    plot_2Dconvexhull(convexhull, ax, color)

Plot a 2D convex hull

!!! note "Julia 1.9 and above"
    This function is implemented in a package extension and requires `using PyPlot`.

# Arguments
- `convexhull`: a polyhedron or convex hull object.
- `ax::PyObject`: an axes object from matplotlib.
- `facecolor::String="blue"`: the color of the area within the convex hull.
- `alpha::Real=0.3`: the transparency of the convex hull.
- `linewidth::Real=3`: the width of the edges.
- `edgecolor::String="black"`: the color of the edges.

# Returns
- `ax::PyObject`: updated `ax` that includes a plot of the convex hull.

# Examples
```julia
ENV["MPLBACKEND"]="qt5agg"
using PyPlot
using SymmetryReduceBZ.Symmetry: calc_bz, calc_ibz
using SymmetryReduceBZ.Plotting: plot_2Dconvexhull
real_latvecs = [1 0; 0 1]
convention="ordinary"
atom_types=[0]
atom_pos = Array([0 0]')
coords = "Cartesian"
makeprim=false
bz = calc_bz(real_latvecs,atom_types,atom_pos,coords,makeprim,convention)
ibz = calc_ibz(real_latvecs,atom_types,atom_pos,coords,makeprim,convention)
ax = plot_2Dconvexhull(bz,facecolor="deepskyblue",linewidth=3,edgecolor="cyan",alpha=0.2)
ax = plot_2Dconvexhull(ibz,ax;facecolor="coral",linewidth=3,edgecolor="magenta",alpha=0.4)
# output
PyObject <AxesSubplot: >
```
"""
function plot_2Dconvexhull end

@doc """
    plot_3Dconvexhull(convexhull,ax;color)

Plot a 3D convex hull

!!! note "Julia 1.9 and above"
    This function is implemented in a package extension and requires `using PyPlot`.

# Arguments
- `convexhull`: a polyhedron or convex hull object.
- `ax::PyObject`: an axes object from matplotlib.
- `facecolors::String="blue"`: the color of the faces of the convex hull.
- `alpha::Real=0.3`: the transparency of the faces of the convex hull.
- `linewidths::Real=1`: the width of the edges of the convex hull.
- `edgecolors::String="black"`: the color of the edges of the convex hull.

# Returns
- `ax::PyObject`: updated `ax` that includes a plot of the convex hull.

# Examples
```julia
ENV["MPLBACKEND"]="qt5agg"
using PyPlot
using SymmetryReduceBZ.Symmetry: calc_bz, calc_ibz
using SymmetryReduceBZ.Plotting: plot_3Dconvexhull
real_latvecs = [1 0 0; 0 1 0; 0 0 1]
convention="ordinary"
atom_types=[0]
atom_pos = Array([0 0 0]')
coords = "Cartesian"
makeprim=false
bz = calc_bz(real_latvecs,atom_types,atom_pos,coords,makeprim,convention)
ibz = calc_ibz(real_latvecs,atom_types,atom_pos,coords,makeprim,convention)
using3D()
fig = figure()
ax = fig.add_subplot(111, projection="3d")
ax = plot_3Dconvexhull(ibz,ax,facecolors="coral",alpha=1,edgecolors="black",linewidths = 1)
ax = plot_3Dconvexhull(bz,ax,facecolors="deepskyblue",edgecolors="white",linewidths=1,alpha=0.2)
# output
PyObject <Axes3DSubplot: >
```
"""
function plot_3Dconvexhull end

@doc """
    plot_convexhulls(real_latvecs,atom_types,atom_pos,coords,makeprim,convention;rtol,atol)

Plot the Brillouin and Irreducible Brillouin zone in 2D or 3D.

!!! note "Julia 1.9 and above"
    This function is implemented in a package extension and requires `using PyPlot`.

# Arguments
- `real_latvecs::AbstractMatrix{<:Real}`: the basis of a real-space lattice as
    columns of a matrix.
- `atom_types:AbstractVector{<:Int}`: a list of atom types as integers.
- `atom_pos::AbstractMatrix{<:Real}`: the positions of atoms in the crystal
    structure as columns of a matrix.
- `coords::String`: indicates the positions of the atoms are in \"lattice\" or
    \"Cartesian\" coordinates.
- `makeprim::Bool=false`: make the unit cell primitive before calculating the
    the IBZ if equal to `true`.
- `convention::String="ordinary"`: the convention used to go between real and
    reciprocal space. The two conventions are ordinary (temporal) frequency and
    angular frequency. The transformation from real to reciprocal space is
    unitary if the convention is ordinary.
- `rtol::Real=sqrt(eps(float(maximum(real_latvecs))))` a relative tolerance for
    floating point comparisons.
- `atol::Real=1e-9`: an absolute tolerance for floating point comparisons.

# Returns
- `ax::PyObject`: an updated `ax` with plots of the BZ and IBZ.

# Examples
```julia
ENV["MPLBACKEND"]="qt5agg"
using PyPlot
using SymmetryReduceBZ
real_latvecs = [1 0; .5 1]
atom_types=[0]
atom_pos = Array([0 0]')
coords = "Cartesian"
makeprim = true
convention = "ordinary"
ax=plot_convexhulls(real_latvecs,atom_types,atom_pos,coords,makeprim,convention)
# output
PyObject <AxesSubplot: >
```
"""
function plot_convexhulls end

end #module
