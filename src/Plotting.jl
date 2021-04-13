module Plotting

ENV["MPLBACKEND"]="qt5agg"
include("Symmetry.jl")
include("Utilities.jl")
import .Symmetry: calc_bz, calc_ibz
import .Utilities: get_uniquefacets, sortpts_perm
import QHull: Chull
import PyCall: PyObject, pyimport
import PyPlot: figaspect, figure, subplots

export plot_convexhulls

@doc """
    plot_2Dconvexhull(convexhull, ax, color)

Plot a 2D convex hull

# Arguments
- `convexhull::Chull{<:Real}`: a convex hull object.
- `ax::PyObject`: an axes object from matplotlib.
- `facecolor::String="blue"`: the color of the area within the convex hull.
- `alpha::Real=0.3`: the transparency of the convex hull.
- `linewidth::Real=3`: the width of the edges.
- `edgecolor::String="black"`: the color of the edges.

# Returns
- `ax::PyObject`: updated `ax` that includes a plot of the convex hull.

# Examples
```jldoctest
import SymmetryReduceBZ.Symmetry: calc_bz, calc_ibz
import SymmetryReduceBZ.Plotting: plot_2Dconvexhull
real_latvecs = [1 0; 0 1]
convention="ordinary"
atom_types=[0]
atom_pos = Array([0 0]')
coords = "Cartesian"
ibzformat = "convex hull"
makeprim=false
bz = calc_bz(real_latvecs,atom_types,atom_pos,coords,ibzformat,makeprim,convention)
ibz = calc_ibz(real_latvecs,atom_types,atom_pos,coords,ibzformat,makeprim,convention)
ax = plot_2Dconvexhull(bz,facecolor="deepskyblue",linewidth=3,edgecolor="cyan",alpha=0.2)
ax = plot_2Dconvexhull(ibz,ax;facecolor="coral",linewidth=3,edgecolor="magenta",alpha=0.4)
# output
PyObject <AxesSubplot:>
```
"""
function plot_2Dconvexhull(convexhull::Chull{<:Real},
    ax::Union{PyObject,Nothing}=nothing;facecolor::String="blue",
    alpha::Real=0.3,linewidth::Real=3,edgecolor::String="black")::PyObject

    patch=pyimport("matplotlib.patches")
    collections=pyimport("matplotlib.collections")

    bzpts=Array(convexhull.points')
    c=[sum(bzpts[i,:])/size(bzpts,2) for i=1:2]
    angles=zeros(size(bzpts,2))
    for i=1:size(bzpts,2)
        (x,y)=bzpts[:,i] - c
        angles[i] = atan(y,x)
    end
    perm = sortperm(angles)
    bzpts = bzpts[:,perm]
    (x,y)=[bzpts[i,:] for i=1:2]

    if ax == nothing fig,ax = subplots() end
    ax.fill(x,y, facecolor=facecolor,edgecolor=edgecolor,
        linewidth=linewidth,alpha=alpha)
    ax.set_aspect(1)
    ax
end

@doc """
    plot_3Dconvexhull(convexhull,ax;color)

Plot a 3D convex hull

# Arguments
- `convexhull::Chull{<:Real}`: a convex hull object.
- `ax::PyObject`: an axes object from matplotlib.
- `facecolors::String="blue"`: the color of the faces of the convex hull.
- `alpha::Real=0.3`: the transparency of the faces of the convex hull.
- `linewidths::Real=1`: the width of the edges of the convex hull.
- `edgecolors::String="black"`: the color of the edges of the convex hull.

# Returns
- `ax::PyObject`: updated `ax` that includes a plot of the convex hull.

# Examples
```jldoctest
import SymmetryReduceBZ.Symmetry: calc_bz, calc_ibz
import SymmetryReduceBZ.Plotting: plot_3Dconvexhull
using PyPlot
real_latvecs = [1 0 0; 0 1 0; 0 0 1]
convention="ordinary"
atom_types=[0]
atom_pos = Array([0 0 0]')
coords = "Cartesian"
bzformat = "convex hull"
makeprim=false
bz = calc_bz(real_latvecs,atom_types,atom_pos,coords,bzformat,makeprim,convention)
ibz = calc_ibz(real_latvecs,atom_types,atom_pos,coords,bzformat,makeprim,convention)
fig = figure()
ax = fig.add_subplot(111, projection="3d")
ax = plot_3Dconvexhull(ibz,ax,facecolors="coral",alpha=1,edgecolors="black",linewidths = 1)
ax = plot_3Dconvexhull(bz,ax,facecolors="deepskyblue",edgecolors="white",linewidths=1,alpha=0.2)
# output
PyObject <Axes3DSubplot:>
```
"""
function plot_3Dconvexhull(convexhull::Chull{<:Real}, ax::Union{PyObject,Nothing}=nothing;
    facecolors::String="blue",alpha::Real=0.3,linewidths::Real=1,edgecolors::String="black")::PyObject

    ϵ=0.1*(convexhull.volume)^1/3
    plotrange=[[minimum(convexhull.points[:,i])-ϵ,
        maximum(convexhull.points[:,i])+ϵ] for i=1:3]

    art3d=pyimport("mpl_toolkits.mplot3d.art3d")

    facesᵢ=get_uniquefacets(convexhull)
    edgesᵢ=deepcopy(facesᵢ)
    for i=1:length(edgesᵢ)
        append!(edgesᵢ[i],edgesᵢ[i][1])
    end

    faces=[convexhull.points[i,:] for i in facesᵢ]
    edges=[convexhull.points[i,:] for i in edgesᵢ]

    # Reshape faces and edges
    faces=[[faces[j][i,:] for i=1:size(faces[j],1)] for j=1:length(faces)]
    edges=[[edges[j][i,:] for i=1:size(edges[j],1)] for j=1:length(edges)]

    p=art3d.Poly3DCollection(faces, alpha=alpha, facecolors=facecolors)
    l=art3d.Line3DCollection(edges, colors=edgecolors,linewidths=linewidths)

    if ax == nothing 
        fig = figure()
        ax = fig.add_subplot(111, projection="3d")
    end

    ax.add_collection3d(p)
    ax.add_collection3d(l)
    ax.auto_scale_xyz(plotrange[1],plotrange[2],plotrange[3])
    ax.set_box_aspect((1, 1, 1))

    return ax
end

@doc """
    plot_convexhulls(real_latvecs,atom_types,atom_pos,coords,makeprim,convention;rtol,atol)

Plot the Brillouin and Irreducible Brillouin zone in 2D or 3D.

# Arguments
- `real_latvecs::AbstractArray{<:Real,2}`: the basis of a real-space lattice as
    columns of an array.
- `atom_types:AbstractArray{<:Int,1}`: a list of atom types as integers.
- `atom_pos::AbstractArray{<:Real,2}`: the positions of atoms in the crystal
    structure as columns of an array.
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
```jldoctest
using SymmetryReduceBZ
real_latvecs = [1 0; .5 1]
atom_types=[0]
atom_pos = Array([0 0]')
coords = "Cartesian"
makeprim = true
convention = "ordinary"
ax=plot_convexhulls(real_latvecs,atom_types,atom_pos,coords,makeprim,convention)
# output
PyObject <AxesSubplot:>
```
"""
function plot_convexhulls(real_latvecs,atom_types,atom_pos,coords,makeprim,
    convention,ax::Union{PyObject,Nothing}=nothing;
    rtol::Real=sqrt(eps(float(maximum(real_latvecs)))),atol::Real=1e-9)::PyObject

    art3d=pyimport("mpl_toolkits.mplot3d.art3d")
    format = "convex hull"
    bz = calc_bz(real_latvecs,atom_types,atom_pos,coords,format,makeprim,
        convention,rtol=rtol,atol=atol)
    ibz = calc_ibz(real_latvecs,atom_types,atom_pos,coords,format,makeprim,
        convention,rtol=rtol,atol=atol)
    dim = size(real_latvecs,1)

    if ax == nothing
        if dim == 2
            fig,ax=subplots(figsize=figaspect(1)*1.5)
        elseif dim == 3
            fig = figure(figsize=figaspect(1)*1.5)
            ax = fig.add_subplot(111, projection="3d")
        else
            throw(ArgumentError("The lattice vectors must be in a 2x2 or 3x3
                array."))
        end
    end

    if dim == 2
        ax = plot_2Dconvexhull(bz,ax,facecolor="deepskyblue")
        ax = plot_2Dconvexhull(ibz,ax,facecolor="coral")
    else # dim == 3
        ax = plot_3Dconvexhull(ibz,ax,facecolors="lightcoral")
        ax = plot_3Dconvexhull(bz,ax,facecolors="deepskyblue")
    end
    ax
end

end #module
