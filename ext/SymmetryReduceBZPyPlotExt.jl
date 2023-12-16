module SymmetryReduceBZPyPlotExt

ENV["MPLBACKEND"]="qt5agg"

import SymmetryReduceBZ.Symmetry: calc_bz, calc_ibz
import SymmetryReduceBZ.Utilities: get_uniquefacets
import SymmetryReduceBZ.Plotting: plot_convexhulls, plot_2Dconvexhull, plot_3Dconvexhull
import QHull: Chull
import PyCall: PyObject, pyimport
import PyPlot: figaspect, figure, subplots

function plot_2Dconvexhull(convexhull::Chull{<:Real},
    ax::Union{PyObject,Nothing}=nothing;facecolor::String="blue",
    alpha::Real=0.3,linewidth::Real=3,edgecolor::String="black")::PyObject

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

    if ax === nothing fig,ax = subplots() end
    ax.fill(x,y, facecolor=facecolor,edgecolor=edgecolor,
        linewidth=linewidth,alpha=alpha)
    ax.set_aspect(1)
    ax
end

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

    if ax === nothing
        fig = figure()
        ax = fig.add_subplot(111, projection="3d")
    end

    ax.add_collection3d(p)
    ax.add_collection3d(l)
    ax.auto_scale_xyz(plotrange[1],plotrange[2],plotrange[3])
    ax.set_box_aspect((1, 1, 1))

    return ax
end

function plot_convexhulls(real_latvecs::AbstractMatrix{<:Real},
    atom_types::AbstractVector{<:Int},atom_pos::AbstractMatrix{<:Real},
    coords::String,makeprim::Bool,
    convention::String,ax::Union{PyObject,Nothing}=nothing;
    rtol::Real=sqrt(eps(float(maximum(real_latvecs)))),atol::Real=1e-9)::PyObject

    format = "convex hull"
    bz = calc_bz(real_latvecs,atom_types,atom_pos,coords,format,makeprim,
        convention,rtol=rtol,atol=atol)
    ibz = calc_ibz(real_latvecs,atom_types,atom_pos,coords,format,makeprim,
        convention,rtol=rtol,atol=atol)
    dim = size(real_latvecs,1)

    if ax === nothing
        if dim == 2
            fig,ax=subplots(figsize=figaspect(1)*1.5)
        elseif dim == 3
            fig = figure(figsize=figaspect(1)*1.5)
            ax = fig.add_subplot(111, projection="3d")
        else
            throw(ArgumentError("The lattice vectors must be in a 2x2 or 3x3
                matrix."))
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
