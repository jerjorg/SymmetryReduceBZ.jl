module SymmetryReduceBZPyPlotExt

using PyPlot: figaspect, figure, subplots, using3D, art3D
using SymmetryReduceBZ.Symmetry: calc_bz, calc_ibz
using SymmetryReduceBZ.Utilities: get_uniquefacets, vertices, volume, sortpts2D
import SymmetryReduceBZ.Plotting: plot_convexhulls, plot_2Dconvexhull, plot_3Dconvexhull

function plot_2Dconvexhull(convexhull,
    ax=nothing;facecolor::String="blue",
    alpha::Real=0.3,linewidth::Real=3,edgecolor::String="black")

    bzpts=reduce(hcat, vertices(convexhull))
    x,y=eachrow(bzpts[:,sortpts2D(bzpts)])

    if ax === nothing fig,ax = subplots() end
    ax.fill(x,y, facecolor=facecolor,edgecolor=edgecolor,
        linewidth=linewidth,alpha=alpha)
    ax.set_aspect(1)
    ax
end

function plot_3Dconvexhull(convexhull, ax=nothing;
    facecolors::String="blue",alpha::Real=0.3,linewidths::Real=1,edgecolors::String="black")

    ϵ=0.1*(volume(convexhull))^1/3
    plotrange=[[minimum(v)-ϵ,maximum(v)+ϵ] for v in vertices(convexhull)]

    faces=get_uniquefacets(convexhull)
    edges=deepcopy(faces)
    for i in eachindex(edges)
        push!(edges[i],first(edges[i]))
    end

    if ax === nothing
        using3D()
        fig = figure()
        ax = fig.add_subplot(111, projection="3d")
    end

    p=art3D.Poly3DCollection(faces; alpha, facecolors)
    l=art3D.Line3DCollection(edges; colors=edgecolors, linewidths)

    ax.add_collection3d(p)
    ax.add_collection3d(l)
    ax.auto_scale_xyz(plotrange[1],plotrange[2],plotrange[3])
    ax.set_box_aspect((1, 1, 1))

    return ax
end

function plot_convexhulls(real_latvecs::AbstractMatrix{<:Real},
    atom_types::AbstractVector{<:Int},atom_pos::AbstractMatrix{<:Real},
    coords::String,makeprim::Bool,
    convention::String,ax=nothing;
    rtol::Real=sqrt(eps(float(maximum(real_latvecs)))),atol::Real=1e-9)

    bz = calc_bz(real_latvecs,atom_types,atom_pos,coords,makeprim,convention; rtol,atol)
    ibz = calc_ibz(real_latvecs,atom_types,atom_pos,coords,makeprim,convention; rtol,atol)
    dim = size(real_latvecs,1)

    if ax === nothing
        if dim == 2
            fig,ax=subplots(figsize=figaspect(1)*1.5)
        elseif dim == 3
            using3D()
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
