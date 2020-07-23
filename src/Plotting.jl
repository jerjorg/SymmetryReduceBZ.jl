module Plotting

include("Symmetry.jl")
include("Lattices.jl")

using .Symmetry
using .Lattices

ENV["MPLBACKEND"]="qt5agg"
using PyCall, PyPlot, QHull, LinearAlgebra

export plot_cells

@doc """
    function sort_points_comparison(a,b,c)

A "less than" function for sorting Cartesian points in 2D.

# Arguments
- `a::AbstractArray{<:Real,1}`: a point in Cartesian coordinates.
- `b::AbstractArray{<:Real,1}`: a point in Cartesian coordinates.
- `c::AbstractArray{<:Real,1}`: the position of the origin in Cartesian
    coordinates.

# Returns
- `Bool`: a boolean that indicates if `a` is less than `b`.


# Examples
```jldoctest
using IBZ
a=[0,0]
b=[1,0]
c=[0.5,0.5]
IBZ.Plotting.sort_points_comparison(a,b,c)
# output
false
```
"""
function sort_points_comparison(a::AbstractArray{<:Real,1},
    b::AbstractArray{<:Real,1},c::AbstractArray{<:Real,1})::Bool
    (ax,ay)=a
    (bx,by)=b
    (cx,cy)=c

    if ((ay - cy) > 0) & ((by - cy) < 0)
        return true
    elseif ((ay - cy) < 0) & ((by - cy) > 0)
        return false
    elseif ((ay - cy) >= 0) & ((by - cy) >= 0)
        return bx > ax
    else #((ay - cy) < 0) & ((by - cy) < 0)
        return ax > bx
    end
end

@doc """
    affine_trans(pts)

Calculate the affine transformation that maps the points to the xy-plane.

# Arguments
- `pts::AbstractArray{<:Real,2}`: an array of Cartesian points as the columns
    of a 2D array. The points must all lie on a plane in 3D.

# Returns
- `M::AbstractArray{<:Real,2}`: the affine transformation matrix that operates
    on points in homogeneous coordinates from the left.

# Examples
```jldoctest
using IBZ
pts = [0.5 0.5 0.5; 0.5 -0.5 0.5; -0.5 0.5 0.5; -0.5 -0.5 0.5]'
IBZ.Plotting.affine_trans(pts)
# output
4×4 Array{Float64,2}:
  0.0  -1.0   0.0  0.5
 -1.0   0.0   0.0  0.5
  0.0   0.0  -1.0  0.5
  0.0   0.0   0.0  1.0
```
"""
function affine_trans(pts::AbstractArray{<:Real,2})::AbstractArray{<:Real,2}
    a,b,c = [pts[:,i] for i=1:3]

    # Create a coordinate system with two vectors lying on the plane the points
    # lie on.
    u = b-a
    v = c-a
    u = u/norm(u)
    v = v - (u⋅v)*u/(u⋅u)
    v = v/norm(v)
    w = u×v

    # Augmented matrix of affine transform
    inv(vcat(hcat([u v w],a),[0 0 0 1]))
end

@doc """
    function mapto_xyplane(pts)

Map Cartesian points embedded in 3D to the xy-plane embedded in 2D.

# Arguments
- `pts::AbstractArray{<:Real,2}`: Cartesian points in 3D as columns of a 2D
    array.

# Returns
- `AbstractArray{<:Real,2}`: Cartesian points in2D as columns of a 2D array.

# Examples
```jldoctest
using IBZ
pts = [0.5 -0.5 0.5; 0.5 -0.5 -0.5; 0.5 0.5 -0.5; 0.5 0.5 0.5]'
IBZ.Plotting.mapto_xyplane(pts)
# output
2×4 Array{Float64,2}:
 0.0  1.0  1.0  0.0
 0.0  0.0  1.0  1.0
```
"""
function mapto_xyplane(pts::AbstractArray{<:Real,2})::AbstractArray{<:Real,2}

    M = affine_trans(pts)
    reduce(hcat,[(M*[pts[:,i]..., 1])[1:2] for i=1:size(pts,2)])
end

@doc """
    function sortslice_perm(xypts,sxypts)

Return the permutation vector that maps Cartesian 2D points `xypts` to `sxypts`.

# Arguments
- `xypts::AbstractArray{<:Real,2}`: 2D Cartesian points as columns of a 2D
    array.
- `sxypts::AbstractArray{<:Real,2}`: sorted 2D Cartesian points as columns of a
    2D array.

# Returns
- `perm::AbstractArray{<:Real,1}`: a permutation vector that sorts an array of
    2D Cartesian coordinates.

# Examples
```jldoctest
using IBZ
xypts = [0 0; 0 1; 1 0; 1 1]'
sxypts = [0 1; 1 1; 1 0; 0 0]'
perm=IBZ.Plotting.sortslice_perm(xypts,sxypts)
xypts[:,perm]
# output
2×4 Array{Int64,2}:
 0  1  1  0
 1  1  0  0
```
"""
function sortslice_perm(xypts::AbstractArray{<:Real,2},
    sxypts::AbstractArray{<:Real,2})::AbstractArray{<:Real,1}
    perm = zeros(Int64,size(xypts,2))
    for i=1:size(xypts,2)
        for j=1:size(xypts,2)
            if isapprox(xypts[:,i],sxypts[:,j])
                perm[j]=i
                continue
            end
        end
    end
    perm
end

@doc """
    function sortpts_perm(pts)

Calculate the permutation vector that sorts Cartesian points embedded in 3D that
    lie on a plane (counter)clockwise with respect to the average of all points.

# Arguments
- `pts::AbstractArray{<:Real,2}`: Cartesian points embedded in 3D that all lie
    on a plane. The points are columns of a 2D array.

# Returns
- `perm::AbstractArray{<:Real,1}`: the permutation vector that orders the points
    clockwise or counterclockwise.

# Examples
```jldoctest
using IBZ
pts = [0.5 -0.5 0.5; 0.5 -0.5 -0.5; 0.5 0.5 -0.5; 0.5 0.5 0.5]'
perm=IBZ.Plotting.sortpts_perm(pts)
pts[:,perm]
# output
3×4 Array{Float64,2}:
 0.5   0.5   0.5   0.5
 0.5   0.5  -0.5  -0.5
 0.5  -0.5  -0.5   0.5
```
"""
function sortpts_perm(pts::AbstractArray{<:Real,2})
    xypts=mapto_xyplane(pts)
    middlept=[sum(xypts[i,:])/size(pts,2) for i=1:2]
    sxypts=sortslices(xypts, lt=(x,y)->sort_points_comparison(x,y,middlept),
        dims=2)
    perm=sortslice_perm(xypts,sxypts)
    return perm
end

@doc """
    get_uniquefacets(ch)

Calculate the unique facets of a convex hull object.
"""
function get_uniquefacets(ch)
    facets = ch.facets
    unique_facets = []
    removed=zeros(Int64,size(facets,1))
    for i=1:size(facets,1)
        if removed[i] == 0
            removed[i]=1
            face=ch.simplices[i]
            for j=i+1:size(facets,1)
                if isapprox(facets[i,:],facets[j,:],rtol=1e-6)
                    removed[j]=1
                    append!(face,ch.simplices[j])
                end
            end
            face = unique(reduce(hcat,face)[:])
            # Order the corners of the face either clockwise or
            # counterclockwise.
            face = face[sortpts_perm(Array(ch.points[face,:]'))]
            append!(unique_facets,[face])
        end
    end
    unique_facets
end

@doc """
    plot_2Dcell(convexhull, ax, color)

Plot a 2D convex hull

# Arguments
- `convexhull::Chull{<:Real}`: a convex hull object.
- `ax::PyObject`: an axes object from matplotlib.
- `color::String`: the face color of the convex hull.

# Returns
- `ax::PyObject`: updated `ax` that includes a plot of the convex hull.

# Examples
```
using IBZ
real_latvecs = [1 0; 0 1]
convention = "ordinary"
bzformat = "convex hull"
bz = calc_bz(real_latvecs,convention,bzformat)

fig,ax=subplots(figsize=figaspect(1)*1.5)
ax = IBZ.Plotting.plot_2Dcell(bz,ax,"deepskyblue");
color="deepskyblue"
IBZ.Plotting.plot_2Dcell(bz,ax,color)
```
"""
function plot_2Dcell(convexhull::Chull{<:Real}, ax::PyObject,
    color::String)::PyObject

    patch=pyimport("matplotlib.patches")
    collections=pyimport("matplotlib.collections")

    bzpts=Array(convexhull.points')
    middlept=[sum(bzpts[i,:])/size(bzpts,2) for i=1:2]
    bzpts=sortslices(bzpts, lt=(x,y)->sort_points_comparison(x,y,middlept),
            dims=2)

    (x,y)=[bzpts[i,:] for i=1:2]
    ax.fill(x,y, facecolor=color,edgecolor="black",linewidth=3)
    ax
end

@doc """
    plot_3Dcell(convexhull, ax, color)

Plot a 3D convex hull

# Arguments
- `convexhull::Chull{<:Real}`: a convex hull object.
- `ax::PyObject`: an axes object from matplotlib.
- `color::String`: the face color of the convex hull.

# Returns
- `ax::PyObject`: updated `ax` that includes a plot of the convex hull.

# Examples
```
using IBZ
real_latvecs = [1 0 0; 0 1 0; 0 0 1]
convention = "ordinary"
bzformat = "convex hull"
bz = calc_bz(real_latvecs,convention,bzformat)
fig = figure(figsize=figaspect(1)*1.5)
ax = fig.add_subplot(111, projection="3d")
ax = IBZ.Plotting.plot_3Dcell(bz,ax,"deepskyblue")
```
"""
function plot_3Dcell(convexhull, ax, color, plotrange=false)

    art3d=pyimport("mpl_toolkits.mplot3d.art3d")

    if plotrange ==false
        ϵ=0.1*convexhull.volume
        plotrange=[[minimum(convexhull.points[:,i])-ϵ,
                maximum(convexhull.points[:,i])+ϵ] for i=1:3]
    end

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

    p=art3d.Poly3DCollection(faces, alpha=0.2, facecolors=color)
    l=art3d.Line3DCollection(edges, colors="black", linewidths=1)

    ax.add_collection3d(p)
    ax.add_collection3d(l)

    #ax.view_init(-45,-45)
    ax.set_xlim3d(plotrange[1]...)
    ax.set_ylim3d(plotrange[2]...)
    ax.set_zlim3d(plotrange[3]...)
    return ax
end

@doc """
    plot_cells(real_latvecs,atom_types,atom_pos,coords,convention,rtol,atol)

Plot the Brillouin and Irreducible Brillouin zone in 2D or 3D.

# Arguments
- `real_latvecs::AbstractArray{<:Real,2}`: the basis of a real-space lattice as
    columns of an array.
- `atom_types:AbstractArray{<:Int,1}`: a list of atom types as integers.
- `atom_pos::AbstractArray{<:Real,2}`: the positions of atoms in the crystal
    structure as columns of an array.
- `coords::String`: indicates the positions of the atoms are in \"lattice\" or
    \"Cartesian\" coordinates.
- `convention::String="ordinary"`: the convention used to go between real and
    reciprocal space. The two conventions are ordinary (temporal) frequency and
    angular frequency. The transformation from real to reciprocal space is
    unitary if the convention is ordinary.
- `rtol::Real=sqrt(eps(float(maximum(real_latvecs))))` a relative tolerance for
    floating point comparisons.
- `atol::Real=0.0`: an absolute tolerance for floating point comparisons.

# Returns
- `(fig,ax)`: the figure and axes Python objects.

# Examples
```
using IBZ
real_latvecs = [1 0; .5 1]
atom_types=[0]
coords = "Cartesian"
atom_pos = Array([0 0]')
convention = "ordinary"
(fig,ax)=plot_cells(real_latvecs,atom_types,atom_pos,coords,convention)
```
"""
function plot_cells(real_latvecs,atom_types,atom_pos,coords,convention,
        rtol::Real=sqrt(eps(float(maximum(real_latvecs)))),atol::Real=0.0)
    mplot3d=pyimport("mpl_toolkits.mplot3d")
    bzformat = "convex hull"
    bz = calc_bz(real_latvecs,convention,bzformat)
    ibzformat = "convex hull"
    ibz = calc_ibz(real_latvecs,atom_types,atom_pos,coords,ibzformat,convention,
        rtol,atol)

    dim = size(real_latvecs,1)
    if dim == 2
        fig,ax=subplots(figsize=figaspect(1)*1.5)
        ax = plot_2Dcell(bz,ax,"deepskyblue")
        ax = plot_2Dcell(ibz,ax,"coral")
    elseif dim == 3
        fig = figure(figsize=figaspect(1)*1.5)
        ax = fig.add_subplot(111, projection="3d")
        ϵ=0.1*bz.volume
        plotrange=[[minimum(bz.points[:,i])-ϵ,
                maximum(bz.points[:,i])+ϵ] for i=1:3]
        ax = plot_3Dcell(bz,ax,"blue",plotrange)
        ax = plot_3Dcell(ibz,ax,"red",plotrange)
    else
        throw(ArgumentError("The lattice vectors must be in a 2x2 or 3x3
            array."))
    end
    (fig,ax)
end

end #module
