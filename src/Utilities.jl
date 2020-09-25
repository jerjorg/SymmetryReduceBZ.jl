module Utilities

import LinearAlgebra: cross, dot, norm
import LinearAlgebra.BLAS: gemv
import Distances: euclidean

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
IBZ.Utilities.affine_trans(pts)
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
    v = v - dot(u,v)*u/dot(u,u)
    v = v/norm(v)
    w = cross(u,v)

    # Augmented matrix of affine transform
    inv(vcat(hcat([u v w],a),[0 0 0 1]))
end

@doc """
    contains(pt,pts,rtol,atol)

Check if an array of points contains a point.

# Arguments
- `pt::AbstractArray{<:Real,1}`: a 1D array of reals
- `pts::AbstractArray{<:Real,2}`: a 2D array of reals. Coordinates of points
    are the columns of the array.
- `rtol::Real=sqrt(eps(float(maximum(pts))))`: a relative tolerance for floating
    point comparisons
- `atol::Real=0.0`: an absolute tolerance for floating point comparisons.

# Returns
- `Bool`: a boolean that indicates the presence or absence of `pt` in `pts`.

# Examples
```jldoctest
pts = Array([1 2; 2 3; 3 4; 4 5]')
pt = [1,2]
contains(pt,pts)
# output
true
```
"""
function contains(pt::AbstractArray{<:Real,1},pts::AbstractArray{<:Real,2},
        rtol::Real=sqrt(eps(float(maximum(pts)))),atol::Real=0.0)::Bool
    any(isapprox(pt,pts[:,i],rtol=rtol,atol=atol) for i=1:size(pts,2))
end

@doc """
    contains(array,arrays,rtol,atol)

Check if an array of arrays contains an array.

# Arguments
- `array`: an array of reals or arbitrary dimension.
- `arrays`: a abstract array of reals of arbitrary dimension.
- `rtol::Real=sqrt(eps(float(maximum(pts))))`: a relative tolerance for floating
    point comparisons.
- `atol::Real=0.0`: an absolute tolerance for floating point comparisons.

# Returns
- `Bool`: a boolean that indicates the presence of absence of `array` in
    `arrays`.

# Examples
```jldoctest
arrays = [[1 2; 2 3], [2 3; 4 5]]
array = [1 2; 2 3]
contains(array, arrays)
# output
true
```
"""
function contains(array,arrays,
    rtol::Real=sqrt(eps(float(maximum(Iterators.flatten(array))))),
    atol::Real=0.0)::Bool
    any(isapprox(array,a,rtol=rtol,atol=atol) for a in arrays)
end

# @doc """
#     contains(pt,pts,rtol,atol)
#
# Check if an array of points contains a point.
#
# # Arguments
# - `item`: a real or abstract array of reals.
# - `pts`:: an abstract array of reals.
# - `rtol::Real=sqrt(eps(float(maximum(pts))))`: a relative tolerance for floating
#     point comparisons
# - `atol::Real=0.0`: an absolute tolerance for floating point comparisons.
#
# # Returns
# - `Bool`: a boolean that indicates the presence of absence of `pt` in `pts`.
#
# # Examples
# ```jldoctest
# pts = Array([1 2; 2 3; 3 4; 4 5]')
# pt = [1,2]
# contains(pt,pts)
# # output
# true
# ```
# """
# contains(item, itr,
#     rtol::Real=sqrt(eps(float(maximum(Iterators.flatten(itr))))),
#     atol::Real=0.0) = any(x -> isapprox(item,x,rtol=rtol,atol=atol),itr)

@doc """
    edgelengths(basis,radius)

Calculate the edge lengths of a parallelepiped circumscribed by a sphere.

# Arguments
- `basis::Array{<:Real,2}`: a 2x2 or 3x3 matrix whose columns give the
    parallelogram or parallelepiped directions, respectively.
- `radius::Real`: the radius of the sphere.
- `rtol::Real=sqrt(eps(float(radius)))`: a relative tolerace for
    floating point comparisons.
- `atol::Real=0.0`: an absolute tolerance for floating point
    comparisons.

# Returns
- `[la,lb,lc]::Array{Float64,1}`: a list of parallelepiped lengths.

# Examples
```jldoctest
using IBZ
basis=Array([1. 0. 0.; 0. 1. 0.; 0. 0. 1.])
radius=3.0
IBZ.Symmetry.edgelengths(basis,radius)
# output
3-element Array{Float64,1}:
 3.0
 3.0
 3.0
```
"""
function edgelengths(basis::Array{<:Real,2}, radius::Real,
        rtol::Real=sqrt(eps(float(radius))), atol::Real=0.0)::Array{Float64,1}

    if radius < 0
        throw(ArgumentError("The radius has to be a positive number."))
    end

    if size(basis) == (2,2)
        (a,b)=[basis[:,i] for i=1:2]
        ax,ay=a
        bx,by=b
        la=2*abs(radius*sqrt(bx^2+by^2)/(ay*bx-ax*by))
        lb=2*abs(radius*sqrt(ax^2+ay^2)/(ay*bx-ax*by))
        return [la,lb]

    elseif size(basis) == (3,3)
        (a,b,c)=[basis[:,i] for i=1:3]
        ax,ay,az=a
        bx,by,bz=b
        cx,cy,cz=c

        la=abs(radius*sqrt((by*cx-bx*cy)^2+(bz*cx-bx*cz)^2+(bz*cy-by*cz)^2)/
            (az*by*cx-ay*bz*cx-az*bx*cy+ax*bz*cy+ay*bx*cz-ax*by*cz))
        lb=abs((radius*sqrt((ay*cx-ax*cy)^2+(az*cx-ax*cz)^2+(az*cy-ay*cz)^2))/
            (az*by*cx-ay*bz*cx-az*bx*cy+ax*bz*cy+ay*bx*cz-ax*by*cz))
        lc=abs(radius*sqrt((ay*bx-ax*by)^2+(az*bx-ax*bz)^2+(az*by-ay*bz)^2)/
            (az*by*cx-ay*bz*cx-az*bx*cy+ax*bz*cy+ay*bx*cz-ax*by*cz))
        return [la,lb,lc]
    else
        throw(ArgumentError("Basis has to be a 2x2 or 3x3 matrix."))
    end
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
    function mapto_xyplane(pts)

Map Cartesian points embedded in 3D on a plane to the xy-plane embedded in 2D.

# Arguments
- `pts::AbstractArray{<:Real,2}`: Cartesian points in 3D as columns of a 2D
    array.

# Returns
- `AbstractArray{<:Real,2}`: Cartesian points in 2D as columns of a 2D array.

# Examples
```jldoctest
using IBZ
pts = [0.5 -0.5 0.5; 0.5 -0.5 -0.5; 0.5 0.5 -0.5; 0.5 0.5 0.5]'
IBZ.Utilities.mapto_xyplane(pts)
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
    sample_circle(basis,radius,offset,rtol,atol)

Sample uniformly within a circle centered about a point.

## Arguments
- `basis::Array{<:Real,2}`: a 2x2 matrix whose columns are the grid generating
    vectors.
- `radius::Real`: the radius of the circle.
- `offset::Array{<:Real,1}=[0.,0.]`: the xy-coordinates of the center of the
    circle.
- `rtol::Real=sqrt(eps(float(radius)))`: a relative tolerace for
    floating point comparisons.
- `atol::Real=0.0`: an absolute tolerance for floating point
    comparisons.

## Returns
- `pts::Array{Float64,2}` a matrix whose columns are sample points in Cartesian
    coordinates.

## Examples
```jldoctest
using IBZ
basis=Array([1. 0.; 0. 1.]')
radius=1.0
offset=[0.,0.]
IBZ.Symmetry.sample_circle(basis,radius,offset)
# output
2×5 Array{Float64,2}:
  0.0  -1.0  0.0  1.0  0.0
 -1.0   0.0  0.0  0.0  1.0
```
"""
function sample_circle(basis::AbstractArray{<:Real,2}, radius::Real,
    offset::AbstractArray{<:Real,1}=[0.,0.],
    rtol::Real=sqrt(eps(float(radius))), atol::Real=0.0)::Array{Float64,2}

    # Put the offset in lattice coordinates and round.
    (o1,o2)=round.(inv(basis)*offset)
    lens=edgelengths(basis,radius)
    n1,n2=round.(lens) .+ 1

    l=0;
    pt=Array{Float64,1}(undef,2)
    pts=Array{Float64,2}(undef,2,Int((2*n1+1)*(2*n2+1)));
    distances=Array{Float64,1}(undef,size(pts,2))
    for (i,j) in Iterators.product((-n1+o1):(n1+o1),(-n2+o2):(n2+o2))
        l+=1
        pt=gemv('N',float(basis),[i,j])
        pts[:,l]=pt
        distances[l]=euclidean(pt,offset)
    end

    return pts[:,distances.<=radius]

end

@doc """
    sample_sphere(basis,radius,offset,rtol,atol)

Sample uniformly within a circle centered about a point.

# Arguments
- `basis::Array{<:Real,2}`: a 3x3 matrix whose columns are the grid generating
    vectors.
- `radius::Real`: the radius of the sphere.
- `offset::Array{<:Real,1}=[0.,0.]`: the xy-coordinates of the center of the
    circle.
- `rtol::Real=sqrt(eps(float(radius)))`: a relative tolerace for
    floating point comparisons.
- `atol::Real=0.0`: an absolute tolerance for floating point
    comparisons.

# Returns
- `pts::Array{Float64,2}` a matrix whose columns are sample points in Cartesian
    coordinates.

# Examples
```jldoctest
using IBZ
basis=Array([1. 0. 0.; 0. 1. 0.; 0. 0. 1.])
radius=1.0
offset=[0.,0.,0.]
IBZ.Symmetry.sample_sphere(basis,radius,offset)
# output
3×7 Array{Float64,2}:
  0.0   0.0  -1.0  0.0  1.0  0.0  0.0
  0.0  -1.0   0.0  0.0  0.0  1.0  0.0
 -1.0   0.0   0.0  0.0  0.0  0.0  1.0
```
"""
function sample_sphere(basis::AbstractArray{<:Real,2}, radius::Real,
    offset::Array{<:Real,1}=[0.,0.,0.], rtol::Real=sqrt(eps(float(radius))),
    atol::Real=0.0)::Array{Float64,2}

    # Put the offset in lattice coordinates and round.
    (o1,o2,o3)=round.(inv(basis)*offset)
    lens=edgelengths(basis,radius)
    n1,n2,n3=round.(lens) .+ 1

    l=0;
    pt=Array{Float64,1}(undef,3)
    pts=Array{Float64,2}(undef,3,Int((2*n1+1)*(2*n2+1)*(2*n3+1)));
    distances=Array{Float64,1}(undef,size(pts,2))
    for (i,j,k) in Iterators.product((-n1+o1):(n1+o1),(-n2+o2):(n2+o2),
                                     (-n3+o3):(n3+o3))
        l+=1
        pt=gemv('N',float(basis),[i,j,k])
        pts[:,l]=pt
        distances[l]=euclidean(pt,offset)
    end

    pts[:,findall(x->(x<radius||isapprox(x,radius,rtol=rtol)),distances)]
end

@doc """
    shoelace(vertices)

Calculate the area of a polygon with the shoelace algorithm.

# Arguments
- `vertices::AbstractArray{<:Real,2}`: the xy-coordinates of the vertices
    of the polygon as the columns of a 2D array.

# Returns
- `<:Real`: the area of the polygon.
"""
function shoelace(vertices)
    xs = vertices[1,:]
    ys = vertices[2,:]
    abs(xs[end]*ys[1] - xs[1]*ys[end] +
        sum([xs[i]*ys[i+1]-xs[i+1]*ys[i] for i=1:(size(vertices,2)-1)]))/2
end

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
IBZ.Utilities.sort_points_comparison(a,b,c)
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
        return ax > bx
    else
        return ax < bx
    end
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
perm=IBZ.Utilities.sortpts_perm(pts)
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
perm=IBZ.Utilities.sortslice_perm(xypts,sxypts)
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
    remove_duplicates(points,rtol,atol)

Remove duplicate points from an array.

# Arguments
- `points::AbstractArray{<:Real,2}`: the points are columns of a 2D array.
- `rtol::Real=sqrt(eps(float(maximum(points))))`: a relative tolerance for
    floating point comparisons.
- `atol::Real=0.0`: an absolume tolerance for floating point comparisons.

# Returns
- `uniquepts::AbstractArray{<:Real,2}`: a 2D array of unique points as columns.

# Examples
```jldoctest
using IBZ
points=Array([1 2; 2 3; 3 4; 1 2]')
IBZ.Utilities.remove_duplicates(points)
# output
2×3 Array{Int64,2}:
 1  2  3
 2  3  4
```
"""
function remove_duplicates(points::AbstractArray{<:Real,2},
    rtol::Real=sqrt(eps(float(maximum(points)))),
    atol::Real=0.0)::AbstractArray{<:Real,2}
    uniquepts=[]
    for i=1:size(points,2)
        pt=points[:,i]
        if !any([isapprox(pt,uniquepts[i],rtol=rtol) for i=1:length(uniquepts)])
            append!(uniquepts,[pt])
        end
    end
    reduce(hcat,uniquepts)
end

end # module
