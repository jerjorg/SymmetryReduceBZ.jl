module Utilities

using LinearAlgebra: cross, dot, norm
using Distances: euclidean
using Base.Iterators: flatten
using Polyhedra: Polyhedra, Polyhedron, removehredundancy!, removevredundancy!, hrep, vrep, halfspaces, points, incidentpointindices

@doc """
    volume(hull)

# Arguments
- `hull::Polyhedron`

# Returns
- `vol::Real`

# Examples
"""
volume(p::Polyhedron) = Polyhedra.volume(p)

@doc """
    vertices(hull)

# Arguments
- `hull::Polyhedron`: a convex hull of a polytope

# Returns
- `pts`: a iterator of points (vectors) at the vertices of the convex hull

# Examples
```jldoctest
using SymmetryReduceBZ.Utilities: vertices
pts = Array([1 2; 2 3; 3 4; 4 5]')
pt = [1,2]
contains(pt,pts)
# output
true
```
"""
vertices(hull::Polyhedron) = points(hull)

@doc """
    affine_trans(pts)

Calculate the affine transformation that maps the points to the xy-plane.

# Arguments
- `pts::AbstractMatrix{<:Real}`: Cartesian points as the columns of a matrix.
    The points must all lie on a plane in 3D.

# Returns
- `M::AbstractMatrix{<:Real}`: the affine transformation matrix that operates
    on points in homogeneous coordinates from the left.

# Examples
```jldoctest
using SymmetryReduceBZ
pts = [0.5 0.5 0.5; 0.5 -0.5 0.5; -0.5 0.5 0.5; -0.5 -0.5 0.5]'
SymmetryReduceBZ.Utilities.affine_trans(pts)
# output
4×4 Matrix{Float64}:
  0.0  -1.0   0.0  0.5
 -1.0   0.0   0.0  0.5
  0.0   0.0  -1.0  0.5
  0.0   0.0   0.0  1.0
```
"""
function affine_trans(pts::AbstractMatrix{<:Real})
    a,b,c, = eachcol(pts)

    # Create a coordinate system with two vectors lying on the plane the points
    # lie on.
    u0 = b-a
    v0 = c-a
    u = u0/norm(u0)
    v1 = v0 - dot(u,v0)*u/dot(u,u)
    v = v1/norm(v1)
    w = cross(u,v)

    # Augmented matrix of affine transform
    inv(hvcat((4,1), u, v, w, a, [0 0 0 1]))
end

@doc """
    contains(pt,pts;rtol,atol)

Check if a point is contained in a matrix of points as columns.

# Arguments
- `pt::AbstractVector{<:Real}`: a point whose coordinates are the components of
    a vector.
- `pts::AbstractMatrix{<:Real}`: coordinates of points as columns of a matrix.
- `rtol::Real=sqrt(eps(float(maximum(pts))))`: a relative tolerance for floating
    point comparisons
- `atol::Real=1e-9`: an absolute tolerance for floating point comparisons.

# Returns
- `Bool`: a boolean that indicates the presence or absence of `pt` in `pts`.

# Examples
```jldoctest
using SymmetryReduceBZ.Utilities: contains
pts = Array([1 2; 2 3; 3 4; 4 5]')
pt = [1,2]
contains(pt,pts)
# output
true
```
"""
function contains(pt::AbstractVector{<:Real},pts::AbstractMatrix{<:Real};
        rtol::Real=sqrt(eps(float(maximum(pts)))),atol::Real=1e-9)::Bool
    any(isapprox(pt,pts[:,i]; rtol, atol) for i=1:size(pts,2))
end

@doc """
    contains(array,arrays;rtol,atol)

Check if an array of arrays contains an array.

# Arguments
- `array::AbstractArray`: an array of reals of arbitrary dimension.
- `arrays::AbstractArray`: a nested array of arrays of arbitrary dimension.
- `rtol::Real=sqrt(eps(float(maximum(pts))))`: a relative tolerance for floating
    point comparisons.
- `atol::Real=1e-9`: an absolute tolerance for floating point comparisons.

# Returns
- `Bool`: a boolean that indicates the presence of absence of `array` in
    `arrays`.

# Examples
```jldoctest
using SymmetryReduceBZ.Utilities: contains
arrays = [[1 2; 2 3], [2 3; 4 5]]
array = [1 2; 2 3]
contains(array, arrays)
# output
true
```
"""
function contains(array::AbstractArray,arrays::AbstractArray;
    rtol::Real=sqrt(eps(float(maximum(Iterators.flatten(array))))),
    atol::Real=1e-9)::Bool
    any(isapprox(array,a; rtol, atol) for a in arrays)
end

@doc """
    edgelengths(basis,radius;rtol,atol)

Calculate the edge lengths of a parallelepiped circumscribed by a sphere.

# Arguments
- `basis::AbstractMatrix{<:Real}`: a 2x2 or 3x3 matrix whose columns give the
    parallelogram or parallelepiped directions, respectively.
- `radius::Real`: the radius of the sphere.
- `rtol::Real=sqrt(eps(float(radius)))`: a relative tolerace for
    floating point comparisons.
- `atol::Real=1e-9`: an absolute tolerance for floating point
    comparisons.

# Returns
- `[la,lb,lc]::AbstractVector{<:Real}`: a list of parallelepiped lengths.

# Examples
```jldoctest
using SymmetryReduceBZ
basis=Array([1. 0. 0.; 0. 1. 0.; 0. 0. 1.])
radius=3.0
SymmetryReduceBZ.Utilities.edgelengths(basis,radius)
# output
3-element Vector{Float64}:
 3.0
 3.0
 3.0
```
"""
function edgelengths(basis::AbstractMatrix{<:Real}, radius::Real;
        rtol::Real=sqrt(eps(float(radius))), atol::Real=1e-9)::AbstractVector{<:Real}

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

get_simplex(v::Vector{<:Vector{<:Real}}, i) = v[i]
get_simplex(v::Matrix{<:Real}, i) = v[i,:]


@doc """
    get_uniquefacets(ch)

Calculate the unique facets of a convex hull.

!!! note "QHull.jl package extension"
    This function is available through a package extension of
    [QHull.jl](https://github.com/JuliaPolyhedra/QHull.jl).
    After installing Python, SciPy, and QHull.jl, do `using QHull` to load it.

# Arguments
- `ch::Polyhedron`: a convex hull in 3D from `QHull`.

# Returns
- `unique_facets::Vector{Vector{Int64}}`: a nested list of the
    indices of points that lie on each face. For example, the points that lie on
    the first face are `ch.points[unique_facets[1],:]`.

# Examples
```jldoctest
using QHull
using SymmetryReduceBZ.Utilities: get_uniquefacets
using SymmetryReduceBZ.Symmetry: calc_bz
real_latvecs = [1 0 0; 0 1 0; 0 0 1]
atom_types = [0]
atom_pos = Array([0 0 0]')
coordinates = "Cartesian"
makeprim = false
convention = "ordinary"
bz = calc_bz(real_latvecs,atom_types,atom_pos,coordinates,makeprim,convention)
get_uniquefacets(bz)
# output
6-element Vector{Vector{Int32}}:
 [1, 2, 3, 4]
 [7, 2, 3, 5]
 [6, 4, 3, 5]
 [7, 2, 1, 8]
 [6, 4, 1, 8]
 [8, 7, 5, 6]
```
"""
function get_uniquefacetsindices(poly::Polyhedron)
    removehredundancy!(poly)
    removevredundancy!(poly)
    map(eachindex(halfspaces(poly))) do i
        id = incidentpointindices(poly, i)
        permute!(id, sortpts_perm(mapreduce(j -> get(poly, j), hcat, id)))
        return id
    end
end

function get_uniquefacets(poly::Polyhedron)
    map(get_uniquefacetsindices(poly)) do j
        map(i -> get(poly, i), j)
    end
end

@doc """
    function mapto_xyplane(pts)

Map Cartesian points embedded in 3D on a plane to the xy-plane embedded in 2D.

# Arguments
- `pts::AbstractMatrix{<:Real}`: Cartesian points embedded in 3D as columns of a
    matrix.

# Returns
- `AbstractMatrix{<:Real}`: Cartesian points in 2D as columns of a matrix.

# Examples
```jldoctest
using SymmetryReduceBZ
pts = [0.5 -0.5 0.5; 0.5 -0.5 -0.5; 0.5 0.5 -0.5; 0.5 0.5 0.5]'
SymmetryReduceBZ.Utilities.mapto_xyplane(pts)
# output
2×4 Matrix{Float64}:
 0.0  1.0  1.0  0.0
 0.0  0.0  1.0  1.0
```
"""
function mapto_xyplane(pts::AbstractMatrix{<:Real})

    M = affine_trans(pts)
    mapslices(pt -> first(M*vcat(pt, 1), 2), pts, dims=1)
end

@doc """
    sample_circle(basis,radius,offset;rtol,atol)

Sample uniformly within a circle centered about a point.

## Arguments
- `basis::AbstractMatrix{<:Real}`: a 2x2 matrix whose columns are the grid
    generating vectors.
- `radius::Real`: the radius of the circle.
- `offset::AbstractVector{<:Real}=[0.,0.]`: the xy-coordinates of the center of
    the circle.
- `rtol::Real=sqrt(eps(float(radius)))`: a relative tolerace for floating point
    comparisons.
- `atol::Real=1e-9`: an absolute tolerance for floating point comparisons.

## Returns
- `pts::AbstractMatrix{<:Real}` a matrix whose columns are sample points in Cartesian
    coordinates.

## Examples
```jldoctest
using SymmetryReduceBZ
basis=Array([1. 0.; 0. 1.]')
radius=1.0
offset=[0.,0.]
SymmetryReduceBZ.Utilities.sample_circle(basis,radius,offset)
# output
2×5 Matrix{Float64}:
  0.0  -1.0  0.0  1.0  0.0
 -1.0   0.0  0.0  0.0  1.0
```
"""
function sample_circle(basis::AbstractMatrix{<:Real}, radius::Real,
    offset::AbstractVector{<:Real}=[0.,0.];
    rtol::Real=sqrt(eps(float(radius))), atol::Real=1e-9)

    # Put the offset in lattice coordinates and round.
    (o1,o2)=round.(inv(basis)*offset)
    lens=edgelengths(basis,radius)
    n1,n2=round.(lens) .+ 1

    l=0;
    pts=Matrix{eltype(basis)}(undef,2,Int((2*n1+1)*(2*n2+1)));
    distances=Vector{float(eltype(basis))}(undef,size(pts,2))
    for (i,j) in Iterators.product((-n1+o1):(n1+o1),(-n2+o2):(n2+o2))
        l+=1
        pt= basis * [i, j]
        pts[:,l] = pt
        distances[l]=euclidean(pt,offset)
    end

    return pts[:,distances.<=radius]

end

@doc """
    sample_sphere(basis,radius,offset;rtol,atol)

Sample uniformly within a sphere centered about a point.

# Arguments
- `basis::AbstractMatrix{<:Real}`: a 3x3 matrix whose columns are the grid generating
    vectors.
- `radius::Real`: the radius of the sphere.
- `offset::AbstractVector{<:Real}=[0.,0.]`: the xy-coordinates of the center of the
    circle.
- `rtol::Real=sqrt(eps(float(radius)))`: a relative tolerace for
    floating point comparisons.
- `atol::Real=1e-9`: an absolute tolerance for floating point
    comparisons.

# Returns
- `pts::AbstractMatrix{<:Real}` a matrix whose columns are sample points in Cartesian
    coordinates.

# Examples
```jldoctest
using SymmetryReduceBZ
basis=Array([1. 0. 0.; 0. 1. 0.; 0. 0. 1.])
radius=1.0
offset=[0.,0.,0.]
SymmetryReduceBZ.Utilities.sample_sphere(basis,radius,offset)
# output
3×7 Matrix{Float64}:
  0.0   0.0  -1.0  0.0  1.0  0.0  0.0
  0.0  -1.0   0.0  0.0  0.0  1.0  0.0
 -1.0   0.0   0.0  0.0  0.0  0.0  1.0
```
"""
function sample_sphere(basis::AbstractMatrix{<:Real}, radius::Real,
    offset::AbstractVector{<:Real}=[0.,0.,0.]; rtol::Real=sqrt(eps(float(radius))),
    atol::Real=1e-9)

    # Put the offset in lattice coordinates and round.
    (o1,o2,o3)=round.(inv(basis)*offset)
    lens=edgelengths(basis,radius)
    n1,n2,n3=round.(lens) .+ 1

    l=0;
    pts=Matrix{eltype(basis)}(undef,3,Int((2*n1+1)*(2*n2+1)*(2*n3+1)));
    distances=Vector{float(eltype(basis))}(undef,size(pts,2))
    for (i,j,k) in Iterators.product((-n1+o1):(n1+o1),(-n2+o2):(n2+o2),
                                     (-n3+o3):(n3+o3))
        l+=1
        pt= basis * [i,j,k]
        pts[:,l]=pt
        distances[l]=euclidean(pt,offset)
    end

    pts[:,findall(x->(x<radius||isapprox(x,radius,rtol=rtol)),distances)]
end

@doc """
    shoelace(vertices)

Calculate the area of a polygon with the shoelace algorithm.

# Arguments
- `vertices::AbstractMatrix{<:Real}`: the xy-coordinates of the vertices
    of the polygon as the columns of a matrix.

# Returns
- `area::Real`: the area of the polygon.

# Examples
```jldoctest
using SymmetryReduceBZ.Utilities: shoelace
pts = [0 0 1; -1 1 0]
shoelace(pts)
# output
1.0
```
"""
function shoelace(vertices::AbstractMatrix{<:Real})
    @assert size(vertices, 1) == 2  "shoelace only works for polygons"
    verts = vertices[:,sortpts2D(vertices)]
    xs = verts[begin,:]
    ys = verts[end,:]
    A = (ys[end]+ys[begin])*(xs[end]-xs[begin])
    for i in axes(verts,2)[begin:end-1]
        A += (ys[i]+ys[i+1])*(xs[i]-xs[i+1])
    end
    abs(A)/2
end

@doc """
    function sortpts2D(pts)

Calculate the permutation vector that sorts 2D Cartesian points counterclockwise with
    respect to the average of the points.

# Arguments
- `pts::AbstractMatrix{<:Real}`: Cartesian points in 2D.

# Returns
- `perm::AbstractVector{<:Real}`: the permutation vector that orders the points
    clockwise or counterclockwise.
```
"""
function sortpts2D(pts::AbstractMatrix{<:Real})
    @assert size(pts,1) == 2 "sortpts2D only accepts polygons"
    cx,cy=map(x->sum(x)/size(pts,2),eachrow(pts))
    angles=map(eachcol(pts)) do v
        vx,vy = v
        atan(vy-cy,vx-cx)
    end
    return sortperm(angles)
end

@doc """
    function sortpts_perm(pts)

Calculate the permutation vector that sorts Cartesian points embedded in 3D that
    lie on a plane (counter)clockwise with respect to the average of all points.

# Arguments
- `pts::AbstractMatrix{<:Real}`: Cartesian points embedded in 3D that all lie
    on a plane. The points are columns of a matrix.

# Returns
- `::AbstractVector{<:Real}`: the permutation vector that orders the points
    clockwise or counterclockwise.

# Examples
```jldoctest
using SymmetryReduceBZ.Utilities: sortpts_perm
pts = [0.5 -0.5 0.5; 0.5 -0.5 -0.5; 0.5 0.5 -0.5; 0.5 0.5 0.5]'
perm=sortpts_perm(pts)
pts[:,perm]
# output
3×4 Matrix{Float64}:
  0.5   0.5   0.5  0.5
 -0.5  -0.5   0.5  0.5
  0.5  -0.5  -0.5  0.5
```
"""
function sortpts_perm(pts::AbstractMatrix{<:Real})
    xypts=mapto_xyplane(pts)
    sortpts2D(xypts)
end

@doc """
    unique_points(points;rtol,atol)

Remove duplicate points.

# Arguments
- `points::AbstractMatrix{<:Real}`: the points are columns of a matrix.
- `rtol::Real=sqrt(eps(float(maximum(flatten(points)))))`: a relative tolerance
    for floating point comparisons.
- `atol::Real=1e-9`: an absolute tolerance for floating point comparisons.

# Returns
- `uniquepts::AbstractMatrix{<:Real}`: the unique points as columns of a matrix.

# Examples
```jldoctest
using SymmetryReduceBZ
points=Array([1 2; 2 3; 3 4; 1 2]')
SymmetryReduceBZ.Utilities.unique_points(points)
# output
2×3 Matrix{Int64}:
 1  2  3
 2  3  4
```
"""
function unique_points(points::AbstractMatrix{<:Real};
    rtol::Real=sqrt(eps(float(maximum(flatten(points))))),
    atol::Real=1e-9)

    uniquepts=similar(points)
    j1 = firstindex(uniquepts, 2)
    numpts = 0
    for i=axes(points,2)
        if !any(j -> isapprox(points[:,i],uniquepts[:,j]; rtol, atol),
                j1:j1+numpts-1)
            numpts += 1
            uniquepts[:,j1+numpts-1] = points[:,i]
        end
    end
    uniquepts[:,j1:j1+numpts-1]
end

@doc """
    remove_duplicates(points;rtol,atol)

Remove duplicates from an array.

# Arguments
- `points::AbstractVector`: items in a vector, which can be floats or arrays.
- `rtol::Real=sqrt(eps(float(maximum(points))))`: relative tolerance.
- `atol::Real=1e-9`: absolute tolerance.

# Returns
- `uniquepts::AbstractVector`: an vector with only unique elements.

# Examples
```jldoctest
using SymmetryReduceBZ.Utilities: remove_duplicates
test = [1.,1.,2,2,]
remove_duplicates(test)
# output
2-element Vector{Float64}:
 1.0
 2.0
```
"""
function remove_duplicates(points::AbstractVector;
    rtol::Real=sqrt(eps(float(maximum(flatten(points))))),
    atol::Real=1e-9)
    uniquepts=similar(points)
    i1 = firstindex(uniquepts)
    npts = 0
    for i=eachindex(points)
        pt=points[i]
        if !any(i -> isapprox(pt,uniquepts[i]; rtol, atol), 1:npts)
            npts += 1
            uniquepts[i1+npts-1] = pt
        end
    end
    uniquepts[i1:i1+npts-1]
end

@doc """
    points_in_ball(points,radius,offset,rtol=sqrt(eps(float(radius))),atol=1e-9)

Calculate the points within a ball (circle, sphere, ...).

# Arguments
- `points::AbstractMatrix{<:Real}`: points in Cartesian coordinates as columns of a matrix.
- `radius::Real`: the radius of the ball.
- `offset::AbstractVector{<:Real}`: the location of the center of the ball in Cartesian coordinates.
- `rtol::Real=sqrt(eps(float(radius)))`: a relative tolerance for floating point comparisons.
- `atol::Real=1e-9`: an absolute tolerance.

# Returns
- `ball_points::AbstractVector{<:Int}`: the indices of points in `points` within the ball.

# Examples
```jldoctest
using SymmetryReduceBZ.Utilities: points_in_ball
points = [0.0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.25 0.3 0.35 0.4 0.45 0.5 0.3 0.35 0.4 0.45 0.5 0.35 0.4 0.45 0.5 0.4 0.45 0.5 0.45 0.5 0.5; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.15 0.15 0.15 0.15 0.15 0.15 0.15 0.15 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.25 0.25 0.25 0.25 0.25 0.25 0.3 0.3 0.3 0.3 0.3 0.35 0.35 0.35 0.35 0.4 0.4 0.4 0.45 0.45 0.5]
radius = 0.1
offset = [0,0]
points_in_ball(points,radius,offset)
# output
4-element Vector{Int64}:
  1
  2
  3
 12
```
"""
function points_in_ball(points::AbstractMatrix{<:Real},radius::Real,
    offset::AbstractVector{<:Real};rtol::Real=sqrt(eps(float(radius))),
    atol::Real=1e-9)

    ball_points = zeros(Int,size(points,2))
    count = 0
    for i=axes(points,2)
        if (norm(points[:,i] - offset) < radius) ||
            isapprox(norm(points[:,i] - offset),radius; rtol, atol)
            count+=1
            ball_points[count] = i
        end
    end
    ball_points[1:count]
end

end # module
