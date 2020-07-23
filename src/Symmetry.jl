module Symmetry

using Distances, LinearAlgebra, CDDLib, QHull, Polyhedra, Combinatorics

include("Lattices.jl")
using .Lattices

export mapto_unitcell, calc_pointgroup, calc_spacegroup, sampleCircle,
    sampleSphere, calc_bz, calc_ibz

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

    if radius < 0 | isapprox(radius, 0, rtol=rtol, atol=atol)
        error("The radius has to be a positive number.")
    end

    if size(basis,1) == 2
        (a,b)=[basis[:,i] for i=1:2]
        ax,ay=a
        bx,by=b
        la=2*abs(radius*sqrt(bx^2+by^2)/(ay*bx-ax*by))
        lb=2*abs(radius*sqrt(ax^2+ay^2)/(ay*bx-ax*by))
        return [la,lb]

    elseif size(basis,1) == 3
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
    sampleCircle(basis,radius,offset,rtol,atol)

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
IBZ.Symmetry.sampleCircle(basis,radius,offset)
# output
2×5 Array{Float64,2}:
  0.0  -1.0  0.0  1.0  0.0
 -1.0   0.0  0.0  0.0  1.0
```
"""
function sampleCircle(basis::AbstractArray{<:Real,2}, radius::Real,
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
        pt=BLAS.gemv('N',float(basis),[i,j])
        pts[:,l]=pt
        distances[l]=euclidean(pt,offset)
    end

    return pts[:,distances.<=radius]

end

@doc """
    sampleSphere(basis,radius,offset,rtol,atol)

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
IBZ.Symmetry.sampleSphere(basis,radius,offset)
# output
3×7 Array{Float64,2}:
  0.0   0.0  -1.0  0.0  1.0  0.0  0.0
  0.0  -1.0   0.0  0.0  0.0  1.0  0.0
 -1.0   0.0   0.0  0.0  0.0  0.0  1.0
```
"""
function sampleSphere(basis::AbstractArray{<:Real,2}, radius::Real,
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
        pt=BLAS.gemv('N',float(basis),[i,j,k])
        pts[:,l]=pt
        distances[l]=euclidean(pt,offset)
    end

    pts[:,findall(x->(x<radius||isapprox(x,radius,rtol=rtol)),distances)]
end

@doc """
    calc_pointgroup(latvecs,rtol,atol)

Calculate the point group of lattice in 2D or 3D.

# Arguments
- `latvecs::AbstractArray{<:Real,2}`: the basis of the lattice as columns of an
    array.
- `rtol::Real=sqrt(eps(float(maximum(real_latvecs))))`: a relative tolerance for
    floating point comparisons. It is used to compare lengths of vectors and the
    volumes of primitive cells.
- `atol::Real=0.0`: an absolute tolerance for floating point comparisons. It is
    used to compare lengths of vectors and the volumes of primitive cells.

# Returns
- `pointgroup::Array{Array{Float64,2},1}`: the point group of the lattice. The
    operators operate on points in Cartesian coordinates.

# Examples
```jldoctest
using IBZ
basis = [1 0; 0 1]
IBZ.Symmetry.calc_pointgroup(basis)
# output
8-element Array{Array{Float64,2},1}:
 [0.0 -1.0; -1.0 0.0]
 [0.0 -1.0; 1.0 0.0]
 [-1.0 0.0; 0.0 -1.0]
 [1.0 0.0; 0.0 -1.0]
 [-1.0 0.0; 0.0 1.0]
 [1.0 0.0; 0.0 1.0]
 [0.0 1.0; -1.0 0.0]
 [0.0 1.0; 1.0 0.0]
```
"""
function calc_pointgroup(latvecs::AbstractArray{<:Real,2},
    rtol::Real=sqrt(eps(float(maximum(latvecs)))),
    atol::Real=0.0)::Array{Array{Float64,2},1}

    dim=size(latvecs,1)
    latvecs = minkowski_reduce(latvecs)
    radius = maximum([norm(latvecs[:,i]) for i=1:dim])
    if dim==2
        pts=sampleCircle(latvecs,radius,[0,0],rtol,atol)
    elseif dim==3
        pts=sampleSphere(latvecs,radius,[0,0,0],rtol,atol)
    else
        throw(ArgumentError("The lattice basis must be a 2x2 or 3x3 array."))
    end

    normsᵢ=[norm(latvecs[:,i]) for i=1:dim]
    sizeᵢ=abs(det(latvecs))
    pointgroup=Array{Array{Float64,2}}([])
    inv_latvecs=inv(latvecs)
    for perm=permutations(1:size(pts,2),dim)
        latvecsᵣ=pts[:,perm]
        normsᵣ=[norm(latvecsᵣ[:,i]) for i=1:dim]
        sizeᵣ=abs(det(latvecsᵣ))
        # Point operation preserves length of lattice vectors and
        # volume on primitive cell.
        if (isapprox(normsᵢ,normsᵣ,rtol=rtol,atol=atol) &
            isapprox(sizeᵢ,sizeᵣ,rtol=rtol,atol=atol))
            op=latvecsᵣ*inv_latvecs
            # Point operators are orthogonal.
            if isapprox(op'*op,I,rtol=rtol,atol=atol)
                append!(pointgroup,[op])
            end
        end
    end
    pointgroup
end


@doc """
    mapto_unitcell(pt,real_latvecs,inv_latvecs,coords,rtol)

Map a point to the first unit (primitive) cell.

# Arguments
- `pt::AbstractArray{<:Real,1}`: a point in lattice or Cartesian coordinates.
- `latvecs::AbstractArray{<:Real,2}`: the basis vectors of the lattice as columns
    of an array.
- `inv_latvecs::AbstractArray{<:Real,2}`: the inverse of the matrix of that
    contains the lattice vectors.
- `coords::String`: indicates whether `pt` is in \"Cartesian\" or \"lattice\"
    coordinates.
- `rtol::Real=sqrt(eps(float(maximum(inv_latvecs))))`: a relative tolerance for
    floating point comparisons. Finite precision errors creep up when `pt` is
    transformed to lattice coordinates because the transformation requires
    calculating a matrix inverse. The components of the point in lattice
    coordinates are checked to ensure that values close to 1 are equal to 1.

# Returns
- `AbstractArray{<:Real,1}`: a translationally equivalent point to `pt` in the
    first unit cell in the same coordinates.

#Examples
```jldoctest
using IBZ
real_latvecs = [0 1 2; 0 -1 1; 1 0 0]
inv_latvecs=inv(real_latvecs)
pt=[1,2,3.2]
coords = "Cartesian"
IBZ.Symmetry.mapto_unitcell(pt,real_latvecs,inv_latvecs,coords)
# output
3-element Array{Float64,1}:
 0.0
 0.0
 0.20000000000000018
```
"""
function mapto_unitcell(pt::AbstractArray{<:Real,1},
    latvecs::AbstractArray{<:Real,2},inv_latvecs::AbstractArray{<:Real,2},
    coords::String,rtol::Real=sqrt(eps(float(maximum(inv_latvecs))))
    )::Array{Float64,1}
    if coords == "Cartesian"
        Array{Float64,1}(latvecs*[isapprox(mod(c,1),1,rtol=rtol) ? 0 : mod(c,1)
            for c=inv_latvecs*pt])
    elseif coords == "lattice"
        Array{Float64,1}([isapprox(mod(c,1),1,rtol=rtol) ? 0 : mod(c,1) for
            c=pt])
    else
        throw(ArgumentError("Allowed coordinates of the point are \"Cartesian\"
                and \"lattice\"."))
    end
end

@doc """
    calc_spacegroup(real_latvecs,atom_types,atom_pos,coords,rtol,atol)

Calculate the space group of a crystal structure.

# Arguments
- `real_latvecs::AbstractArray{<:Real,2}`: the basis of the lattice as columns
    of an array.
- `atom_types::AbstractArray{<:Int,1}`: a list of atom types as integers.
- `atom_pos::AbstractArray{<:Real,2}`: the positions of atoms in the crystal
    structure as columns of an array.
- `coords::String`: indicates the positions of the atoms are in \"lattice\" or
    \"Cartesian\" coordinates.
- `rtol::Real=sqrt(eps(float(maximum(real_latvecs))))` a relative tolerance for
    floating point comparisons.
- `atol::Real=0.0`: an absolute tolerance for floating point comparisons.

# Returns
- `spacegroup`: the space group of the crystal structure. The first element of
    `spacegroup` is a list of fractional translations, and the second element is
    a list of point operators. The translations are in Cartesian coordinates,
    and the operators operate on points in Cartesian coordinates.

# Examples
```jldoctest
using IBZ
real_latvecs = Array([1 0; 2 1]')
atom_types = [0, 1]
atom_pos = Array([0 0; 0.5 0.5]')
coords = "Cartesian"
IBZ.Symmetry.calc_spacegroup(real_latvecs,atom_types,atom_pos,coords)
# output
(Any[[0.0, 0.0], [0.0, 0.0]], Any[[1.0 0.0; 0.0 1.0], [0.0 1.0; 1.0 0.0]])
```
"""
function calc_spacegroup(real_latvecs::AbstractArray{<:Real,2},
    atom_types::AbstractArray{<:Int,1},atom_pos::AbstractArray{<:Real,2},
    coords::String,rtol::Real=sqrt(eps(float(maximum(real_latvecs)))),
    atol::Real=0.0)

    if length(atom_types) != size(atom_pos,2)
        throw(ArgumentError("The number of atom types and positions must be the
            same."))
    end

    real_latvecs = minkowski_reduce(real_latvecs)
    inv_latvecs = inv(real_latvecs)
    pointgroup = calc_pointgroup(real_latvecs,rtol,atol)
    numatoms = length(atom_types)

    # Map points to the primitive cell.
    atom_pos = reduce(hcat,[mapto_unitcell(atom_pos[:,i],real_latvecs,
                inv_latvecs,coords,rtol) for i=1:numatoms])

    # Place atom positions in Cartesian coordinates.
    if coords == "lattice"
        atom_pos = real_latvecs*atom_pos
    end

    ops_spacegroup=[]
    trans_spacegroup=[]
    kindᵣ=atom_types[1]

    # Atoms of the same type as the first.
    same_atoms=findall(x->x==kindᵣ,atom_types)
    for op in pointgroup
        # Use the first atom to calculate allowed fractional translations.
        posᵣ=op*atom_pos[:,1]

        for atomᵢ=same_atoms
            # Calculate a candidate fractional translation.
            ftrans = mapto_unitcell(atom_pos[:,atomᵢ]-posᵣ,real_latvecs,
                inv_latvecs,"Cartesian")

            # Check that all atoms land on the lattice when rotated and
            # translated.
            mapped = false
            for atomⱼ = 1:numatoms
                kindⱼ = atom_types[atomⱼ]
                posⱼ = op*atom_pos[:,atomⱼ] + ftrans
                mapped = any([kindⱼ == atom_types[i] &&
                        isapprox(posⱼ,atom_pos[:,i],rtol=rtol,atol=atol) ?
                            true : false for i=1:numatoms])
                if !mapped continue end
            end
            if mapped
                append!(ops_spacegroup,[op])
                append!(trans_spacegroup,[ftrans])
            end
        end
    end
    (trans_spacegroup,ops_spacegroup)
end

@doc """
    calc_bz(real_latvecs, convention, vertsOrHrep=true, eps=1e-9)

Calculate the Brillouin zone for the given real-spcae lattice.

# Arguments
- `real_latvecs::AbstractArray{<:Real,2}`: the real-space lattice vectors or
    primitive translation vectors as columns of a 3x3 array.
- `convention::String="ordinary"`: the convention used to go between real and
    reciprocal space. The two conventions are ordinary (temporal) frequency and
    angular frequency. The transformation from real to reciprocal space is
    unitary if the convention is ordinary.
- `bzformat::String`: the format of the Brillouin zone. Options include
    \"convex-hull\" and \"half-space\".

# Returns
- `bz`: the vertices or half-space representation of the Brillouin zone
    depending on the value of `vertsOrHrep`.

# Examples
```jldoctest
using IBZ
real_latvecs=[1 0; 0 1]
convention="ordinary"
bzformat = "convex hull"
IBZ.Symmetry.calc_bz(real_latvecs,convention,bzformat)
# output
Convex Hull of 4 points in 2 dimensions
Hull segment vertex indices:
[3, 2, 1, 4]
Points on convex hull in original order:

[0.5 0.5; 0.5 -0.5; -0.5 -0.5; -0.5 0.5]
```
"""
function calc_bz(real_latvecs::AbstractArray{<:Real,2},convention::String,
    bzformat::String)

    recip_latvecs = minkowski_reduce(get_recip_latvecs(real_latvecs,convention))

    if size(real_latvecs) == (2,2)
        latpts = reduce(hcat,[recip_latvecs*[i,j] for
            (i,j)=Iterators.product(-2:2,-2:2)])
    else
        latpts = reduce(hcat,[recip_latvecs*[i,j,k] for
            (i,j,k)=Iterators.product(-2:2,-2:2,-2:2)])
    end

    distances = [norm(latpts[:,i]) for i=1:size(latpts,2)] .^2 ./2
    bz = HalfSpace(latpts[:,2],distances[2])
    for i=3:size(latpts,2)
        bz = bz ∩ HalfSpace(latpts[:,i],distances[i])
    end
    verts = reduce(hcat,points(polyhedron(bz,CDDLib.Library())))
    bzvolume = chull(Array(verts')).volume

    if !(bzvolume ≈ abs(det(recip_latvecs)))
        error("The area or volume of the Brillouin zone is incorrect.")
    end

    if bzformat == "half-space"
        bz
    elseif bzformat == "convex hull"
        chull(Array(verts'))
    else
        throw(ArgumentError("Formats for the BZ include \"half-space\" and
            \"convex hull\"."))
    end
end

@doc """
    calc_ibz(real_latvecs,atom_types,atom_pos,coords,ibzformat,convention,rtol,
        atol)

Calculate the irreducible Brillouin zone of a crystal structure in 2D or 3D.

# Arguments
- `real_latvecs::AbstractArray{<:Real,2}`: the basis of a real-space lattice as
    columns of an array.
- `atom_types:AbstractArray{<:Int,1}`: a list of atom types as integers.
- `atom_pos::AbstractArray{<:Real,2}`: the positions of atoms in the crystal
    structure as columns of an array.
- `coords::String`: indicates the positions of the atoms are in \"lattice\" or
    \"Cartesian\" coordinates.
- `ibzformat::String`: the format of the irreducible Brillouin zone. Options
    include \"convex-hull\" and \"half-space\".
- `convention::String="ordinary"`: the convention used to go between real and
    reciprocal space. The two conventions are ordinary (temporal) frequency and
    angular frequency. The transformation from real to reciprocal space is
    unitary if the convention is ordinary.
- `rtol::Real=sqrt(eps(float(maximum(real_latvecs))))` a relative tolerance for
    floating point comparisons.
- `atol::Real=0.0`: an absolute tolerance for floating point comparisons.

# Returns
- `ibz`: the irreducible Brillouin zone as a convex hull or intersection of
    half-spaces.

# Examples
```jldoctest
using IBZ
real_latvecs = [1 0; 0 1]
convention="ordinary"
atom_types=[0]
atom_pos = Array([0 0]')
coords = "Cartesian"
ibzformat = "convex hull"
IBZ.Symmetry.calc_ibz(real_latvecs,atom_types,atom_pos,coords,ibzformat,
    convention)
# output
Convex Hull of 3 points in 2 dimensions
Hull segment vertex indices:
[1, 2, 3]
Points on convex hull in original order:

[0.0 0.0; 0.5 0.0; 0.5 0.5]
```
"""
function calc_ibz(real_latvecs::AbstractArray{<:Real,2},
        atom_types::AbstractArray{<:Int,1},atom_pos::AbstractArray{<:Real,2},
        coords::String,ibzformat::String,convention::String="ordinary",
        rtol::Real=sqrt(eps(float(maximum(real_latvecs)))),
        atol::Real=0.0)::Chull{<:Real}
    pointgroup = calc_spacegroup(real_latvecs,atom_types,atom_pos,coords)[2]
    copy_pg = deepcopy(pointgroup)
    bzformat = "half-space"
    bz = calc_bz(real_latvecs,convention,bzformat)
    bz_vertices = collect(points(polyhedron(bz,CDDLib.Library())))
    ibz = bz
    while length(copy_pg) > 0
        op = pop!(copy_pg)
        for v=bz_vertices
            vʳ=op*v
            if !isapprox(vʳ,v,rtol=rtol,atol=atol)
                a = vʳ-v
                ibz = ibz ∩ HalfSpace(a,0)
                break
            end
        end
    end
    ibz_vertices = reduce(hcat,points(polyhedron(ibz,CDDLib.Library())))
    bz_vertices = reduce(hcat,bz_vertices)
    bzvolume = chull(Array(bz_vertices')).volume
    ibzvolume = chull(Array(ibz_vertices')).volume
    reduction = bzvolume/ibzvolume
    if !(ibzvolume ≈ bzvolume/size(pointgroup,1))
        error("The area or volume of the irreducible Brillouin zone is
            incorrect.")
    end
    if ibzformat == "half-space"
        ibz
    elseif ibzformat == "convex hull"
        chull(Array(ibz_vertices'))
    else
        throw(ArgumentError("Formats for the IBZ include \"half-space\" and
            \"convex hull\"."))
    end
end

num_ops = [8, 4, 12, 2, 2, 4, 8, 12, 16, 24, 48]
lat_types = ["square", "rectangular", "triangular", "oblique", "triclinic",
    "monoclinic", "orthorhombic", "rhombohedral", "tetragonal", "hexagonal",
    "cubic"]

"""Give the size of the point group of a Bravais lattice."""
pointgroup_size = Dict{String,Integer}(i=>j for (i,j)=zip(lat_types, num_ops))

end #module
