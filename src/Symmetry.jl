module Symmetry

import CDDLib: Library
import Polyhedra: HalfSpace, polyhedron, points
import LinearAlgebra: norm, det, I
import QHull: Chull, chull
import Combinatorics: permutations

include("Lattices.jl")
include("Utilities.jl")
import .Lattices: get_recip_latvecs, minkowski_reduce
import .Utilities: sample_circle, sample_sphere, contains, remove_duplicates

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
    if size(latvecs) == (2,2)
        pts=sample_circle(latvecs,radius,[0,0],rtol,atol)
    elseif size(latvecs) == (3,3)
        pts=sample_sphere(latvecs,radius,[0,0,0],rtol,atol)
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
    mapto_unitcell(pt,latvecs,inv_latvecs,coords,rtol,atol)

Map a point to the first unit (primitive) cell.

# Arguments
- `pt::AbstractArray{<:Real,1}`: a point in lattice or Cartesian coordinates.
- `latvecs::AbstractArray{<:Real,2}`: the basis vectors of the lattice as
    columns of an array.
- `inv_latvecs::AbstractArray{<:Real,2}`: the inverse of the matrix of that
    contains the lattice vectors.
- `coords::String`: indicates whether `pt` is in \"Cartesian\" or \"lattice\"
    coordinates.
- `rtol::Real=sqrt(eps(float(maximum(inv_latvecs))))`: a relative tolerance for
    floating point comparisons. Finite precision errors creep up when `pt` is
    transformed to lattice coordinates because the transformation requires
    calculating a matrix inverse. The components of the point in lattice
    coordinates are checked to ensure that values close to 1 are equal to 1.
- `atol::Real=0.0`: an absolute tolerance for floating point comparisons.

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
    coords::String,rtol::Real=sqrt(eps(float(maximum(inv_latvecs)))),
    atol::Real=0.0)::Array{Float64,1}
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
(Any[[0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0]], Any[[0.0 -1.0; -1.0 0.0], [0.0 -1.0; 1.0 0.0], [-1.0 0.0; 0.0 -1.0], [1.0 0.0; 0.0 -1.0], [-1.0 0.0; 0.0 1.0], [1.0 0.0; 0.0 1.0], [0.0 1.0; -1.0 0.0], [0.0 1.0; 1.0 0.0]])
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

    # Place atom positions in Cartesian coordinates.
    if coords == "lattice"
        atom_pos = real_latvecs*atom_pos
    end

    real_latvecs = minkowski_reduce(real_latvecs)
    inv_latvecs = inv(real_latvecs)
    pointgroup = calc_pointgroup(real_latvecs,rtol,atol)
    numatoms = length(atom_types)

    # Map points to the primitive cell.
    atom_pos = reduce(hcat,[mapto_unitcell(atom_pos[:,i],real_latvecs,
                inv_latvecs,"Cartesian",rtol) for i=1:numatoms])

    ops_spacegroup=[]
    trans_spacegroup=[]
    kindᵣ=atom_types[1]
    opts = []

    # Atoms of the same type as the first atom.
    same_atoms=findall(x->x==kindᵣ,atom_types)

    for op in pointgroup
        # Use the first atom to calculate allowed fractional translations.
        posᵣ=op*atom_pos[:,1]

        for atomᵢ=same_atoms
            # Calculate a candidate fractional translation.
            ftrans = mapto_unitcell(atom_pos[:,atomᵢ]-posᵣ,real_latvecs,
                inv_latvecs,"Cartesian",rtol,atol)

            if !contains(ftrans, opts)
                append!(opts, [ftrans])
            end

            # Check that all atoms land on the lattice when rotated and
            # translated.
            mapped = false
            for atomⱼ = 1:numatoms
                kindⱼ = atom_types[atomⱼ]
                posⱼ = mapto_unitcell(op*atom_pos[:,atomⱼ] + ftrans,
                    real_latvecs,inv_latvecs,"Cartesian",rtol,atol)

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
    calc_bz(real_latvecs, convention, bzformat)

Calculate the Brillouin zone for the given real-space lattice basis.

# Arguments
- `real_latvecs::AbstractArray{<:Real,2}`: the real-space lattice vectors or
    primitive translation vectors as columns of a 2x2 or 3x3 array.
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
    verts = reduce(hcat,points(polyhedron(bz,Library())))
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
        atol::Real=0.0)

    (prim_types,prim_pos,prim_latvecs) = make_primitive(real_latvecs,atom_types,
        atom_pos,coords)
    pointgroup = calc_spacegroup(prim_latvecs,prim_types,prim_pos,coords)[2]
	sizepg = size(pointgroup,1)
    bzformat = "half-space"
    bz = calc_bz(prim_latvecs,convention,bzformat)
    bz_vertices = collect(points(polyhedron(bz,Library())))
    ibz = bz
    for v=bz_vertices
		for i=length(pointgroup):-1:1
			op = pointgroup[i]
            vʳ=op*v
            if !isapprox(vʳ,v,rtol=rtol,atol=atol) #&& rpg[i] == 0
                a = vʳ-v
                ibz = ibz ∩ HalfSpace(a,0)
				deleteat!(pointgroup,i)
            end
		if length(pointgroup) == 0
			break
		end
        end
    end

    ibz_vertices = reduce(hcat,points(polyhedron(ibz,Library())))
    bz_vertices = reduce(hcat,bz_vertices)
    bzvolume = chull(Array(bz_vertices')).volume
    ibzvolume = chull(Array(ibz_vertices')).volume
    reduction = bzvolume/ibzvolume

    if !(ibzvolume ≈ bzvolume/sizepg)
		@show ibzvolume
		@show bzvolume/sizepg
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

@doc """
    make_primitive(real_latvecs,atom_types,atom_pos,coords,rtol,atol)

Make a given unit cell primitive.

This is a Julia translation of the function by the same in
    https://github.com/msg-byu/symlib.

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
- `prim_types::AbstractArray{<:Int,1}`: a list of atom types as integers in the
    primitive unit cell.
- `prim_pos::AbstractArray{<:Real,2}`: the positions of the atoms in in the 
    crystal structure as columns of an array.
- `prim_latvecs::AbstractArray{<:Real,2}`: the primitive lattice vectors as
    columns of an array.

# Examples
```jldoctest
import IBZ.Lattices: genlat_CUB
import IBZ.Symmetry: make_primitive
a = 1.0
real_latvecs = genlat_CUB(a)
atom_types = [0,0]
atom_pos = Array([0 0 0; 0.5 0.5 0.5]')
ibzformat = "convex hull"
coords = "Cartesian"
convention = "ordinary"
make_primitive(real_latvecs, atom_types, atom_pos, coords)
# output
([0], [0.0; 0.0; 0.0], [1.0 0.0 0.5; 0.0 1.0 0.5; 0.0 0.0 0.5])
```
"""
function make_primitive(real_latvecs::AbstractArray{<:Real,2},
    atom_types::AbstractArray{<:Int,1},atom_pos::AbstractArray{<:Real,2},
    coords::String,rtol::Real=sqrt(eps(float(maximum(real_latvecs)))),
    atol::Real=0.0)

    (ftrans, pg) = calc_spacegroup(real_latvecs,atom_types,atom_pos,coords)
    # Unique fractional translations
    ftrans = remove_duplicates(reduce(hcat, ftrans))

    # No need to consider the origin.
    dim = size(real_latvecs,1)
    if dim==2 origin=[0,0] else origin=[0,0,0] end

    for i = 1:size(ftrans,2)
        if isapprox(ftrans[:,i],origin)
            ftrans  = ftrans[:,1:end .!= i]
            break
        end
    end

    prim_latvecs = real_latvecs
    inv_latvecs = inv(real_latvecs)

    if size(ftrans,2)!=0
        # Number of lattice vector options
        nopts = size(ftrans,2)+3

        # Lattice vector options
        opts = hcat(real_latvecs, ftrans)
        if dim == 2
            for i=1:nopts-2, j=i+1:nopts-1
                prim_latvecs = opts[:,[i,j]]
                inv_latvecs = inv(prim_latvecs)
                # Move all translations into the potential primitive unit cell.
                test = [mapto_unitcell(ftrans[:,i],real_latvecs,
                            inv_latvecs,"Cartesian") for i=1:size(ftrans,2)]

                # Check if all coordinates of all fractional translations are
                # integers in lattice coordinates.
                test = [inv_latvecs*t for t=test]
                test = isapprox(mod.(collect(Iterators.flatten(test)),1)
                            ,zeros(2*length(test)))
                if test
                    break
                end
            end
        else
            for i=1:nopts-2, j=i+1:nopts-1, k=j+1:nopts
                prim_latvecs = opts[:,[i,j,k]]
                inv_latvecs = inv(prim_latvecs)
                # Move all translations into the potential primitive unit cell.
                test = [mapto_unitcell(ftrans[:,i],real_latvecs,
                            inv_latvecs,"Cartesian") for i=1:size(ftrans,2)]

                # Check if all coordinates of all fractional translations are
                # integers in lattice coordinates.
                test = [inv_latvecs*t for t=test]
                test = isapprox(mod.(collect(Iterators.flatten(test)),1)
                            ,zeros(3*length(test)))
                if test
                    break
                end
            end
        end
    end

    # Map all atoms into the primitive cell.
    all_prim_pos = atom_pos
    if coords == "lattice"
        all_prim_pos = real_latvecs*all_prim_pos
    end
    all_prim_pos = reduce(hcat,[mapto_unitcell(all_prim_pos[:,i], prim_latvecs,
        inv_latvecs,"Cartesian") for i=1:length(atom_types)])

    # Remove atoms at the same positions that are the same type.
    tmp = remove_duplicates(vcat(all_prim_pos,atom_types'))
    prim_types = Int.(tmp[dim + 1,:])
    prim_pos = tmp[1:dim,:]

    (prim_types,prim_pos,prim_latvecs)
end

end #module
