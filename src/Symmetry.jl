module Symmetry

import CDDLib
import Polyhedra: HalfSpace, polyhedron, points, Intersection, Library, DefaultLibrary
import LinearAlgebra: norm, det, I, dot, checksquare
import QHull: chull, Chull
import Combinatorics: permutations
import Base.Iterators: product, flatten

import ..Lattices: get_recip_latvecs, minkowski_reduce
import ..Utilities: sample_circle, sample_sphere, contains, unique_points, remove_duplicates

@doc """
    calc_pointgroup(latvecs;rtol,atol)

Calculate the point group of a lattice in 2D or 3D.

# Arguments
- `latvecs::AbstractMatrix{<:Real}`: the basis of the lattice as columns of a
    matrix.
- `rtol::Real=sqrt(eps(float(maximum(real_latvecs))))`: a relative tolerance for
    floating point comparisons. It is used to compare lengths of vectors and the
    volumes of primitive cells.
- `atol::Real=1e-9`: an absolute tolerance for floating point comparisons. It is
    used to compare lengths of vectors and the volumes of primitive cells.

# Returns
- `pointgroup::Vector{Matrix{Float64}}`: the point group of the lattice. The
    operators operate on points in Cartesian coordinates from the right.

# Examples
```jldoctest
using SymmetryReduceBZ
basis = [1 0; 0 1]
SymmetryReduceBZ.Symmetry.calc_pointgroup(basis)
# output
8-element Vector{Matrix{Float64}}:
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
function calc_pointgroup(latvecs::AbstractMatrix{<:Real};
    rtol::Real=sqrt(eps(float(maximum(latvecs)))),
    atol::Real=1e-9)

    dim=checksquare(latvecs)
    latvecs = minkowski_reduce(latvecs,rtol=rtol,atol=atol)
    radius = maximum([norm(latvecs[:,i]) for i=axes(latvecs,2)])*1.001
    if dim == 2
        pts=sample_circle(latvecs,radius,[0,0],rtol=rtol,atol=atol)
    elseif dim == 3
        pts=sample_sphere(latvecs,radius,[0,0,0],rtol=rtol,atol=atol)
    else
        throw(ArgumentError("The lattice basis must be a 2x2 or 3x3 matrix."))
    end

    normsᵢ=[norm(latvecs[:,i]) for i=axes(latvecs,2)]
    sizeᵢ=abs(det(latvecs))
    pointgroup=Array{Float64,2}[]
    inv_latvecs=inv(latvecs)
    for perm=permutations(1:size(pts,2),dim)
        latvecsᵣ=pts[:,perm]
        normsᵣ=[norm(latvecsᵣ[:,i]) for i=axes(latvecsᵣ,2)]
        sizeᵣ=abs(det(latvecsᵣ))
        # Point operation preserves length of lattice vectors and
        # volume of primitive cell.
        if (isapprox(normsᵢ,normsᵣ,rtol=rtol,atol=atol) &&
            isapprox(sizeᵢ,sizeᵣ,rtol=rtol,atol=atol))
            op=latvecsᵣ*inv_latvecs
            # Point operators are orthogonal.
            if isapprox(op'*op,I,rtol=rtol,atol=atol)
                push!(pointgroup,op)
            end
        end
    end
    pointgroup
end


@doc """
    mapto_unitcell(pt,latvecs,inv_latvecs,coordinates;rtol,atol)

Map a point to the first unit cell.

# Arguments
- `pt::AbstractVector{<:Real}`: a point in lattice or Cartesian coordinates.
- `latvecs::AbstractMatrix{<:Real}`: the basis vectors of the lattice as
    columns of a matrix.
- `inv_latvecs::AbstractMatrix{<:Real}`: the inverse of the matrix of that
    contains the lattice vectors.
- `coordinates::String`: indicates whether `pt` is in \"Cartesian\" or
    \"lattice\" coordinates.
- `rtol::Real=sqrt(eps(float(maximum(inv_latvecs))))`: a relative tolerance for
    floating point comparisons. Finite precision errors creep up when `pt` is
    transformed to lattice coordinates because the transformation requires
    calculating a matrix inverse. The components of the point in lattice
    coordinates are checked to ensure that values close to 1 are equal to 1.
- `atol::Real=1e-9`: an absolute tolerance for floating point comparisons.

# Returns
- `AbstractVector{<:Real}`: a translationally equivalent point to `pt` in the
    first unit cell in the same coordinates.

# Examples
```jldoctest
using SymmetryReduceBZ
real_latvecs = [0 1 2; 0 -1 1; 1 0 0]
inv_latvecs=inv(real_latvecs)
pt=[1,2,3.2]
coordinates = "Cartesian"
SymmetryReduceBZ.Symmetry.mapto_unitcell(pt,real_latvecs,inv_latvecs,
    coordinates)
# output
3-element Vector{Float64}:
 0.0
 0.0
 0.20000000000000018
```
"""
function mapto_unitcell(pt::AbstractVector{<:Real},
    latvecs::AbstractMatrix{<:Real},inv_latvecs::AbstractMatrix{<:Real},
    coordinates::String;rtol::Real=sqrt(eps(float(maximum(inv_latvecs)))),
    atol::Real=1e-9)
    T = promote_type(eltype(pt), eltype(latvecs), eltype(inv_latvecs))
    if coordinates == "Cartesian"
        latvecs*T[isapprox(mod(c,1),1,rtol=rtol,atol=atol) ? 0 : mod(c,1) for c=inv_latvecs*pt]
    elseif coordinates == "lattice"
        T[isapprox(mod(c,1),1,rtol=rtol,atol=atol) ? 0 : mod(c,1) for c=pt]
    else
        throw(ArgumentError("Allowed coordinates of the point are \"Cartesian\"
                and \"lattice\"."))
    end
end

@doc """
    mapto_unitcell(points,latvecs,inv_latvecs,coordinates;rtol,atol)

Map points as columns of a matrix to the unitcell.
"""
function mapto_unitcell(points::AbstractMatrix{<:Real},
    latvecs::AbstractMatrix{<:Real},inv_latvecs::AbstractMatrix{<:Real},
    coordinates::String;rtol::Real=sqrt(eps(float(maximum(inv_latvecs)))),
    atol::Real=1e-9)

    reduce(hcat,[mapto_unitcell(points[:,i],latvecs,inv_latvecs,coordinates,
        rtol=rtol,atol=atol) for i=axes(points,2)])
end

@doc """
    mapto_bz(kpoint,recip_latvecs,inv_latvecs,coordinates;rtol,atol)

Map a k-point to a translationally equivalent point within the Brillouin zone.

# Arguments
- `kpoint::AbstractVector{<:Real}`: a single *k*-point in lattice or Cartesian
    coordinates.
- `recip_latvecs::AbstractMatrix{<:Real}`: the reciprocal lattice vectors as
    columns of a matrix.
- `inv_latvecs::AbstractMatrix{<:Real}`: the inverse matrix of the reciprocal
    lattice vectors.
- `coordinates::String`: the coordinates of the given point, either \"lattice\"
    or \"Cartesian\". The point returned will be in the same coordinates.
- `rtol::Real=sqrt(eps(float(maximum(recip_latvecs))))`: a relative tolerance
    for floating point comparisons. Finite precision errors creep in when `pt`
    is transformed to lattice coordinates because the transformation requires
    calculating a matrix inverse. The components of the *k*-point in lattice
    coordinates are checked to ensure that values close to 1 are equal to 1.
- `atol::Real=1e-9`: an absolute tolerance for floating point comparisons.

# Returns
- `bz_point::AbstractVector{<:Real}`: the symmetrically equivalent *k*-point
    within the Brillouin zone in either lattice or Cartesian coordinates,
    depending on the coordinates specified.

# Examples
```jldoctest
import SymmetryReduceBZ.Symmetry: mapto_bz
import LinearAlgebra: inv
recip_latvecs = [1 0 0; 0 1 0; 0 0 1]
inv_latvecs = inv(recip_latvecs)
kpoint = [2, 3, 2]
coordinates = "Cartesian"
mapto_bz(kpoint, recip_latvecs, inv_latvecs, coordinates)
# output
3-element Vector{Float64}:
 0.0
 0.0
 0.0
```
"""
function mapto_bz(kpoint::AbstractVector{<:Real},
    recip_latvecs::AbstractMatrix{<:Real},
    inv_rlatvecs::AbstractMatrix{<:Real},coordinates::String;
    rtol::Real=sqrt(eps(float(maximum(recip_latvecs)))),
    atol::Real=1e-9)

    uc_point = mapto_unitcell(kpoint,recip_latvecs,inv_rlatvecs,coordinates,
        rtol=rtol,atol=atol)
    
    coordinates == "lattice" && (uc_point = recip_latvecs*uc_point)

    bz_dist = Inf
    bz_point = similar(uc_point)
    
    if checksquare(recip_latvecs) == 2
        shifts = collect.(product(0:-1:-1,0:-1:-1))
    else
        shifts = collect.(product(0:-1:-1,0:-1:-1,0:-1:-1))
    end

    for shift = shifts
        tpoint = uc_point + recip_latvecs*shift
        if norm(tpoint) < bz_dist
            bz_point = tpoint
            bz_dist = norm(tpoint)
        end
    end

    if coordinates == "Cartesian"
        bz_point
    elseif coordinates == "lattice"
        inv_rlatvecs*bz_point
    else
        throw(ArgumentError("Allowed coordinates of the point are \"Cartesian\"
        and \"lattice\"."))
    end
end

@doc """
    mapto_bz(kpoints,recip_latvecs,inv_rlatvecs,coordinates;rtol,atol)

Map points as columns of a matrix to the Brillouin zone.
"""
function mapto_bz(kpoints::AbstractMatrix{<:Real},
    recip_latvecs::AbstractMatrix{<:Real},
    inv_rlatvecs::AbstractMatrix{<:Real},coordinates::String;
    rtol::Real=sqrt(eps(float(maximum(recip_latvecs)))),
    atol::Real=1e-9)

    reduce(hcat,[mapto_bz(kpoints[:,i],recip_latvecs,inv_rlatvecs,coordinates,
        rtol=rtol, atol=atol) for i=axes(kpoints,2)])
end

@doc """
    inhull(point, chull; rtol, atol)

Check if a point lies within a convex hull (including the boundaries).

# Arguments
- `point::AbstractVector{<:Real}`: a point in Cartesian coordinates.
- `chull::Chull{Float64}`: a convex hull in 2D or 3D.
- `rtol::Real=sqrt(eps(float(maximum(flatten(chull.points)))))`: a relative
    tolerance for floating point comparisons. Needed when a point is on the
    boundary of the convex hull.
- `atol::Real=1e-9`: an absolute tolerance for floating point comparisons.

# Returns
- `inside::Bool`: if true, the point lies within the convex hull.

# Examples
```jldoctest
import QHull: chull
import SymmetryReduceBZ.Symmetry: inhull
pts = [0.5 0.25; 0.5 -0.25; -0.5 -0.25; -0.5 0.25]
pt = [0,0]
convexhull = chull(pts)
inhull(pt,convexhull)
# output
true
```
"""
function inhull(point::AbstractVector{<:Real}, chull::Chull{<:Real};
    rtol::Real=sqrt(eps(float(maximum(flatten(chull.points))))),
    atol::Real=1e-9)

    # hullpts = Array(chull.points')
    distances = chull.facets[:,end]
    norms = chull.facets[:,begin:end-1]'
    inside = true
    for (dist,i)=zip(distances, axes(norms,2))
        s = dot(point + dist*norms[:,i], norms[:,i])

        if !(s <= 0 || isapprox(s,0.0,rtol=rtol,atol=atol))
            inside = false
            break
        end
    end
    inside
end

"""
    mapto_ibz(kpoint,recip_latvecs,inv_rlatvecs,ibz,pointgroup,coordinates;rtol,
        atol)

Map a point to a symmetrically equivalent point within the IBZ.

# Arguments
- `kpoint::AbstractVector{<:Real}`: a *k*-point in 2D or 3D in Cartesian
    coordinates.
- `recip_latvecs::AbstractMatrix{<:Real}`: the reciprocal lattice vectors as
    columns of a a matrix.
- `inv_rlatvecs::AbstractMatrix{<:Real}`: the inverse of the square matrix
    `recip_latvecs`.
- `ibz::Chull{Float64}`: the irreducible Brillouin zone as as a convex hull
    objects from `QHull`.
- `pointgroup::Vector{Matrix{Float64}}`: a list of point symmetry operators
    in matrix form that operate on points from the left.
- `coordinates::String`: the coordinates the *k*-point is in. Options are
    \"lattice\" and \"Cartesian\". The *k*-point within the IBZ is returned in
    the same coordinates.
- `rtol::Real=sqrt(eps(float(maximum(recip_latvecs))))`: a relative tolerance
    for floating point comparisons. The *k*-point is first mapped the unit cell
    and `rtol` is used when comparing components of the *k*-point to 1. It is
    also used for comparing floats to zero when checking if the point lies
    within `ibz`.
- `atol::Real=1e-9`: an absolute tolerance for floating point comparisons. This
    is used everywhere `rtol` is used.

# Returns
- `rot_point::AbstractVector{<:Real}`: a symmetrically equivalent *k*-point to
    `kpoint` within the irreducible Brillouin zone in the same coordinates as
    `coordinates`.

# Examples
```jldoctest
import SymmetryReduceBZ.Lattices: get_recip_latvecs
import SymmetryReduceBZ.Symmetry: calc_spacegroup, mapto_ibz
import LinearAlgebra: inv
import QHull: chull
real_latvecs = [1 0; 0 2]
atom_types=[0]
atom_pos=Array([0 0]')
coordinates="Cartesian"
convention="ordinary"
recip_latvecs = get_recip_latvecs(real_latvecs,convention)
inv_rlatvecs = inv(recip_latvecs)
(ftrans,pg) = calc_spacegroup(real_latvecs,atom_types,atom_pos,coordinates)
ibz = chull([0.0 0.25; 0.0 0.0; 0.5 0.0; 0.5 0.25])
kpoint = [2,3]
ibz_point = mapto_ibz(kpoint,recip_latvecs,inv_rlatvecs,ibz,pg,coordinates)
# output
2-element Vector{Float64}:
 0.0
 0.0
```
"""
function mapto_ibz(kpoint::AbstractVector{<:Real},
        recip_latvecs::AbstractMatrix{<:Real},
        inv_rlatvecs::AbstractMatrix{<:Real}, ibz::Chull{Float64},
        pointgroup::Vector{Matrix{Float64}}, coordinates::String;
        rtol::Real=sqrt(eps(float(maximum(recip_latvecs)))),
        atol::Real=1e-9)

    bzpoint = mapto_bz(kpoint,recip_latvecs,inv_rlatvecs,coordinates,
        rtol=rtol,atol=atol)
    
    coordinates == "lattice" && (bzpoint = recip_latvecs * bzpoint)

    for op in pointgroup
        rot_point = op*bzpoint
        if inhull(rot_point,ibz,rtol=rtol,atol=atol)
            if coordinates == "lattice"
                rot_point = inv_rlatvecs*rot_point
            end
            return rot_point
        end
    end
    error("Failed to map the k-point to the IBZ.")
end

"""
    mapto_ibz(kpoints,recip_latvecs,inv_rlatvecs,ibz,pointgroup,coordinates;
        rtol,atol)

Map points as columns of a matrix to the IBZ and then remove duplicate points.
"""
function mapto_ibz(kpoints::AbstractMatrix{<:Real},
        recip_latvecs::AbstractMatrix{<:Real},
        inv_rlatvecs::AbstractMatrix{<:Real}, ibz::Chull{Float64},
        pointgroup::Vector{Matrix{Float64}}, coordinates::String;
        rtol::Real=sqrt(eps(float(maximum(recip_latvecs)))),
        atol::Real=1e-9)

    ibzpts = reduce(hcat,[mapto_ibz(kpoints[:,i],recip_latvecs,
        inv_rlatvecs,ibz,pointgroup,coordinates) for i=axes(kpoints,2)])

    unique_points(ibzpts,rtol=rtol,atol=atol)
end

@doc """
    calc_spacegroup(real_latvecs,atom_types,atom_pos,coordinates;rtol,atol)

Calculate the space group of a crystal structure.

# Arguments
- `real_latvecs::AbstractMatrix{<:Real}`: the basis of the lattice as columns
    of a matrix.
- `atom_types::AbstractVector{<:Int}`: a list of atom types as integers.
- `atom_pos::AbstractMatrix{<:Real}`: the positions of atoms in the crystal
    structure as columns of a matrix.
- `coordinates::String`: indicates the positions of the atoms are in \"lattice\"
    or \"Cartesian\" coordinates.
- `rtol::Real=sqrt(eps(float(maximum(real_latvecs))))` a relative tolerance for
    floating point comparisons.
- `atol::Real=1e-9`: an absolute tolerance for floating point comparisons.

# Returns
- `spacegroup::Tuple`: the space group of the crystal structure. The first element of
    `spacegroup` is a list of fractional translations, and the second element is
    a list of point operators. The translations are in Cartesian coordinates,
    and the operators operate on points in Cartesian coordinates.

# Examples
```jldoctest
using SymmetryReduceBZ
real_latvecs = Array([1 0; 2 1]')
atom_types = [0, 1]
atom_pos = Array([0 0; 0.5 0.5]')
coordinates = "Cartesian"
SymmetryReduceBZ.Symmetry.calc_spacegroup(real_latvecs,atom_types,atom_pos,
    coordinates)
# output
([[0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0]], [[0.0 -1.0; -1.0 0.0], [0.0 -1.0; 1.0 0.0], [-1.0 0.0; 0.0 -1.0], [1.0 0.0; 0.0 -1.0], [-1.0 0.0; 0.0 1.0], [1.0 0.0; 0.0 1.0], [0.0 1.0; -1.0 0.0], [0.0 1.0; 1.0 0.0]])
```
"""
function calc_spacegroup(real_latvecs::AbstractMatrix{<:Real},
    atom_types::AbstractVector{<:Int},atom_pos::AbstractMatrix{<:Real},
    coordinates::String;rtol::Real=sqrt(eps(float(maximum(real_latvecs)))),
    atol::Real=1e-9)

    if length(atom_types) != size(atom_pos,2)
        throw(ArgumentError("The number of atom types and positions must be the
            same."))
    end

    # Place atom positions in Cartesian coordinates.
    atom_pos_1 = coordinates == "lattice" ? real_latvecs*atom_pos : one(real_latvecs)*atom_pos

    real_latvecs = minkowski_reduce(real_latvecs,rtol=rtol,atol=atol)
    inv_latvecs = inv(real_latvecs)
    pointgroup = calc_pointgroup(real_latvecs,rtol=rtol,atol=atol)
    numatoms = length(atom_types)

    # Map points to the primitive cell.
    atom_pos = mapto_unitcell(atom_pos_1[:,1:numatoms],real_latvecs,inv_latvecs,"Cartesian",rtol=rtol,atol=atol)

    ops_spacegroup=Array{Float64,2}[]
    trans_spacegroup=Array{Float64,1}[]
    kindᵣ=atom_types[1]
    # opts = Array{Float64,1}[]

    # Atoms of the same type as the first atom.
    same_atoms=findall(x->x==kindᵣ,atom_types)

    for op in pointgroup
        # Use the first atom to calculate allowed fractional translations.
        posᵣ=op*atom_pos[:,1]

        for atomᵢ=same_atoms
            # Calculate a candidate fractional translation.
            ftrans = mapto_unitcell(atom_pos[:,atomᵢ]-posᵣ,real_latvecs,
                inv_latvecs,"Cartesian",rtol=rtol,atol=atol)

            # if !contains(ftrans, opts)
            #     push!(opts, ftrans)
            # end

            # Check that all atoms land on the lattice when rotated and
            # translated.
            mapped = false
            for atomⱼ = 1:numatoms
                kindⱼ = atom_types[atomⱼ]
                posⱼ = mapto_unitcell(op*atom_pos[:,atomⱼ] + ftrans,
                    real_latvecs,inv_latvecs,"Cartesian",rtol=rtol,atol=atol)

                mapped = any([kindⱼ == atom_types[i] &&
                        isapprox(posⱼ,atom_pos[:,i],rtol=rtol,atol=atol) ?
                            true : false for i=1:numatoms])
                if !mapped continue end
            end
            if mapped
                push!(ops_spacegroup,op)
                push!(trans_spacegroup,ftrans)
            end
        end
    end

    (trans_spacegroup,ops_spacegroup)
end

@doc """
    calc_bz(real_latvecs,atom_types,atom_pos,coordinates,bzformat,makeprim,
        convention,library;rtol,atol)

Calculate the Brillouin zone for the given real-space lattice basis.

# Arguments
- `real_latvecs::AbstractMatrix{<:Real}`: the real-space lattice vectors or
    primitive translation vectors as columns of a 2x2 or 3x3 matrix.
- `atom_typesAbstractVector{<:Int}`: a list of atom types as integers.
- `atom_pos::AbstractMatrix{<:Real}`: the positions of atoms in the crystal
    structure as columns of a matrix.
- `coordinates::String`: indicates the positions of the atoms are in \"lattice\"
    or \"Cartesian\" coordinates.
- `bzformat::String`: the format of the Brillouin zone. Options include
    \"convex hull\" and \"half-space\".
- `makeprim::Bool=false`: make the unit cell primitive before calculating the
    the BZ if equal to `true`.
- `convention::String="ordinary"`: the convention used to go between real and
    reciprocal space. The two conventions are ordinary (temporal) frequency and
    angular frequency. The transformation from real to reciprocal space is
    unitary if the convention is ordinary.
- `library::Polyhedra.Library=CDDLib.Library()`: a polyhedron manipulation library
- `rtol::Real=sqrt(eps(float(maximum(real_latvecs))))` a relative tolerance for
    floating point comparisons.
- `atol::Real=1e-9`: an absolute tolerance for floating point comparisons.

# Returns
- `bz`: the vertices or half-space representation of the Brillouin zone
    depending on the value of `vertsOrHrep`.

# Examples
```jldoctest
using SymmetryReduceBZ
real_latvecs = [1 0; 0 1]
convention="ordinary"
atom_types=[0]
atom_pos = Array([0 0]')
coordinates = "Cartesian"
bzformat = "convex hull"
makeprim=false
SymmetryReduceBZ.Symmetry.calc_bz(real_latvecs,atom_types,atom_pos,coordinates,
    bzformat,makeprim,convention)
# output
Convex Hull of 4 points in 2 dimensions
Hull segment vertex indices:
Int32[3, 2, 1, 4]
Points on convex hull in original order:

[0.5 0.5; 0.5 -0.5; -0.5 -0.5; -0.5 0.5]
```
"""
function calc_bz(real_latvecs::AbstractMatrix{<:Real},
    atom_types::AbstractVector{<:Int},atom_pos::AbstractMatrix{<:Real},
    coordinates::String,bzformat::String,makeprim::Bool=false,
    convention::String="ordinary", library::Library=CDDLib.Library();
    rtol::Real=sqrt(eps(float(maximum(real_latvecs)))),atol::Real=1e-9)

    if makeprim
        (prim_latvecs,prim_types,prim_pos) = make_primitive(real_latvecs,
            atom_types,atom_pos,coordinates,rtol=rtol,atol=atol)
    else
        (prim_types,prim_pos,prim_latvecs)=(atom_types,float(atom_pos),float(real_latvecs))
    end

    recip_latvecs = minkowski_reduce(get_recip_latvecs(prim_latvecs,convention),
        rtol=rtol,atol=atol)
    latpts_ = Vector{eltype(recip_latvecs)}[]
    if size(real_latvecs) == (2,2)
        append!(latpts_,[recip_latvecs*[i,j] for (i,j)=product(-2:2,-2:2)])
    elseif size(real_latvecs) == (3,3)
        append!(latpts_,[recip_latvecs*[i,j,k] for (i,j,k)=product(-2:2,-2:2,-2:2)])
    else
        throw(ArgumentError("The lattice vectors must be a 2x2 or 3x3 matrix."))
    end
    latpts = reduce(hcat, latpts_)

    distances = [norm(latpts[:,i]) for i=1:size(latpts,2)] .^2 ./2
    bz = HalfSpace(latpts[:,2],distances[2]) ∩ HalfSpace(latpts[:,3],distances[3])
    for i=4:size(latpts,2)
        bz = bz ∩ HalfSpace(latpts[:,i],distances[i])
    end
    verts = reduce(hcat,collect(points(polyhedron(bz, library))))
    bzvolume = chull(Array(verts')).volume

    if !(bzvolume ≈ abs(det(recip_latvecs)))
        error("The area or volume of the Brillouin zone is incorrect.")
    end

    if bzformat == "half-space"
        bz
    elseif bzformat == "convex hull"
        chull(Array(verts'))    # chull is a type-unstable function
    else
        throw(ArgumentError("Formats for the BZ include \"half-space\" and
            \"convex hull\"."))
    end
end

@doc """
    calc_ibz(real_latvecs,atom_types,atom_pos,coordinates,ibzformat,makeprim,
        convention,library;rtol,atol)

Calculate the irreducible Brillouin zone of a crystal structure in 2D or 3D.

# Arguments
- `real_latvecs::AbstractMatrix{<:Real}`: the basis of a real-space lattice as
    columns of a matrix.
- `atom_types:AbstractVector{<:Int}`: a list of atom types as integers.
- `atom_pos::AbstractMatrix{<:Real}`: the positions of atoms in the crystal
    structure as columns of a matrix.
- `coordinates::String`: indicates the positions of the atoms are in \"lattice\"
    or \"Cartesian\" coordinates.
- `ibzformat::String`: the format of the irreducible Brillouin zone. Options
    include \"convex hull\" and \"half-space\".
- `makeprim::Bool=false`: make the unit cell primitive before calculating the
    IBZ if true.
- `convention::String="ordinary"`: the convention used to go between real and
    reciprocal space. The two conventions are ordinary (temporal) frequency and
    angular frequency. The transformation from real to reciprocal space is
    unitary if the convention is ordinary.
- `library::Polyhedra.Library=CDDLib.Library()`: a polyhedron manipulation library
- `rtol::Real=sqrt(eps(float(maximum(real_latvecs))))` a relative tolerance for
    floating point comparisons.
- `atol::Real=1e-9`: an absolute tolerance for floating point comparisons.

# Returns
- `ibz`: the irreducible Brillouin zone as a convex hull or intersection of
    half-spaces.

# Examples
```jldoctest
using SymmetryReduceBZ
real_latvecs = [1 0; 0 1]
convention="ordinary"
atom_types=[0]
atom_pos = Array([0 0]')
coordinates = "Cartesian"
ibzformat = "convex hull"
makeprim=false
SymmetryReduceBZ.Symmetry.calc_ibz(real_latvecs,atom_types,atom_pos,coordinates,
    ibzformat,makeprim,convention)
# output
Convex Hull of 3 points in 2 dimensions
Hull segment vertex indices:
Int32[1, 2, 3]
Points on convex hull in original order:

[0.0 0.0; 0.5 0.0; 0.5 0.5]
```
"""
function calc_ibz(real_latvecs::AbstractMatrix{<:Real},
    atom_types::AbstractVector{<:Int},atom_pos::AbstractMatrix{<:Real},
    coordinates::String,ibzformat::String,makeprim::Bool=false,
    convention::String="ordinary", library::Library=CDDLib.Library();
    rtol::Real=sqrt(eps(float(maximum(real_latvecs)))),atol::Real=1e-9)

    if makeprim
        (prim_latvecs,prim_types,prim_pos) = make_primitive(real_latvecs,
            atom_types,atom_pos,coordinates,rtol=rtol,atol=atol)
    else
        (prim_types,prim_pos,prim_latvecs)=(atom_types,float(atom_pos),float(real_latvecs))
        if coordinates == "lattice"
            prim_pos = real_latvecs*prim_pos
        end
    end

    pointgroup = remove_duplicates(calc_spacegroup(prim_latvecs,prim_types,
        prim_pos,"Cartesian",rtol=rtol,atol=atol)[2],rtol=rtol,atol=atol)
    sizepg = size(pointgroup,1)
    # enforce type here due to instability in calc_bz
    bz::Intersection{eltype(prim_latvecs), Vector{eltype(prim_latvecs)}, Int64} =
        calc_bz(prim_latvecs,prim_types,prim_pos,"Cartesian","half-space",false,
        convention,rtol=rtol,atol=atol)
    bz_vertices = collect(points(polyhedron(bz, library)))
    ibz = bz
    for v=bz_vertices
        for i=length(pointgroup):-1:1
            op = pointgroup[i]
            vʳ=op*v
            if !isapprox(vʳ,v,rtol=rtol,atol=atol)
                a = vʳ-v
                ibz = ibz ∩ HalfSpace(a,0)
                deleteat!(pointgroup,i)
            end
        if length(pointgroup) == 0
            break
        end
        end
    end

    ibz_vertices = reduce(hcat,collect(points(polyhedron(ibz, library))))
    bz_vertices = reduce(hcat,bz_vertices)
    bzvolume = chull(Array(bz_vertices')).volume
    ibzvolume = chull(Array(ibz_vertices')).volume
    # reduction = bzvolume/ibzvolume

    if !(ibzvolume ≈ bzvolume/sizepg)
        @show ibzvolume
        @show bzvolume/sizepg
        error("The area or volume of the irreducible Brillouin zone is
            incorrect.")
    end
    if ibzformat == "half-space"
        ibz
    elseif ibzformat == "convex hull"
        chull(Array(ibz_vertices')) # chull is type-unstable
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
    make_primitive(real_latvecs,atom_types,atom_pos,coordinates;rtol,atol)

Make a given unit cell primitive.

This is a Julia translation of the function by the same in
    https://github.com/msg-byu/symlib.

# Arguments
- `real_latvecs::AbstractMatrix{<:Real}`: the basis of the lattice as columns
    of a matrix.
- `atom_types::AbstractVector{<:Int}`: a list of atom types as integers.
- `atom_pos::AbstractMatrix{<:Real}`: the positions of atoms in the crystal
    structure as columns of a matrix.
- `coordinates::String`: indicates the positions of the atoms are in \"lattice\"
    or \"Cartesian\" coordinates.
- `rtol::Real=sqrt(eps(float(maximum(real_latvecs))))` a relative tolerance for
    floating point comparisons.
- `atol::Real=1e-9`: an absolute tolerance for floating point comparisons.

# Returns
- `prim_latvecs::AbstractMatrix{<:Real}`: the primitive lattice vectors as
    columns of a matrix.
- `prim_types::AbstractVector{<:Int}`: a list of atom types as integers in the
    primitive unit cell.
- `prim_pos::AbstractMatrix{<:Real}`: the positions of the atoms in in the
    crystal structure as columns of a matrix in Cartesian coordinates.

# Examples
```jldoctest
import SymmetryReduceBZ.Lattices: genlat_CUB
import SymmetryReduceBZ.Symmetry: make_primitive
a = 1.0
real_latvecs = genlat_CUB(a)
atom_types = [0,0]
atom_pos = Array([0 0 0; 0.5 0.5 0.5]')
ibzformat = "convex hull"
coordinates = "Cartesian"
convention = "ordinary"
make_primitive(real_latvecs, atom_types, atom_pos, coordinates)
# output
([1.0 0.0 0.5; 0.0 1.0 0.5; 0.0 0.0 0.5], [0], [0.0; 0.0; 0.0;;])
```
"""
function make_primitive(real_latvecs::AbstractMatrix{<:Real},
    atom_types::AbstractVector{<:Integer},atom_pos::AbstractMatrix{<:Real},
    coordinates::String;rtol::Real=sqrt(eps(float(maximum(real_latvecs)))),
    atol::Real=1e-9)

    if length(atom_types) == 1
        return (float(real_latvecs), atom_types, float(atom_pos)) # fractional translations will be floating-point
    end

    (ftrans_, pg) = calc_spacegroup(real_latvecs,atom_types,atom_pos,coordinates,
        rtol=rtol,atol=atol)
    # Unique fractional translations
    ftrans = unique_points(reduce(hcat, ftrans_),rtol=rtol,atol=atol)

    # No need to consider the origin.
    dim = size(real_latvecs,1)
    if dim==2
        origin=[0.,0.]
    elseif dim==3
        origin=[0.,0.,0.]
    else
        throw(ArgumentError("Can only make 2D or 3D unit cells primitive."))
    end
    for i = axes(ftrans,2)
        if isapprox(ftrans[:,i],origin,rtol=rtol,atol=atol)
            ftrans  = ftrans[:,1:end .!= i]
            break
        end
    end

    prim_latvecs = similar(real_latvecs, promote_type(eltype(real_latvecs), eltype(ftrans)))
    inv_latvecs = inv(real_latvecs)
    mapped = false
    if size(ftrans,2)!=0
        # Number of lattice vector options
        nopts = size(ftrans,2)+dim
        # Lattice vector options
        opts = hcat(real_latvecs, ftrans)
        if dim == 2
            iters = [[i,j] for i=1:nopts-2 for j=i+1:nopts-1]
        else
            iters = [[i,j,k] for i=1:nopts-2 for j=i+1:nopts-1 for k=j+1:nopts]
        end

        for iter=iters
            prim_latvecs .= opts[:,iter]
            if isapprox(det(prim_latvecs),0,atol=atol)
                continue
            end
            inv_latvecs .= inv(prim_latvecs)
            # Move all lattice points into the potential primitive unit cell.
            tests = mapto_unitcell(opts[:,1:nopts],prim_latvecs,
                    inv_latvecs,"Cartesian",rtol=rtol,atol=atol)
            # Check if all coordinates of all fractional translations are
            # integers in lattice coordinates.
            test = [inv_latvecs*tests[:,i] for i=axes(tests,2)]
            if isapprox(mod.(collect(flatten(test)),1)
                ,zeros(dim*length(test)),rtol=rtol,atol=atol)
                mapped = true
                break
            end
        end
    end

    if !mapped
        prim_latvecs .= real_latvecs
        inv_latvecs .= inv(real_latvecs)
    end

    # Map all atoms into the primitive cell.
    all_prim_pos_ = similar(atom_pos, promote_type(eltype(real_latvecs), eltype(atom_pos)))
    all_prim_pos_ .= atom_pos
    if coordinates == "lattice"
        all_prim_pos_ .= real_latvecs*all_prim_pos_
    end
    all_prim_pos = mapto_unitcell(all_prim_pos_[:,1:length(atom_types)], prim_latvecs,
        inv_latvecs,"Cartesian",rtol=rtol,atol=atol)

    # Remove atoms at the same positions that are the same type.
    tmp = unique_points(vcat(all_prim_pos,atom_types'),rtol=rtol,atol=atol)
    prim_types = Int.(tmp[dim + 1,:])
    prim_pos = tmp[1:dim,:]

    (prim_latvecs,prim_types,prim_pos)
end

@doc """
    complete_orbit(pt,pointgroup,rtol=sqrt(eps(float(maximum(pt)))),atol=1e-9)

Complete the orbit of a point.

# Arguments
- `pt::AbstractVector{<:Real}`: the Cartesian coordinates of a point.
- `pointgroup::Vector{Matrix{Float64}}`: the point group operators
    in a nested list. The operators operate on points in Cartesian coordinates
    from the right.
- `rtol::Real=sqrt(eps(float(maximum(pt))))`: a relative tolerance.
- `atol::Real=1e-9`: an absolute tolerance.

# Returns
- `::AbstractMatrix{<:Real}`: the points of the orbit in Cartesian coordinates as
    columns of a matrix.

# Examples
```jldoctest
import SymmetryReduceBZ.Symmetry: complete_orbit
pt = [0.05, 0.0]
pointgroup = [[0.0 -1.0; -1.0 0.0], [0.0 -1.0; 1.0 0.0], [-1.0 0.0; 0.0 -1.0], [1.0 0.0; 0.0 -1.0], [-1.0 0.0; 0.0 1.0], [1.0 0.0; 0.0 1.0], [0.0 1.0; -1.0 0.0], [0.0 1.0; 1.0 0.0]]
complete_orbit(pt,pointgroup)
# output
2×4 Matrix{Float64}:
  0.0   0.0   -0.05  0.05
 -0.05  0.05   0.0   0.0
```
"""
function complete_orbit(pt::AbstractVector{<:Real},
        pointgroup::Vector{Matrix{Float64}};
        rtol::Real=sqrt(eps(float(maximum(pt)))),
        atol::Real=1e-9)
    unique_points(reduce(hcat,map(op->op*pt,pointgroup)),rtol=rtol,atol=atol)
end

@doc """
    complete_orbit(pt,pointgroup,rtol=sqrt(eps(float(maximum(pt)))),atol=1e-9)

Complete the orbits of multiple points.

# Arguments
- `pt::AbstractMatrix{<:Real}`: the Cartesian coordinates of a points as columns of
    a matrix.
- `pointgroup::Vector{Matrix{Float64}}`: the point group operators
    in a nested list. The operators operate on points in Cartesian coordinates from the right.
- `rtol::Real=sqrt(eps(float(maximum(pt))))`: a relative tolerance.
- `atol::Real=1e-9`: an absolute tolerance.

# Returns
- `::AbstractMatrix{<:Real}`: the unique points of the orbits in Cartesian coordinates as
    columns of a matrix.

# Examples
```jldoctest
import SymmetryReduceBZ.Symmetry: complete_orbit
pts = [0.0 0.05 0.1; 0.0 0.0 0.0]
pointgroup = [[0.0 -1.0; -1.0 0.0], [0.0 -1.0; 1.0 0.0], [-1.0 0.0; 0.0 -1.0], [1.0 0.0; 0.0 -1.0], [-1.0 0.0; 0.0 1.0], [1.0 0.0; 0.0 1.0], [0.0 1.0; -1.0 0.0], [0.0 1.0; 1.0 0.0]]
complete_orbit(pts,pointgroup)
# output
2×9 Matrix{Float64}:
 0.0   0.0   0.0   -0.05  0.05   0.0  0.0  -0.1  0.1
 0.0  -0.05  0.05   0.0   0.0   -0.1  0.1   0.0  0.0
```
"""
function complete_orbit(pts::AbstractMatrix{<:Real},
    pointgroup::Vector{Matrix{Float64}};
    rtol::Real=sqrt(eps(float(maximum(pts)))),
    atol::Real=1e-9)
    unique_points(reduce(hcat,reduce(hcat,[map(op->op*pts[:,i],pointgroup) for i=axes(pts,2)])))
 end
end #module
