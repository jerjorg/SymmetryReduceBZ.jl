module Lattices

using LinearAlgebra, Combinatorics

export get_recip_latvecs, minkowski_reduce, check_reduced, get_latparams
# 2D Bravais lattices
export genlat_SQR, genlat_HXG, genlat_REC, genlat_RECI, genlat_OBL
# 3D Bravais lattices
export genlat_CUB, genlat_FCC, genlat_BCC, genlat_TET, genlat_BCT, genlat_ORC,
    genlat_ORCF, genlat_ORCI, genlat_ORCC, genlat_HEX, genlat_RHL, genlat_MCL,
    genlat_MCLC, genlat_TRI

@doc """
    get_recip_latvecs(real_latvecs, convention)

Calculate the reciprocal lattice vectors.

# Arguments
- `real_latvecs::AbstractArray{<:Real,2}`: the real-space lattice vectors or
    primitive translation vectors as columns of a 2x2 or 3x3 array.
- `convention::String="ordinary"`: the convention used to go between real and
    reciprocal space. The two conventions are ordinary (temporal) frequency and
    angular frequency. The transformation from real to reciprocal space is
    unitary if the convention is ordinary.

# Returns
- `recip_latvecs::Array{<:Real,2}` the reciprocal lattice vectors (reciprocal
    primitive translation vectors) as columns of a 2x2 or 3x3 array.

# Examples
```jldoctest
using IBZ
real_latvecs=[1 0 0; 0 1 0; 0 0 1]
convention="angular"
IBZ.Lattices.get_recip_latvecs(real_latvecs,convention)
# output
3×3 Array{Float64,2}:
 6.28319  0.0      0.0
 0.0      6.28319  0.0
 0.0      0.0      6.28319
```
"""
function get_recip_latvecs(real_latvecs::AbstractArray{<:Real,2},
        convention::String="ordinary")::Array{Float64,2}
    if convention == "ordinary"
        recip_latvecs = Array(inv(real_latvecs)')
    elseif convention == "angular"
        recip_latvecs = Array(2π*inv(real_latvecs)')
    else
        throw(ArgumentError("The allowed conventions are \"ordinary\" and
        \"angular\"."))
    end
    recip_latvecs
end

@doc """
    get_latparams(latvecs)

Calculate the lattice constants and angles of a lattice basis.

# Arguments
- `latvecs::AbstractArray{<:Real,2}`: the lattice basis as columns of an array.

# Returns
- A list where the first element is a list lattice constants `(a,b,c)` and second
    lattice angles in radians `(α,β,γ)`.

# Examples
```jldoctest
using IBZ
latvecs = [1 0; 0 1]
IBZ.Lattices.get_latparams(latvecs)
# output
2-element Array{Array{Float64,1},1}:
 [1.0, 1.0]
 [1.5707963267948966, 1.5707963267948966]
```
"""
function get_latparams(latvecs)

    if size(latvecs) == (2,2)
        (a,b)=[latvecs[:,i] for i=1:2]
        θ=acos(dot(a,b)/(norm(a)*norm(b)))
        [[norm(a),norm(b)],[θ,θ]]
    elseif size(latvecs) == (3,3)
        (a,b,c)=[latvecs[:,i] for i=1:3]
        α=acos(dot(b,c)/(norm(b)*norm(c)))
        β=acos(dot(a,c)/(norm(a)*norm(c)))
        γ=acos(dot(a,b)/(norm(a)*norm(b)))
        [[norm(a),norm(b),norm(c)],[α,β,γ]]
    else
        throw(ArgumentError("The lattice vectors must be a 2x2 or 3x3 array."))
    end
end

"""
    reduce_basis!(basis,k)

Reduces the `k`th lattice vector. This is accomplished by locating the
lattice point closest to the projection of the `k`th lattice vector onto
the line or plane given by the other lattice vector(s), subtracting the
closest lattice point from the `k`th lattice vector, and reordering the
lattice vectors by increasing Euclidean norms.

# Arguments
- `basis::AbstractArray{<:Real,2}`: the lattice basis as columns of an array.
- `k::Int`: Keeps track of which lattice vector needs to be reduced.

# Returns
- `basis::AbstractArray{<:Real,2}`: the partially reduced lattice basis as
    columns of an array.
- `k::Int`: The index of the lattice vector that needs to be reduced next.

# Examples
```jldoctest
using IBZ
basis = Array([1 2 0; 0 1 0; 3 2 1]')
k=2
IBZ.Lattices.reduce_basis!(basis,k)
basis
# output
3×3 Array{Int64,2}:
 0  1  3
 1  2  2
 0  0  1
```
"""
function reduce_basis!(basis::AbstractArray{<:Real,2},k::Int)::Int
    if k == 2
        v1,v2=[basis[:,i] for i=1:k]
        i = round(Int64,dot(v1,v2)/dot(v1,v1))
        vecs = [v2-j*v1 for j=i-1:i+1]
        v2 = vecs[sortperm(norm.(vecs))[1]]
        if norm(v1) <= norm(v2)
            basis[:,k] = v2
            k=3
        else
            k=2
            basis[:,2] = v1
            basis[:,1] = v2
        end

    elseif k==3
        v1,v2,v3=[basis[:,i] for i=1:k]

        i = round(Int64,dot(v3,v1)/dot(v1,v1))
        j = round(Int64,dot(v3,v2)/dot(v2,v2))

        vecs = [[v3-m*v2-l*v1 for m=j-1:j+1,l=i-1:i+1]...]
        v3 = vecs[sortperm(norm.(vecs))[1]]
        if norm(v2) <= norm(v3)
            basis[:,k] = v3
            k=4
        else
            if norm(v1) <= norm(v3)
                k=3
                basis[:,2] = v3
                basis[:,3] = v2
            else
                k=2
                basis[:,1] = v3
                basis[:,2] = v1
                basis[:,3] = v2
            end
        end
    else
        throw(ArgumentError("Argument `k` can be 2 or 3."))
    end
    k
end

@doc """
        minkowski_reduce(basis)

Minkowski reduce a lattice basis. Follows the logic of Fig. 4 in
\"Low-Dimensional Lattice Basis Reduction Revisited\" by Nguyen, 2009.

# Arguments
- `basis::AbstractArray{<:Real,2}`: the lattice basis given by the columns
    of a 2x2 or 3x3 array.

# Returns
- `rbasis:: the Minkowski reduced lattice basis as columns of an array.

# Examples
```jldoctest
using IBZ
basis = [1 2 0; 0 1 0; 0 0 1]
IBZ.Lattices.minkowski_reduce(basis)
# output
3×3 Array{Int64,2}:
 0  1  0
 0  0  1
 1  0  0
```
"""
function minkowski_reduce(
    basis::AbstractArray{<:Real,2})::AbstractArray{<:Real,2}

    # Sort the lattice vectors by increasing norm.
    order = sortperm([basis[:,i] for i=1:size(basis,1)])
    rbasis = basis[:,order]

    k=2
    while k <= size(rbasis,1)
        k=reduce_basis!(rbasis,k)
    end
    rbasis
end

@doc """
    check_reduced(basis)

Verify a lattice basis is Minkowski reduced

# Arguments
- `basis::AbstractArray{<:Real,2}`: the lattice basis given by the columns
    of a 2x2 or 3x3 matrix.

# Returns
- `Bool`: a boolean that indicates if the lattice basis is reduced.

# Examples
``` jldoctest
using IBZ
basis = [1 0; 0 1]
IBZ.Lattices.check_reduced(basis)
# output
true
```
"""
function check_reduced(basis::AbstractArray{<:Real,2})::Bool
    if size(basis,2) == 2
        (a,b) = [basis[:,i] for i=1:2]
        all([norm(a) <= norm(b),
            norm(b) <= norm(b+a),
            norm(b) <= norm(b-a)])
    elseif size(basis,1) == 3
        (a,b,c) = [basis[:,i] for i=1:3]
        all([norm(a) <= norm(b),
        norm(b) <= norm(b+a),
        norm(b) <= norm(b-a),
        norm(b) <= norm(c),
        norm(c) <= norm(c+a),
        norm(c) <= norm(c-a),
        norm(c) <= norm(c+b),
        norm(c) <= norm(c-b),
        norm(c) <= norm(c+b-a),
        norm(c) <= norm(c-b-a),
        norm(c) <= norm(c+b+a),
        norm(c) <= norm(c-b+a)])
    else
       throw(ArgumentError("Can only verify Minkowski reduction in 2D and 3D."))
    end
end


@doc """
    genlat_SQR(a)

Generate a square lattice.

# Arguments
- `a::Real`: the lattice constant

# Returns
- `AbstractArray{<:Real,2}`: the basis of the lattice as
    columns of an array.

# Examples
```jldoctest
using IBZ
a=1
IBZ.Lattices.genlat_SQR(a)
# output
2×2 Array{Int64,2}:
 1  0
 0  1
```
"""
function genlat_SQR(a::Real)::AbstractArray{<:Real,2}
    [a 0; 0 a]
end

@doc """
    genlat_HXG(a)

Generate a 2D hexagonal lattice.

# Arguments
- `a::Real`: the lattice constant

# Returns
- `AbstractArray{<:Real,2}`: the basis of the lattice as
    columns of an array.

# Examples
```jldoctest
using IBZ
a=1
IBZ.Lattices.genlat_HXG(a)
# output
2×2 Array{Float64,2}:
 1.0  -0.5
 0.0   0.866025
```
"""
function genlat_HXG(a::Real)::AbstractArray{<:Real,2}
    Array([a 0; -a/2 a*√3/2]')
end

@doc """
    genlat_REC(a,b)

Generate a rectangular lattice.

# Arguments
- `a::Real`: the lattice constant

# Returns
- `AbstractArray{<:Real,2}`: the basis of the lattice as
    columns of an array.

# Examples
```jldoctest
using IBZ
a=1
b=1.2
IBZ.Lattices.genlat_REC(a,b)
# output
2×2 Array{Float64,2}:
 1.0  0.0
 0.0  1.2
```
"""
function genlat_REC(a::Real,b::Real)::AbstractArray{<:Real,2}
    Array([a 0; 0 b]')
end

@doc """
    genlat_RECI(a,b)

Generate a body-centered rectangular lattice.

# Arguments
- `a::Real`: the lattice constant

# Returns
- `AbstractArray{<:Real,2}`: the basis of the lattice as
    columns of an array.

# Examples
```jldoctest
using IBZ
a=1
b=1.2
IBZ.Lattices.genlat_RECI(a,b)
# output
2×2 Array{Float64,2}:
  0.5  0.5
 -0.6  0.6
```
"""
function genlat_RECI(a::Real,b::Real)::AbstractArray{<:Real,2}
    if a ≈ b
        throw(ArgumentError("The lattice constant must be different sizes for a
            rectangular lattice."))
    end
    Array([a/2 -b/2; a/2 b/2]')
end

@doc """
    genlat_OBL(a,b)

Generate a body-centered rectangular lattice.

# Arguments
- `a::Real`: the lattice constant

# Returns
- `AbstractArray{<:Real,2}`: the basis of the lattice as
    columns of an array.

# Examples
```jldoctest
using IBZ
a=1
b=1.2
θ=π/3
IBZ.Lattices.genlat_OBL(a,b,θ)
# output
2×2 Array{Float64,2}:
 1.0  0.6
 0.0  1.03923
```
"""
function genlat_OBL(a::Real,b::Real,θ::Real)::AbstractArray{<:Real,2}
    if a ≈ b
        throw(ArgumentError("The lattice constant must be different sizes for an
            oblique lattice."))
    end
    if θ ≈ π/2
        throw(ArgumentError("The lattice angle must not equal π/2 for an oblique
            lattice."))
    end
    Array([a 0; b*cos(θ) b*sin(θ)]')
end

@doc """
    genlat_CUB(a)

Generate a simple cubic lattice.

# Arguments
- `a::Real`: the lattice constant

# Returns
- `AbstractArray{<:Real,2}`: the basis of the lattice as
    columns of an array.

# Examples
```jldoctest
using IBZ
a=1
IBZ.Lattices.genlat_CUB(a)
# output
3×3 Array{Int64,2}:
 1  0  0
 0  1  0
 0  0  1
```
"""
function genlat_CUB(a::Real)::AbstractArray{<:Real,2}
    [a 0 0; 0 a 0; 0 0 a]
end

@doc """
    genlat_FCC(a)

Generate a face-centered cubic lattice.

# Arguments
- `a::Real`: the lattice constant

# Returns
- `AbstractArray{<:Real,2}`: the basis of the lattice as
    columns of an array.

# Examples
```jldoctest
using IBZ
a=1
IBZ.Lattices.genlat_FCC(a)
# output
3×3 Array{Float64,2}:
 0.0  0.5  0.5
 0.5  0.0  0.5
 0.5  0.5  0.0
```
"""
function genlat_FCC(a::Real)::AbstractArray{<:Real,2}
    Array([0 a/2 a/2; a/2 0 a/2; a/2 a/2 0]')
end

@doc """
    genlat_BCC(a)

Generate a body-centered cubic lattice.

# Arguments
- `a::Real`: the lattice constant

# Returns
- `AbstractArray{<:Real,2}`: the basis of the lattice as
    columns of an array.

# Examples
```jldoctest
using IBZ
a=1
IBZ.Lattices.genlat_BCC(a)
# output
3×3 Array{Float64,2}:
 -0.5   0.5   0.5
  0.5  -0.5   0.5
  0.5   0.5  -0.5
```
"""
function genlat_BCC(a::Real)::AbstractArray{<:Real,2}
    a/2*Array([-1 1 1; 1 -1 1; 1 1 -1]')
end

@doc """
    genlat_TET(a)

Generate a simple tetragonal lattice.

# Arguments
- `a::Real`: a lattice constant
- `c::Real`: a lattice constant

# Returns
- `AbstractArray{<:Real,2}`: the basis of the lattice as
    columns of an array.

# Examples
```jldoctest
using IBZ
a=1
c=1.2;
IBZ.Lattices.genlat_TET(a,c)
# output
3×3 Array{Float64,2}:
 1.0  0.0  0.0
 0.0  1.0  0.0
 0.0  0.0  1.2
```
"""
function genlat_TET(a::Real,c::Real)::AbstractArray{<:Real,2}
    if a ≈ c
        throw(ArgumentError("The lattices constants must be different for a
            tetragonal lattice."))
    end
    Array([a 0 0; 0 a 0; 0 0 c]')
end

@doc """
    genlat_BCT(a)

Generate a body-centered tetragonal lattice.

# Arguments
- `a::Real`: a lattice constant
- `c::Real`: a lattice constant

# Returns
- `AbstractArray{<:Real,2}`: the basis of the lattice as
    columns of an array.

# Examples
```jldoctest
using IBZ
a=1
c=1.2;
IBZ.Lattices.genlat_BCT(a,c)
# output
3×3 Array{Float64,2}:
 -0.5   0.5   0.5
  0.5  -0.5   0.5
  0.6   0.6  -0.6
```
"""
function genlat_BCT(a::Real,c::Real)::AbstractArray{<:Real,2}
    if a ≈ c
        throw(ArgumentError("The lattices constants must be different for a
            tetragonal lattice."))
    end
    Array([-a/2 a/2 c/2; a/2 -a/2 c/2; a/2 a/2 -c/2]')
end

@doc """
    genlat_ORC(a,b,c)

Generate an orthorhombic lattice.

# Arguments
- `a::Real`: a lattice constant
- `b::Real`: a lattice constant
- `c::Real`: a lattice constant

# Returns
- `AbstractArray{<:Real,2}`: the basis of the lattice as
    columns of an array.

# Examples
```jldoctest
using IBZ
a=1
b=1.4;
c=1.2;
IBZ.Lattices.genlat_ORC(a,b,c)
# output
3×3 Array{Float64,2}:
 1.0  0.0  0.0
 0.0  1.4  0.0
 0.0  0.0  1.2
```
"""
function genlat_ORC(a::Real,b::Real,c::Real)::AbstractArray{<:Real,2}
    if a ≈ b || a ≈ c || b ≈ c
        throw(ArgumentError("The lattices constants must all be different for an
            orthorhombic lattice."))
    end
    Array([a 0 0; 0 b 0; 0 0 c]')
end

@doc """
    genlat_ORCF(a,b,c)

Generate a face-centered orthorhombic lattice.

# Arguments
- `a::Real`: a lattice constant
- `b::Real`: a lattice constant
- `c::Real`: a lattice constant

# Returns
- `AbstractArray{<:Real,2}`: the basis of the lattice as
    columns of an array.

# Examples
```jldoctest
using IBZ
a=1
b=1.4;
c=1.2;
IBZ.Lattices.genlat_ORCF(a,b,c)
# output
3×3 Array{Float64,2}:
 0.0  0.5  0.5
 0.7  0.0  0.7
 0.6  0.6  0.0
```
"""
function genlat_ORCF(a::Real,b::Real,c::Real)::AbstractArray{<:Real,2}
    if a ≈ b || a ≈ c || b ≈ c
        throw(ArgumentError("The lattices constants must all be different for an
            orthorhombic lattice."))
    end
    Array([0 b/2 c/2; a/2 0 c/2; a/2 b/2 0]')
end

@doc """
    genlat_ORCI(a,b,c)

Generate a body-centered orthorhombic lattice.

# Arguments
- `a::Real`: a lattice constant
- `b::Real`: a lattice constant
- `c::Real`: a lattice constant

# Returns
- `AbstractArray{<:Real,2}`: the basis of the lattice as
    columns of an array.

# Examples
```jldoctest
using IBZ
a=1
b=1.4;
c=1.2;
IBZ.Lattices.genlat_ORCI(a,b,c)
# output
3×3 Array{Float64,2}:
 -0.5   0.5   0.5
  0.7  -0.7   0.7
  0.6   0.6  -0.6
```
"""
function genlat_ORCI(a::Real,b::Real,c::Real)::AbstractArray{<:Real,2}
    if a ≈ b || a ≈ c || b ≈ c
        throw(ArgumentError("The lattices constants must all be different for an
            orthorhombic lattice."))
    end
    Array([-a/2 b/2 c/2; a/2 -b/2 c/2; a/2 b/2 -c/2]')
end

@doc """
    genlat_ORCC(a,b,c)

Generate a base-centered orthorhombic lattice.

# Arguments
- `a::Real`: a lattice constant
- `b::Real`: a lattice constant
- `c::Real`: a lattice constant

# Returns
- `AbstractArray{<:Real,2}`: the basis of the lattice as
    columns of an array.

# Examples
```jldoctest
using IBZ
a=1
b=1.2;
c=1.4;
IBZ.Lattices.genlat_ORCC(a,b,c)
# output
3×3 Array{Float64,2}:
  0.5  0.5  0.0
 -0.6  0.6  0.0
  0.0  0.0  1.4
```
"""
function genlat_ORCC(a::Real,b::Real,c::Real)::AbstractArray{<:Real,2}
    if a ≈ b || a ≈ c || b ≈ c
        throw(ArgumentError("The lattices constants must all be different for an
            orthorhombic lattice."))
    end
    Array([a/2 -b/2 0; a/2 b/2 0; 0 0 c]')
end

@doc """
    genlat_HEX(a,c)

Generate a hexagonal lattice.

# Arguments
- `a::Real`: a lattice constant
- `c::Real`: a lattice constant

# Returns
- `AbstractArray{<:Real,2}`: the basis of the lattice as
    columns of an array.

# Examples
```jldoctest
using IBZ
a=1
c=1.2;
IBZ.Lattices.genlat_HEX(a,c)
# output
3×3 Array{Float64,2}:
  0.5       0.5       0.0
 -0.866025  0.866025  0.0
  0.0       0.0       1.2
```
"""
function genlat_HEX(a::Real,c::Real)::AbstractArray{<:Real,2}
    if a ≈ c
        throw(ArgumentError("The lattices constants must be different for a
            hexagonal lattice."))
    end
    Array([a/2 -a*√3/2 0; a/2 a*√3/2 0; 0 0 c]')
end

@doc """
    genlat_RHL(a,α)

Generate a rhombohedral lattice.

# Arguments
- `a::Real`: a lattice constant
- `α::Real`: a lattice angle in radians

# Returns
- `AbstractArray{<:Real,2}`: the basis of the lattice as
    columns of an array.

# Examples
```jldoctest
using IBZ
a=1
α=π/6;
IBZ.Lattices.genlat_RHL(a,α)
# output
3×3 Array{Float64,2}:
  0.965926  0.965926  0.896575
 -0.258819  0.258819  0.0
  0.0       0.0       0.442891
```
"""
function genlat_RHL(a::Real,α::Real)::AbstractArray{<:Real,2}
    if α ≈ π/2
        throw(ArgumentError("The lattice angle cannot be π/2 for a rhombohedral
            lattice."))
    end
    Array([a*cos(α/2) -a*sin(α/2) 0; a*cos(α/2) a*sin(α/2) 0;
        a*cos(α)/cos(α/2) 0 a*√(1-(cos(α)/cos(α/2))^2)]')
end

@doc """
    genlat_MCL(a,b,c,α)

Generate a monoclinic lattice.

# Arguments
- `a::Real`: a lattice constant
- `b::Real`: a lattice constant
- `c::Real`: a lattice constant
- `α::Real`: a lattice angle in radians less than π/2

# Returns
- `AbstractArray{<:Real,2}`: the basis of the lattice as
    columns of an array.

# Examples
```jldoctest
using IBZ
a=1
b=1.2
c=1.4
α=π/6;
IBZ.Lattices.genlat_MCL(a,b,c,α)
# output
3×3 Array{Float64,2}:
 1.0  0.0  0.0
 0.0  1.2  1.21244
 0.0  0.0  0.7
```
"""
function genlat_MCL(a::Real,b::Real,c::Real,α::Real)::AbstractArray{<:Real,2}
    if a ≈ b || a ≈ c || b ≈ c
        throw(ArgumentError("The lattices constants must all be different for a
            monoclinic lattice."))
    end
    if α ≈ π/2 || α > π/2
        throw(ArgumentError("The lattice angle must be less than π/2."))
    end
    Array([a 0 0; 0 b 0; 0 c*cos(α) c*sin(α)]')
end

@doc """
    genlat_MCLC(a,b,c,α)

Generate a base-centered monoclinic lattice.

# Arguments
- `a::Real`: a lattice constant
- `b::Real`: a lattice constant
- `c::Real`: a lattice constant
- `α::Real`: a lattice angle in radians less than π/2

# Returns
- `AbstractArray{<:Real,2}`: the basis of the lattice as
    columns of an array.

# Examples
```jldoctest
using IBZ
a=1
b=1.2
c=1.4
α=π/6;
IBZ.Lattices.genlat_MCLC(a,b,c,α)
# output
3×3 Array{Float64,2}:
 0.5  -0.5  0.0
 0.6   0.6  1.21244
 0.0   0.0  0.7
```
"""
function genlat_MCLC(a::Real,b::Real,c::Real,α::Real)::AbstractArray{<:Real,2}
    if a ≈ b || a ≈ c || b ≈ c
        throw(ArgumentError("The lattices constants must all be different for a
            monoclinic lattice."))
    end
    if α ≈ π/2 || α > π/2
        throw(ArgumentError("The lattice angle must be less than π/2."))
    end
    Array([a/2 b/2 0; -a/2 b/2 0; 0 c*cos(α) c*sin(α)]')
end

@doc """
    genlat_TRI(a,b,c,α,β,γ)

Generate a triclinic lattice.

# Arguments
- `a::Real`: a lattice constant
- `b::Real`: a lattice constant
- `c::Real`: a lattice constant
- `α::Real`: a lattice angle in radians less than π/2

# Returns
- `AbstractArray{<:Real,2}`: the basis of the lattice as
    columns of an array.

# Examples
```jldoctest
using IBZ
a=1
b=1.2
c=1.4
α=π/6;
β=π/3;
γ=π/4;
IBZ.Lattices.genlat_TRI(a,b,c,α,β,γ)
# output
3×3 Array{Float64,2}:
 1.0  0.848528  0.7
 0.0  0.848528  1.01464
 0.0  0.0       0.663702
```
"""
function genlat_TRI(a::Real,b::Real,c::Real,α::Real,β::Real,
    γ::Real)::AbstractArray{<:Real,2}
    if a ≈ b || a ≈ c || b ≈ c
        throw(ArgumentError("The lattices constants must all be different for a
            triclinic lattice."))
    end
    if α ≈ β || α ≈ γ || β ≈ γ
        throw(ArgumentError("No lattice angles can be the same for a triclinic
            lattice."))
    end
    Array([a 0 0; b*cos(γ) b*sin(γ) 0;
            c*cos(β) c/sin(γ)*(cos(α)-cos(β)*cos(γ)) c/sin(γ)*
            √(sin(γ)^2-cos(α)^2-cos(β)^2+2*cos(α)*cos(β)*cos(γ))]')
end

end
