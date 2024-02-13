module SymmetryReduceBZQHullExt

using Base.Iterators: flatten
using LinearAlgebra: dot

using QHull: Chull
using SymmetryReduceBZ.Utilities: get_simplex, sortpts_perm
import SymmetryReduceBZ.Utilities: get_uniquefacetsindices, get_uniquefacets, volume, vertices
import SymmetryReduceBZ.Symmetry: inhull

volume(ch::Chull) = ch.volume
vertices(ch::Chull) = eachrow(ch.points)

function get_uniquefacetsindices(ch::Chull)
    facets = ch.facets
    unique_facets = Vector{eltype(ch.vertices)}[]
    removed=zeros(Bool,size(facets,1))
    for i=1:size(facets,1)
        removed[i] && continue
        removed[i]=true
        face=get_simplex(ch.simplices, i)
        for j=i+1:size(facets,1)
            if isapprox(facets[i,:],facets[j,:],rtol=1e-6)
                removed[j]=true
                append!(face,get_simplex(ch.simplices, j))
            end
        end
        unique!(face)
        # Order the corners of the face either clockwise or counterclockwise.
        permute!(face, sortpts_perm(ch.points[face,:]'))
        push!(unique_facets,face)
    end
    unique_facets
end

function get_uniquefacets(ch::Chull)
    map(get_uniquefacetsindices(ch)) do j
        map(i -> vec(ch.points[i,:]), j)
    end
end

function inhull(point::AbstractVector{<:Real}, chull::Chull;
    rtol::Real=sqrt(eps(float(maximum(flatten(chull.points))))),
    atol::Real=1e-9)

    # hullpts = Array(chull.points')
    distances = chull.facets[:,end]
    norms = chull.facets[:,begin:end-1]'
    inside = true
    for (dist,i)=zip(distances, axes(norms,2))
        s = dot(point + dist*norms[:,i], norms[:,i])

        if !(s <= 0 || isapprox(s,0.0; rtol, atol))
            inside = false
            break
        end
    end
    inside
end


end
