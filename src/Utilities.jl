module Utilities

export unique

@doc """
    unique(points,rtol,atol)

Remove duplicate points from an array.

# Arguments
-`points::AbstractArray{<:Real,2}`: the points are columns of a 2D array.
-`rtol::Real=sqrt(eps(float(maximum(points))))`: a relative tolerance for floating
    point comparisons.
-`atol::Real=0.0`: an absolume tolerance for floating point comparisons.

# Returns
-`uniquepts::AbstractArray{<:Real,2}`: a 2D array of unique points as columns.

# Examples
```jldoctest
using IBZ
points=Array([1 2; 2 3; 3 4; 1 2]')
IBZ.Utilities.unique(points)
# output
2×3 Array{Int64,2}:
 1  2  3
 2  3  4
```
"""
function unique(points::AbstractArray{<:Real,2},
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
