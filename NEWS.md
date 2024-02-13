# SymmetryReduceBZ.jl v0.2

## `bzformat` is removed
Instead, a Polyhedra.jl polyhedron is returned. See that package for further
details on the polyhedral interface, which is more feature complete. To recover
the previous formats from the new one, use the following recipe:
- `bzformat="half-space"`: `poly = calc_bz(...); using Polyhedra; hrep(poly)`
- `bzformat="convex hull"`: `poly = calc_bz(...); using Polyhedra, QHull; chull(permutedims(reduce(hcat, points(poly))))`
Any `SymmetryReduceBZ.Utilities` or `SymmetryReduceBZ.Plotting` functions that
used to work with `bzformat="convex hull"` should still work with the old format.

## Python dependencies moved to extensions
In Julia versions supporting package extensions (v1.9 and above)
SymmetryReduceBZ.jl no longer has Python dependencies by default, simplifying
the package installation. This means QHull.jl and PyPlot.jl are package
extensions, and that the plotting functions of the library require `using
PyPlot` before calling the functions.

## `get_uniquefacets` now returns facets
Previously, this function would return indices of facets, which is now handled
by a `get_uniquefacetsindices` routine.