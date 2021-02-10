[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](
https://jerjorg.github.io/SymmetryReduceBZ.jl/)
[![Build Status](
https://travis-ci.com/jerjorg/SymmetryReduceBZ.jl.svg?branch=master)](
https://travis-ci.com/jerjorg/SymmetryReduceBZ.jl)
[![Coverage Status](
https://coveralls.io/repos/github/jerjorg/SymmetryReduceBZ.jl/badge.svg?branch=master)](
https://coveralls.io/github/jerjorg/SymmetryReduceBZ.jl?branch=master)


# SymmetryReduceBZ
The primary purpose of `SymmetryReduceBZ` is to calculate the irreducible
Brillouin zone (IBZ) for crystal structures in 2D or 3D provided the real-space
lattice vectors, atomic positions, and atom types. It also contains methods for
making unit cells primitive and lattice reduction. See the User Guide in the
[documentation](https://jerjorg.github.io/SymmetryReduceBZ.jl/) for more details
and usage examples.

## Installation

`SymmetryReduceBZ` is a registered Julia package and can be installed using
Julia's package manager `Pkg`.
```
using Pkg
Pkg.add("SymmetryReduceBZ")
```

## Examples

To calculate the irreducible Brillouin zone, provide the lattice and atomic
basis to `calc_ibz`. The IBZ will be returned as either a convex hull or
intersection of half spaces.
```@example
import SymmetryReduceBZ.Lattices: genlat_CUB
import SymmetryReduceBZ.Symmetry: calc_ibz
a = 2.0
real_latvecs = genlat_CUB(a)
atom_types = [0,0]
atom_pos = Array([0 0 0; 0.5 0.5 0.5]')
ibzformat = "convex hull"
coordinates = "Cartesian"
makeprim = false
convention = "ordinary"
ibz = calc_ibz(real_latvecs,atom_types,atom_pos,coordinates,ibzformat,
  makeprim,convention)
```
The arguments for `calc_ibz` are as follows:
- `real_latvecs`: the real-space lattice vectors as columns of an array.
- `atom_types`: a list of atom types as integers.
- `atom_pos`: the positions of atoms in the crystal structure as columns of an
	array.
- `coords`: the positions of the atoms in \"lattice\" or \"Cartesian\"
	coordinates.
- `ibzformat`: the format of the irreducible Brillouin zone. Options include
	\"convex hull\" and \"half-space\".
- `convention`: the convention used to go between real and reciprocal space. The
	two conventions are \"ordinary\" (temporal) frequency and \"angular\"
	frequency.
- `primitive`: make the unit cell primitive before calculating the the IBZ if
	true.
	
The vertices of the ibz are accessed with `ibz.points[ibz.vertices,:]`. The
vertices of the IBZ and `ibz.points` should be the same. The rows of the array
are the vertices of the IBZ in Cartesian coordinates. Other attributes of the
IBZ are accessible, such as the volume `ibz.volume`. The faces of the IBZ are
calculated with
```
import SymmetryReduceBZ.Utilities: get_uniquefacets
indices = get_uniquefacets(ibz)
facets = [ibz.points[ind] for ind=indices]
```
`facets` is a list of points at the corners of each facet.The function
`get_uniquefacets` returns the indices of points that lie on the same facet. The
facets are available to through `ibz` as simplices, but often multiple simplices
lie on the same facet.

The functions `plot_convexhulls` is useful for visualizing the Brillouin zone
and irreducible Brillouin zone. The arguments are the same as those from
`calc_ibz`.
```@example
import SymmetryReduceBZ.Plotting: plot_convexhulls
import SymmetryReduceBZ.Lattices: genlat_CUB
a = 2.0
real_latvecs = genlat_CUB(a)
atom_types = [0,0]
atom_pos = Array([0 0 0; 0.5 0.5 0.5]')
coordinates = "Cartesian"
makeprim = false
convention = "ordinary"
(fig,ax)=plot_convexhulls(real_latvecs,atom_types,atom_pos,coordinates,
  makeprim,convention)
```
