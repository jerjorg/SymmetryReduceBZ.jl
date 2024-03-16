[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](
https://jerjorg.github.io/SymmetryReduceBZ.jl/)
[![Build Status](https://github.com/jerjorg/SymmetryReduceBZ.jl/actions/workflows/CI.yml/badge.svg?branch=master)](
https://github.com/jerjorg/SymmetryReduceBZ.jl/actions/?query=workflow:CI)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)


# SymmetryReduceBZ
The primary purpose of `SymmetryReduceBZ` is to calculate the irreducible
Brillouin zone (IBZ) for crystal structures in 2D or 3D provided the real-space
lattice vectors, atom positions, and atom types. It also contains methods for
making unit cells primitive and lattice reduction. See the User Guide in the
[documentation](https://jerjorg.github.io/SymmetryReduceBZ.jl/) for more details
and usage examples in Python. Details of the algorithm are explained [here](https://arxiv.org/abs/2104.05856).

> **Breaking changes in v0.2**:
See the NEWS.md file for a description of the breaking changes in v0.2.

## Installation

`SymmetryReduceBZ` is a registered Julia package and can be installed using
Julia's package manager `Pkg`.
```
using Pkg
Pkg.add("SymmetryReduceBZ")
```
In order to use the plotting functionality,
you will also need to have the Python Matplotlib library installed, which
[PyCall.jl](https://github.com/JuliaPy/PyCall.jl) can setup automatically via [Conda.jl](https://github.com/JuliaPy/Conda.jl).
If you are installing PyCall, PyPlot, and SymmetryReduceBZ for the first time, just do
`ENV["PYTHON"]=""` before running `Pkg.add(["SymmetryReduceBZ", "PyPlot"])`. Otherwise you
can reconfigure PyCall to use Conda via:
```
ENV["PYTHON"]=""
Pkg.build("PyCall")
```

## Examples

To calculate the irreducible Brillouin zone, provide the lattice and atomic
basis to `calc_ibz`. The IBZ will be returned as a polyhedron from Polyhedra.jl,
which can be viewed either as a convex hull or an intersection of half spaces.
```@example
import SymmetryReduceBZ.Lattices: genlat_CUB
import SymmetryReduceBZ.Symmetry: calc_ibz
a = 2.0
real_latvecs = genlat_CUB(a)
atom_types = [0,0]
atom_pos = Array([0 0 0; 0.5 0.5 0.5]')
coordinates = "Cartesian"
makeprim = false
convention = "ordinary"
ibz = calc_ibz(real_latvecs,atom_types,atom_pos,coordinates,
  makeprim,convention)
```
The arguments for `calc_ibz` are as follows:
- `real_latvecs`: the real-space lattice vectors as columns of a matrix.
- `atom_types`: a vector of atom types as integers.
- `atom_pos`: the positions of atoms in the crystal structure as columns of a matrix.
- `coordinates`: indicates the atoms are in `"lattice"` or `"Cartesian"`
	coordinates.
- `makeprim`: make the unit cell primitive before calculating the IBZ if
	`true`.
- `convention`: the convention used to go between real and reciprocal space. The
	two conventions are `"ordinary"` (temporal) frequency and `"angular"`
	frequency.
- `library::Polyhedra.Library=CDDLib.Library()`: a polyhedron manipulation library
- `rtol=sqrt(eps(float(maximum(real_latvecs))))`: (optional) a relative tolerance for floating-point comparisons.
- `atol=1e-9`: (optional) an absolute tolerance for floating-point comparisons.	
	
The vertices of the ibz are accessed with
`SymmetryReduceBZ.Utilities.vertices(ibz)` as an iterator of vectors. Other
attributes are accessible such as `SymmetryReduceBZ.Utilities.volume(ibz)`. The
faces of the IBZ are calculated with
```
import SymmetryReduceBZ.Utilities: get_uniquefacets
facets = get_uniquefacets(ibz)
```
`facets` is a list of points at the corners of each facet. The function
`get_uniquefacets` returns a list of the points that lie on each facet.
See the [documentation](https://juliapolyhedra.github.io/Polyhedra.jl/stable/)
for more details about the polyhedral interface.

The function `plot_convexhulls` is useful for visualizing the Brillouin zone
and irreducible Brillouin zone. The arguments are the same as those from
`calc_ibz`.
```@example
ENV["MPLBACKEND"]="qt5agg"
using PyPlot
import SymmetryReduceBZ.Plotting: plot_convexhulls
import SymmetryReduceBZ.Lattices: genlat_CUB
a = 2.0
real_latvecs = genlat_CUB(a)
atom_types = [0,0]
atom_pos = Array([0 0 0; 0.5 0.5 0.5]')
coordinates = "Cartesian"
makeprim = false
convention = "ordinary"
ax=plot_convexhulls(real_latvecs,atom_types,atom_pos,coordinates,
  makeprim,convention)
```
![IBZ](https://github.com/jerjorg/SymmetryReduceBZ.jl/blob/master/plots/ibz.png)

The functions `plot_2Dconvexhull` and `plot_3Dconvexhull` allow greater customization of 
the appearance of the convex hull.

```@example
ENV["MPLBACKEND"]="qt5agg"
using PyPlot
import SymmetryReduceBZ.Symmetry: calc_bz, calc_ibz
import SymmetryReduceBZ.Plotting: plot_2Dconvexhull
real_latvecs = [1 0; 0 1]
convention="ordinary"
atom_types=[0]
atom_pos = Array([0 0]')
coords = "Cartesian"
makeprim=false
bz = calc_bz(real_latvecs,atom_types,atom_pos,coords,makeprim,convention)
ibz = calc_ibz(real_latvecs,atom_types,atom_pos,coords,makeprim,convention)
ax = plot_2Dconvexhull(bz,facecolor="deepskyblue",linewidth=3,edgecolor="cyan",alpha=0.2)
ax = plot_2Dconvexhull(ibz,ax;facecolor="coral",linewidth=3,edgecolor="magenta",alpha=0.4)
axis("off")
```
![IBZ](https://github.com/jerjorg/SymmetryReduceBZ.jl/blob/master/plots/ibz-3.png)

```@example
ENV["MPLBACKEND"]="qt5agg"
using PyPlot
import SymmetryReduceBZ.Symmetry: calc_bz, calc_ibz
import SymmetryReduceBZ.Plotting: plot_3Dconvexhull
real_latvecs = [1 0 0; 0 1 0; 0 0 1]
convention="ordinary"
atom_types=[0]
atom_pos = Array([0 0 0]')
coords = "Cartesian"
makeprim=false
bz = calc_bz(real_latvecs,atom_types,atom_pos,coords,makeprim,convention)
ibz = calc_ibz(real_latvecs,atom_types,atom_pos,coords,makeprim,convention)
fig = figure()
ax = fig.add_subplot(111, projection="3d")
ax = plot_3Dconvexhull(ibz,ax,facecolors="pink",alpha=1,edgecolors="black",linewidths = 1)
ax = plot_3Dconvexhull(bz,ax,facecolors="deepskyblue",edgecolors="white",linewidths=1,alpha=0.2)
axis("off")
```
![IBZ](https://github.com/jerjorg/SymmetryReduceBZ.jl/blob/master/plots/ibz-2.png)
