# User guide

## Julia

To calculate the IBZ, simply provide the lattice and atomic basis to `calc_ibz`.
The IBZ will be returned as either a convex hull or intersection of half spaces.
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
The columns of `real_latvecs` are the lattice generating vectors, the columns
of `atom_pos` are the positions of the atoms in Cartesian coordinates (in this
case), `coordinates` are the coordinates of the atom positions, and `convention`
gives the convention for going from real to reciprocal space (whether or not to
multiply by 2π). There is a simple function for visualizing the IBZ along with
the Brillouin zone (BZ).
```@example
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

## Python
It is possible to call `SymmetryReduceBZ` from Python using 
[PyJulia](https://pyjulia.readthedocs.io/en/latest/). In a
Jupyter notebook, IPython magic can be used to directly evaluate Julia commands
in notebook cells.
```
from julia import Julia
jl = Julia(runtime="/usr/local/bin/julia")
```
```
%%julia
using Pkg
Pkg.add("SymmetryReduceBZ")
Pkg.add("PyPlot")
using PyPlot
import SymmetryReduceBZ.Plotting: plot_convexhulls
import SymmetryReduceBZ.Lattices: genlat_CUB

a = 1.0
real_latvecs = genlat_CUB(a)
atom_types = [0,1]
atom_pos = Array([0 0 0; 0.5 0.5 0.5]')
coordinates = "Cartesian"
primitive = true
makeprim = false
convention = "ordinary"
ax=plot_convexhulls(real_latvecs,atom_types,atom_pos,coordinates,makeprim,
   convention)
```

Without IPython magic, Julia commands can be called with the `eval` method of
`PyJulia`.
```
from julia import Julia
jl = Julia(runtime="/usr/local/bin/julia")

jl.eval("""
using Pkg
Pkg.add(\"SymmetryReduceBZ\")
Pkg.add(\"PyPlot\")
using PyPlot
import SymmetryReduceBZ.Plotting: plot_convexhulls
import SymmetryReduceBZ.Lattices: genlat_CUB
a = 1.0
real_latvecs = genlat_CUB(a)
atom_types = [0,1]
atom_pos = Array([0 0 0; 0.5 0.5 0.5]')
coordinates = "Cartesian"
makeprim = true
convention = "ordinary"
ax=plot_convexhulls(real_latvecs,atom_types,atom_pos,coordinates,makeprim,
   convention)
""")
```
