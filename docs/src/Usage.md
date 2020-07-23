# User guide

`IBZ.jl` calculates the irreducible Brillouin zone (IBZ) of a crystal structure
in 2D and 3D. It also contains functions related to the symmetry of lattices and
lattice reduction. Later on, it will be able to calculate high symmetry points
and paths within the IBZ.

To calculate the IBZ, simply provide the lattice and atomic basis to `calc_ibz`.
The IBZ will be returned as either a convex hull or intersection of half spaces.
```@example
using IBZ # hide
a = 1.0;
real_latvecs = genlat_CUB(a)
atom_types = [0,0]
atom_pos = Array([0 0 0; 0.5 0.5 0.5]')
ibzformat = "convex hull"
coordinates = "Cartesian"
convention = "ordinary"
ibz = calc_ibz(real_latvecs,atom_types,atom_pos,coordinates,ibzformat,
  convention)
```
The columns of `real_latvecs` are the lattice generating vectors, the columns
of `atom_pos` are the positions of the atoms in Cartesian coordinates (in this
case), `coordinates` are the coordinates of the atom positions, and `convention`
gives the convention for going from real to reciprocal space (whether or not to
multiply by 2π). There is a simple function for visualizing the IBZ along with
the Brillouin zone (BZ).
```@example
using IBZ # hide
a = 1.0 # hide
real_latvecs = genlat_CUB(a) # hide
atom_types = [0,0] # hide
atom_pos = Array([0 0 0; 0.5 0.5 0.5]') # hide
coordinates = "Cartesian" # hide
convention = "ordinary" # hide
plot_cells(real_latvecs,atom_types,atom_pos,coordinates,convention)
```
