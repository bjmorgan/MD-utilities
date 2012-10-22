# `polyhedra`

This program reads in a set of coordinates for an optimised close-packed lattice, and outputs the list of ion numbers that define the vertices of all octahedra and tetrahedra, and the coordinates of the centres of these polyhedra.

## Input files:

* `polyhedra.inpt`

	Example:
	
	    Cu_lattice.pos
	    864                  natoms
	    152.477600000000   46.6869000000000    40.4324250000000 
	    10.0                 rcut
	    1.0 0.0 0.0          
	    0.0 1.0 0.0
	    0.0 0.0 1.0
	    1.0 0.0 0.0          close-packed vector`
* Coordinates file
	
	In this case `Cu_lattice.pos` with the format
	
	    2.95694091770000        15.5615063227000        8.98311445680000
	    2.95694091770000        15.5615063227000        22.4605894568000
	    2.95694091770000        15.5615063227000        35.9380644568000
	    …

## Output files:

* `tet1.list`
* `tet2.list`
* `oct.list`
* `tet1_c.out`
* `tet2_c.out`
* `oct_c.out`

