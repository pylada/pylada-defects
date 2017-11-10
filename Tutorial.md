## Defect automation code *pylada_defects.py*

---
## Functions specific to Pylada

#### Structure manipulation functions

> More details about the functions below can be found on the [documentation page](http://pylada.github.io/pylada/userguide/crystal.html)

```
pylada.crystal.Structure
```
> - Class *Structure* in Module *crystal*
> - Defines structure

```
pylada.crystal.primitive
```
> - Returns the primitive unit-cell lattice of a input lattice
> - :param lattice: `Structure`, :param tolerance: `tolerance` when comparing position. Default to 1e-8

```
pylada.crystal.read
```
> - Methods to read structures from file
> - structure = read.poscar(file=’POSCAR’)

```
pylada.crystal.write
```
> - Methods to write structures from file
> - write.poscar(structure, file, vasp5=True)

```
pylada.crystal.neighbors
```
> - Returns list of neighbors to input position
> - :param lattice: `Structure`, :param nmax: `Integer` number of first neighbors to search for, 
          :param center: `Position for which to determine first neighbors`, :param tolerance: `tolerance`
> - :returns: `A list of 3-tuples`. The first item is a refence to the
          neighboring atom, the second is the position of its
          relevant periodic image *relative* to the center, the
          third is its distance from the center
```
pylada.crystal.supercell
```
> - Creates a supercell of an input lattice
> - :param lattice:`Structure`, :param cell: `Cell in cartesian coordinates of the supercell`
> - :returns: `Structure` representing the supercell

```
pylada.crystal.defects.reindex_sites
```
> - Reindexes atoms of structure according to lattice sites
> - reindex_sites(structure, lattice, tolerance=0.5)

---
#### Functions handling VASP output

> More details about the functions below can be found on the [documentation page](http://pylada.github.io/pylada/pyapi/vasp/extract.html)

```
pylada.vasp.Extract
```
> - Class *Extract* in Module *vasp*
> - Initializes extraction object to extract DFT data from OUTCAR
> - :param directory: `Directory` where the OUTCAR resides
> - extract=Extract(’path/to/the/folder/’)
> - Other Properties of the Class *Extract*

```
Extract.structure
Extract.fermi_energy
Extract.electropot
```
> - Greps average atomic electrostatic potentials from OUTCAR
> - reutrns a numpy array

```
Extract.eigenvalues
```
> - The eigenvalues are return as the k×n numpy array with k being the number of k-points and n the number of bands


---
#### Functions handling units
```
pylada.physics.a0
pylada.physics.Ry
pylada.ewald.ewald
```
---
#### Specific functions
```
pylada.crystal.third_order_cc
```
