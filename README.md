# quotate
Rotate Lattices so that Quantum Espresso recognises their symmetries.

The spglib library (when it symmetrises a structure) tends to write the primitive {B,F}CC lattices with the a vector pointing along the Cartesian x-axis:
e.g.
```
%BLOCK LATTICE_CART
    2.74336000000000    0.00000000000000    0.00000000000000
   -0.91445330469095    2.58646462244118    0.00000000000000
   -0.91445330469095   -1.29323225046092    2.23994410410340
%ENDBLOCK LATTICE_CART
```
Quantum Espresso expects these lattices to look like, otherwise it won't recognise the symmetry:
```
  2          cubic F (fcc)
      v1 = (a/2)(-1,0,1),  v2 = (a/2)(0,1,1), v3 = (a/2)(-1,1,0)

  3          cubic I (bcc)
      v1 = (a/2)(1,1,1),  v2 = (a/2)(-1,1,1),  v3 = (a/2)(-1,-1,1)
```

Quotate rotates any given lattice so it will be have the appropriate symmetry in quantum espresso.

Since spglib (and codes which leverage it) can deal with symmetry for arbitrary rotations, we're better off using the quantum-espresso convention everywhere...

## Usage
```
usage: quotate.py [-h] [--ibrav IBRAV] [-F OFORMAT] input output

positional arguments:
  input                 Input file name (xyz format)
  output                Output file name (xyz format)

optional arguments:
  -h, --help            show this help message and exit
  --ibrav IBRAV         pw.x ibrav value
  -F OFORMAT, --oformat OFORMAT
                        The output file format
```

## Examples

```
./quotate.py tests/H3S-in.cell tests/H3S-out.cell --ibrav 3
```

```
./quotate.py tests/SiO2-in.cell tests/SiO2-out.cell --ibrav 2
```

```
./quotate.py tests/Cd6Yb-in.poscar tests/Cd6Yb-out.cell --ibrav 3
```

```
./quotate.py tests/SiO2-in.cell tests/SiO2.scf.in --ibrav 2 -F espresso-in
```

## Test
```
./run_test.sh
```
and inspect outputs in /tests/
