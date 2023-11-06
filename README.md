# quotate
Rotate Lattices so that Quantum Espresso recognises their symmetries

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
