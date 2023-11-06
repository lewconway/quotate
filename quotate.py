#! /usr/bin/env python3.9
from ase.io import read, write
from ase.build import make_supercell
from ase.geometry import get_duplicate_atoms
from ase import Atoms
import argparse
import numpy as np
import spglib


def ase_atoms_to_spglib_cell(atoms):
    """Convert ase atoms object to spglib cell tuple
    cell = (lattice, positions, numbers, magmoms)

    :atoms: ase.atoms
    :returns: cell

    """
    return (atoms.get_cell(), atoms.get_scaled_positions(), atoms.get_atomic_numbers())


def spglib_cell_to_ase_atoms(cell):
    """Convert spglib cell tuple to ase atoms object

    :cell: cell
    :returns: ase.Atoms

    """
    return Atoms(symbols=cell[2], scaled_positions=cell[1], cell=cell[0], pbc=True)


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="Input file name (xyz format)")
    parser.add_argument("output", help="Output file name (xyz format)")
    parser.add_argument("--setting", default=1, type=int,
                        help="""
                        Spacegroup setting
                        (ie. 1 = hexagonal basis, 2 = rhombohedral)
                        """)
    parser.add_argument("--ibrav", default=-1, type=int,
                        help="pw.x ibrav value")
    parser.add_argument('-F', '--oformat', default='', help='The output file format')
    return parser.parse_args()


# Parse the command line arguments
args = parse_arguments()

# Read in the structure
structure = make_supercell(read(args.input), np.eye(3)*3, wrap=False)

structure_spglib = ase_atoms_to_spglib_cell(structure)

structure_spglib = spglib.standardize_cell(structure_spglib, symprec=0.1)

structure = spglib_cell_to_ase_atoms(structure_spglib)

if args.ibrav == 2:
    a = structure.cell.cellpar()[0]
    M = a/2 * np.array([[-1, 0, 1], [0, 1, 1], [-1, 1, 0]])
    M2 = structure.get_cell()

elif args.ibrav == 3:
    a = structure.cell.cellpar()[0]
    structure = make_supercell(structure, np.eye(3)*3, wrap=False)
    M = a/2 * np.array([[1, 1, 1], [-1, 1, 1], [-1, -1, 1]])
    M2 = structure.get_cell()

else:
    print("ibrav not supported")
    exit()

pos_tmp = structure.get_positions()@np.linalg.inv(M)
structure.set_cell(M)
structure.set_scaled_positions(pos_tmp)
structure.wrap()
get_duplicate_atoms(structure, cutoff=0.01, delete=True)

if args.oformat == '':
    write(args.output, structure, parallel=False)
elif args.oformat == 'res' and 'res' in args.input:
    write(args.output, structure, write_info=False)
else:
    write(args.output, structure, parallel=False, format=args.oformat)
