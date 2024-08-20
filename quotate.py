#! /usr/bin/env python
from ase.io import read, write
from ase.build import make_supercell
from ase.geometry import get_duplicate_atoms
from ase.spacegroup import get_spacegroup
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

def spg2ibrav(spg):

    family = no2family(spg.no)

    lattice = spg.lattice

    ibrav_dict = {
            # SC
            1: ['cubic', 'P'],
            # FCC
            2: ['cubic', 'F'],
            # BCC
            3: ['cubic', 'I'],
            # Hexagonal or Trigonal P
            4: ['hexagonal', 'P'],
            # Trigonal R
            5: ['trigonal', 'P'],
            # ST
            6: ['tetragonal', 'P'],
            # BCT
            7: ['tetragonal', 'I'],
            # Orthorhombic
            8: ['orthorhombic', 'P'],
            # Base Centred Orthorhombic
            9: ['orthorhombic', 'A'],
            # Face Centred Orthorhombic
            10: ['orthorhombic', 'F'],
            # Body Centred Orthorhombic
            11: ['orthorhombic', 'I'],
            # Monoclinic
            12: ['monoclinic', 'P'],
            # Monoclinic Base Centred
            13: ['monoclinic', 'A'],
            # Triclinic
            14: ['triclinic', 'P']}

    print(family, lattice)
    for ibrav in ibrav_dict.keys():
        pair = ibrav_dict[ibrav]
        if family == pair[0] and lattice == pair[1]:
            return ibrav
    return -1



def no2family(no):
    space_group_to_lattice = {
            'cubic':        range(195, 231),
            'hexagonal':    range(168, 195),
            'trigonal':     range(143, 168),
            'tetragonal':   range(75, 143),
            'orthorhombic': range(16, 75),
            'monoclinic':   range(3, 16),
            'triclinic':    range(1, 3),
    }
    for i in space_group_to_lattice.keys():
        if no in space_group_to_lattice[i]:
            return i


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="Input file name (any ASE supported format)")
    parser.add_argument("output", help="Output file name (any ASE supported format)")
#    parser.add_argument("--setting", default=1, type=int,
#                        help="""
#                        Spacegroup setting
#                        (ie. 1 = hexagonal basis, 2 = rhombohedral)
#                        """)
    parser.add_argument("--ibrav", default=None, type=int,
                        help="pw.x ibrav value")
    parser.add_argument('-F', '--oformat', default='',
                        help='''
                        The output file format (If ASE can not deduce from <input>)
                        '''
                        )
    return parser.parse_args()


# Parse the command line arguments
args = parse_arguments()

# Read in the structure
structure = make_supercell(read(args.input), np.eye(3)*3, wrap=False)

structure_spglib = ase_atoms_to_spglib_cell(structure)

structure_spglib = spglib.standardize_cell(structure_spglib, symprec=0.1)

structure = spglib_cell_to_ase_atoms(structure_spglib)

spg = get_spacegroup(structure, symprec=0.1)

if args.ibrav is None:
    lattice_dict = {'F': 2, 'I': 3, 'R': 5}
    try:
        args.ibrav = lattice_dict[spg.lattice]
        args.ibrav = spg2ibrav(spg)
        print(spg2ibrav(spg), spg.no)
    except KeyError:
        args.ibrav = -1

if args.ibrav == 2:
    a = structure.cell.cellpar()[0]
    M = a/2 * np.array([[-1, 0, 1], [0, 1, 1], [-1, 1, 0]])
    M2 = structure.get_cell()
    pos_tmp = structure.get_positions()@np.linalg.inv(M)

elif args.ibrav == 3:
    a = structure.cell.cellpar()[0]
    structure = make_supercell(structure, np.eye(3)*3, wrap=False)
    M = a/2 * np.array([[1, 1, 1], [-1, 1, 1], [-1, -1, 1]])
    M2 = structure.get_cell()
    pos_tmp = structure.get_positions()@np.linalg.inv(M)
elif args.ibrav == 5:
    '''
      5          Trigonal R, 3fold axis c        celldm(4)=cos(gamma)
      The crystallographic vectors form a three-fold star around
      the z-axis, the primitive cell is a simple rhombohedron:
      v1 = a(tx,-ty,tz),   v2 = a(0,2ty,tz),   v3 = a(-tx,-ty,tz)
      where c=cos(gamma) is the cosine of the angle gamma between
      any pair of crystallographic vectors, tx, ty, tz are:
        tx=sqrt((1-c)/2), ty=sqrt((1-c)/6), tz=sqrt((1+2c)/3)
    '''
    structure_spglib = ase_atoms_to_spglib_cell(structure)

    structure_spglib = spglib.standardize_cell(
        structure_spglib, symprec=0.1, to_primitive=True)

    structure_prim = spglib_cell_to_ase_atoms(structure_spglib)

    a = structure_prim.cell.cellpar()[0]
    gamma = structure_prim.cell.cellpar()[5]
    c = np.cos(gamma*np.pi/180)

    tx = np.sqrt((1-c)/2)
    ty = np.sqrt((1-c)/6)
    tz = np.sqrt((1+2*c)/3)

    structure = make_supercell(structure_prim, np.eye(3), wrap=False)
    structure.wrap()
    M = a * np.array([[tx, -ty, tz], [0, 2*ty, tz], [-tx, -ty, tz]])
    M2 = structure.get_cell()

    # I'm not sure why this is M2!!
    pos_tmp = structure.get_positions()@np.linalg.inv(M2)
elif args.ibrav == 7:
    '''
      7          Tetragonal I (bct)              celldm(3)=c/a
      v1=(a/2)(1,-1,c/a),  v2=(a/2)(1,1,c/a),  v3=(a/2)(-1,-1,c/a)
    '''
    structure_spglib = ase_atoms_to_spglib_cell(structure)

    structure_spglib = spglib.standardize_cell(
        structure_spglib, symprec=0.1, to_primitive=False)

    structure_conv = spglib_cell_to_ase_atoms(structure_spglib)

    a = structure_conv.cell.cellpar()[0]
    b = structure_conv.cell.cellpar()[1]
    c = structure_conv.cell.cellpar()[2]
    structure = make_supercell(structure, np.eye(3)*3, wrap=False)
    M = np.array([[1, -1, c/a], [1, 1, c/a], [-1, -1, c/a]]) * a/2
    M2 = structure.get_cell()
    pos_tmp = structure.get_positions()@np.linalg.inv(M)
elif args.ibrav == 11:
    '''
       11          Orthorhombic body-centered      celldm(2)=b/a
                                             celldm(3)=c/a
      v1=(a/2,b/2,c/2),  v2=(-a/2,b/2,c/2),  v3=(-a/2,-b/2,c/2)

    '''
    structure_spglib = ase_atoms_to_spglib_cell(structure)

    structure_spglib = spglib.standardize_cell(
        structure_spglib, symprec=0.1, to_primitive=False)

    structure_conv = spglib_cell_to_ase_atoms(structure_spglib)

    a = structure_conv.cell.cellpar()[0]
    b = structure_conv.cell.cellpar()[1]
    c = structure_conv.cell.cellpar()[2]
    structure = make_supercell(structure, np.eye(3)*3, wrap=False)
    M = np.array([[a, b, c], [-a, b, c], [-a, -b, c]]) * 0.5
    M2 = structure.get_cell()
    pos_tmp = structure.get_positions()@np.linalg.inv(M)
else:
    print("ibrav not supported")
    exit(1)

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
