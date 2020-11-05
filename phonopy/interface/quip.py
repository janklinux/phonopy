# quip.py - IO routines for phonopy using QUIP
#
# Copyright (C) 2009-2011 Jan Kloppenburg
# All rights reserved.
#
# This file is part of phonopy.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#
# * Redistributions of source code must retain the above copyright
#   notice, this list of conditions and the following disclaimer.
#
# * Redistributions in binary form must reproduce the above copyright
#   notice, this list of conditions and the following disclaimer in
#   the documentation and/or other materials provided with the
#   distribution.
#
# * Neither the name of the phonopy project nor the names of its
#   contributors may be used to endorse or promote products derived
#   from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
# FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
# COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
# ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#

import sys
import numpy as np
from phonopy.structure.atoms import PhonopyAtoms
from phonopy.interface.vasp import check_forces, get_drift_forces


def read_quip(filename):
    with open(filename, 'r') as f:
        lnr = 0
        pos = []
        name = []
        for line in f:
            if lnr == 1:
                bits = line.split()
                for ib, bit in enumerate(bits):
                    if 'Lattice' in bit:
                        cell = [[float(bit.split('=')[1][1:]), float(bits[ib+1]), float(bits[ib+2])],
                                [float(bits[ib+3]), float(bits[ib+4]), float(bits[ib+5])],
                                [float(bits[ib+6]), float(bits[ib+7]), float(bits[ib+8][:-1])]]
            if lnr >= 2:
                name.append(line.split()[0])
                pos.append([float(x) for x in line.split()[1:4]])
            lnr += 1
    return PhonopyAtoms(symbols=name, cell=cell, positions=pos)


def write_quip(filename, atoms):
    with open(filename, 'w') as f:
        f.write('{}\n'.format(len(atoms.symbols)))
        f.write('Lattice="{} {} {} {} {} {} {} {} {}"\n'.format(atoms.cell[0][0], atoms.cell[0][1], atoms.cell[0][2],
                                                                atoms.cell[1][0], atoms.cell[1][1], atoms.cell[1][2],
                                                                atoms.cell[2][0], atoms.cell[2][1], atoms.cell[2][2]))
        for n, a in zip(atoms.symbols, atoms.positions):
            f.write('{} {} {} {}\n'.format(n, a[0], a[1], a[2]))

class Atoms_with_forces(PhonopyAtoms):
    # Hack to phonopy.atoms to maintain ASE compatibility for forces
    def get_forces(self):
        return self.forces


def read_quip_output(filename):
    # Read QUIP output and return geometry, energy and forces from prediction
    with open(filename, 'r') as f:
        lnr = 0
        name = []
        pos = []
        forces = []
        for line in f:
            if line.startswith('AT'):
                if lnr == 0:
                    n_atoms = int(line.split()[1])
                if lnr == 1:
                    bits = line.split()
                    for ib, bit in enumerate(bits):
                        if 'Lattice' in bit:
                            cell = [[float(bit.split('=')[1][1:]), float(bits[ib + 1]), float(bits[ib + 2])],
                                    [float(bits[ib + 3]), float(bits[ib + 4]), float(bits[ib + 5])],
                                    [float(bits[ib + 6]), float(bits[ib + 7]), float(bits[ib + 8][:-1])]]
                if lnr >= 2:
                    name.append(line.split()[1])
                    pos.append([float(x) for x in line.split()[2:5]])
                    forces.append([float(x) for x in line.split()[10:13]])
                lnr += 1

    atoms = Atoms_with_forces(cell=cell, symbols=name, positions=pos)
    atoms.forces = forces

    return atoms


def write_supercells_with_displacements(supercell,
                                        cells_with_disps,
                                        ids,
                                        pre_filename="geometry.xyz",
                                        width=3):
    """Writes perfect supercell and supercells with displacements

    Args:
        supercell: perfect supercell
        cells_with_disps: supercells with displaced atoms
        filename: root-filename
    """

    # original cell
    write_quip(pre_filename + ".supercell", supercell)

    # displaced cells
    for i, cell in zip(ids, cells_with_disps):
        filename = "{pre_filename}-{0:0{width}}".format(
            i, pre_filename=pre_filename, width=width)
        write_quip(filename, cell)


def parse_set_of_forces(num_atoms, forces_filenames, verbose=True):
    """parse the forces from output files in `forces_filenames`"""
    is_parsed = True
    force_sets = []
    for i, filename in enumerate(forces_filenames):
        if verbose:
            sys.stdout.write("%d. " % (i + 1))

        atoms = read_quip_output(filename)
        forces = atoms.forces
        if check_forces(forces, num_atoms, filename, verbose=verbose):
            drift_force = get_drift_forces(forces, filename=filename, verbose=verbose)
            force_sets.append(np.array(forces) - drift_force)
        else:
            is_parsed = False

    if is_parsed:
        return force_sets
    else:
        return []
