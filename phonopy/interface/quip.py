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
from phonopy.structure.atoms import PhonopyAtoms as Atoms
from phonopy.interface.vasp import check_forces, get_drift_forces


def read_quip(filename):
    pass

def write_quip(filename, atoms):
    pass


class Atoms_with_forces(Atoms):
    # Hack to phonopy.atoms to maintain ASE compatibility for forces
    def get_forces(self):
        return self.forces


def read_quip_output(filename):
    # Read QUIP output and return geometry, energy and forces from prediction

    lines = open(filename, 'r').readlines()

    l = 0
    N = 0
    while l < len(lines):
        line = lines[l]
        if "| Number of atoms" in line:
            N = int(line.split()[5])
        elif "| Unit cell:" in line:
            cell = []
            for i in range(3):
                l += 1
                vec = lmap(float, lines[l].split()[1:4])
                cell.append(vec)
        elif ("Atomic structure:" in line) or ("Updated atomic structure:" in line):
            if "Atomic structure:" in line:
                i_sym = 3
                i_pos_min = 4
                i_pos_max = 7
            elif "Updated atomic structure:" in line:
                i_sym = 4
                i_pos_min = 1
                i_pos_max = 4
            l += 1
            symbols = []
            positions = []
            for n in range(N):
                l += 1
                fields = lines[l].split()
                sym = fields[i_sym]
                pos = map(float, fields[i_pos_min:i_pos_max])
                symbols.append(sym)
                positions.append(pos)
        elif "Total atomic forces" in line:
            forces = []
            for i in range(N):
                l += 1
                force = map(float, lines[l].split()[-3:])
                forces.append(force)
        l += 1

    atoms = Atoms_with_forces(cell=cell, symbols=symbols, positions=positions)
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
