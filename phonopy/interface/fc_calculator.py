# Copyright (C) 2019 Atsushi Togo
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

from phonopy.interface.alm import get_fc2 as get_fc2_alm


def get_fc2(supercell,
            primitive,
            displacements,
            forces,
            fc_calculator=None,
            atom_list=None,
            fc_options=None,
            log_level=0):
    """Supercell 2nd order force constants (fc2) are calculated.

    The expected shape of supercell fc2 to be returned is
        (len(atom_list), num_atoms, 3, 3).

    Parameters
    ----------
    supercell : PhonopyAtoms
        Supercell
    primitive : Primitive
        Primitive cell
    displacements : array_like
        Displacements of atoms in supercell.
        shape=(num_snapshots, num_atoms, 3), dtype='double', order='C'
    forces : array_like
        Forces of atoms in supercell.
        shape=(num_snapshots, num_atoms, 3), dtype='double', order='C'
    fc_calculator : str, optional
        Currently only 'alm' is supported. Default is None, meaning invoking
        'alm'.
    atom_list : array_like of int or None, optional
        List of supercell atomic indices that represent the first indices of
        force constant matrix. The default is None, which means all atoms in
        supercell. Two shapes of force constant matrix are readable in phonopy,
        atom_list == [all atomic indices in supercell] (full) or atom_list ==
        [all atoms in primitive cell in supercell atomic indices] (compact).
        The later is the same as primitive.p2s_map.
    fc_options : str, optional
        This is arbitrary string.
    log_level : integer or bool, optional
        Verbosity level. False or 0 means quiet. True or 1 means normal level
        of log to stdout. 2 gives verbose mode.

    Returns
    -------
    fc2 : ndarray
        2nd order force constants.
        shape=(len(atom_list), num_atoms, 3, 3), dtype='double', order='C'.

    """

    if fc_calculator == 'alm' or fc_calculator is None:
        return get_fc2_alm(supercell,
                           primitive,
                           displacements,
                           forces,
                           atom_list=atom_list,
                           alm_options=fc_options,
                           log_level=log_level)
