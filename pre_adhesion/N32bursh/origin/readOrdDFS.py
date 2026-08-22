#!/usr/bin/env python3
"""
Convert .ord file to LAMMPS .data file (Python port of readOrdDFS.m).

Usage:
    python readOrdDFS.py [input.ord] [output.data] [-k <value>]
    python readOrdDFS.py [input.ord] -k <value>
    python readOrdDFS.py -k <value> [input.ord] [output.data]

Default input:  CB2_0.30B1_0.15B_0.02k_0.1000T300N1_32n1_2880N2_0n2_0N3_0n3_0z_95NAll_92160crosslinkfinal.ord
Default output: pdms_DFS_<input_basename>.data
Default k:      4 (distance scaling factor)
"""

import sys
import os
import argparse
import numpy as np
from collections import defaultdict


def dfs_paths(adj, max_len):
    """
    Find all simple paths in graph `adj` up to length `max_len + 1`.
    Returns list of paths (each path is a list of node indices).

    Equivalent to MATLAB dfs_paths(): for maxLen=2, returns paths of
    length 2 and 3; for maxLen=3, returns paths of length 2, 3, and 4.
    """
    N = len(adj)  # adj is 1-indexed (index 0 unused or used as placeholder)
    paths = []

    for start in range(1, N):  # 1-indexed, skip dummy index 0 if present
        stack = [[start]]

        while stack:
            path = stack.pop()

            if len(path) >= 2:
                paths.append(path)

            if len(path) == max_len + 1:
                continue

            last = path[-1]
            for nei in adj[last]:
                if nei not in path:  # simple path
                    stack.append(path + [nei])

    return paths


def lexicographically_smaller(a, b):
    """Return True if list `a` is lexicographically smaller than list `b`."""
    for ai, bi in zip(a, b):
        if ai < bi:
            return True
        elif ai > bi:
            return False
    return False  # equal


def read_ord_to_data(ord_filename, output_filename=None, k=4, r=0.2):
    """
    Read .ord file and write LAMMPS .data file.

    Parameters
    ----------
    ord_filename : str
        Path to input .ord file.
    output_filename : str or None
        Path to output .data file. If None, auto-generate.
    k : float
        Distance scaling factor (default 4).
    r : float
        Random spherical perturbation radius (default 0.2).
    """
    if output_filename is None:
        base = os.path.splitext(os.path.basename(ord_filename))[0]
        output_filename = f'pdms_DFS_{base}.data'

    print(f'Reading: {ord_filename}')
    print(f'Output:  {output_filename}')

    # Read the .ord file.
    # The first line is a header with fewer columns than the data rows.
    # We read all lines as text, skip the header, and parse the rest.
    with open(ord_filename, 'r') as f:
        lines = f.readlines()

    # Parse header line to get Natom (may not match data line count)
    header_parts = lines[0].strip().split()
    print(f'Header: {header_parts}')

    # Parse data lines (skip header)
    data_rows = []
    for line in lines[1:]:
        line = line.strip()
        if not line:
            continue
        parts = line.split()
        data_rows.append([float(x) for x in parts])

    data = np.array(data_rows)
    Natom, column = data.shape
    print(f'Data matrix: {Natom} atoms x {column} columns')

    Tatom = 1
    group = int((column - 5) / 2)
    print(f'group = {group}')

    # --- Coordinate processing ---
    xyz = data[:, 0:3].copy()

    # Random spherical perturbation
    rand_vec = np.random.randn(Natom, 3)
    norms = np.sqrt(np.sum(rand_vec ** 2, axis=1, keepdims=True))
    unit_vec = rand_vec / norms
    displacement = r * unit_vec
    xyz = xyz + displacement

    # Scale coordinates
    xyz = k * xyz

    # Box dimensions
    xlo = 0.0
    xhi = k * 300
    ylo = 0.0
    yhi = k * 32
    zlo = np.min(xyz[:, 2])
    zhi = np.max(xyz[:, 2])

    # --- Atom data ---
    atomID = np.arange(1, Natom + 1, dtype=int)
    atomType = np.ones(Natom, dtype=int)
    atomCharge = np.zeros(Natom, dtype=int)
    molemerID = data[:, -1].astype(int)

    Atoms = np.column_stack([
        atomID,
        atomType,
        atomType,
        atomCharge,
        xyz,
    ])

    # --- Bond extraction ---
    # Connectivity is in columns: 5 + group to 5 + 2*group - 1 (0-indexed)
    bondList = []
    bondID = 1

    for i in range(Natom):
        # Column indices for connectivity (0-indexed)
        conn_start = 5 + group - 1  # 0-indexed: column index
        for g_idx in range(group):
            col = conn_start + g_idx
            j_val = int(data[i, col])
            j = j_val  # 1-indexed atom reference
            if j > 0 and j > (i + 1):  # prevent duplicates (i+1 because 1-indexed)
                bondList.append([bondID, 1, i + 1, j])
                bondID += 1

    bondList = np.array(bondList) if bondList else np.empty((0, 4), dtype=int)
    Nbond = bondList.shape[0]
    print(f'Bonds: {Nbond}')

    # --- Build adjacency list ---
    adj = defaultdict(list)
    # Use 1-indexed adjacency
    for b in range(Nbond):
        a_i = int(bondList[b, 2])
        a_j = int(bondList[b, 3])
        adj[a_i].append(a_j)
        adj[a_j].append(a_i)

    # Convert to list-of-lists for 1-indexed access
    max_node = max(adj.keys()) if adj else 0
    adj_list = [[] for _ in range(max_node + 1)]
    for node, neighbors in adj.items():
        adj_list[node] = neighbors

    print(f'Adjacency built: {max_node} nodes')

    # --- Angle extraction via DFS ---
    print('Finding angles (DFS depth=2)...')
    paths2 = dfs_paths(adj_list, 2)

    angleSet = []
    for p in paths2:
        if len(p) == 3:
            if p[0] < p[2]:  # dedup: keep i < k
                angleSet.append(p)

    if angleSet:
        angleSet = np.unique(np.array(angleSet), axis=0)
    else:
        angleSet = np.empty((0, 3), dtype=int)

    Nangle = angleSet.shape[0]
    print(f'Angles: {Nangle}')

    if Nangle > 0:
        angleList = np.column_stack([
            np.arange(1, Nangle + 1),
            np.ones(Nangle, dtype=int),
            angleSet,
        ])
    else:
        angleList = np.empty((0, 5), dtype=int)

    # --- Dihedral extraction via DFS ---
    print('Finding dihedrals (DFS depth=3)...')
    paths3 = dfs_paths(adj_list, 3)

    dihedralSet = []
    for p in paths3:
        if len(p) == 4:
            p_rev = list(reversed(p))
            if lexicographically_smaller(p, p_rev):
                dihedralSet.append(p)

    if dihedralSet:
        dihedralSet = np.unique(np.array(dihedralSet), axis=0)
    else:
        dihedralSet = np.empty((0, 4), dtype=int)

    Ndihedral = dihedralSet.shape[0]
    print(f'Dihedrals: {Ndihedral}')

    if Ndihedral > 0:
        dihedralList = np.column_stack([
            np.arange(1, Ndihedral + 1),
            np.ones(Ndihedral, dtype=int),
            dihedralSet,
        ])
    else:
        dihedralList = np.empty((0, 6), dtype=int)

    # --- Write LAMMPS data file ---
    print(f'Writing {output_filename}...')
    with open(output_filename, 'w') as f:
        # Header
        f.write('# Generated by Python (readOrdDFS.py)\n')
        f.write(f'{Tatom} atom types\n')
        f.write(f'{Natom} atoms\n')
        f.write(f'1 bond types\n')
        f.write(f'{Nbond} bonds\n')
        f.write(f'1 angle types\n')
        f.write(f'{Nangle} angles\n')
        f.write(f'1 dihedral types\n')
        f.write(f'{Ndihedral} dihedrals\n\n')

        # Box
        f.write(f'{xlo:.4f} {xhi:.4f} xlo xhi\n')
        f.write(f'{ylo:.4f} {yhi:.4f} ylo yhi\n')
        f.write(f'{zlo:.4f} {zhi:.4f} zlo zhi\n\n')

        # Masses
        f.write('Masses\n\n')
        f.write('1 72\n\n')

        # Atoms
        f.write('Atoms  # full\n\n')
        np.savetxt(f, Atoms, fmt='%d %d %d %d %.6f %.6f %.6f')

        # Bonds
        f.write('\nBonds\n\n')
        if Nbond > 0:
            np.savetxt(f, bondList, fmt='%d %d %d %d')

        # Angles
        f.write('\nAngles\n\n')
        if Nangle > 0:
            np.savetxt(f, angleList, fmt='%d %d %d %d %d')

        # Dihedrals
        f.write('\nDihedrals\n\n')
        if Ndihedral > 0:
            np.savetxt(f, dihedralList, fmt='%d %d %d %d %d %d')

        f.write('\n')

    print('Done!')
    return output_filename


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Convert .ord file to LAMMPS .data file'
    )
    parser.add_argument(
        'ord_file', nargs='?',
        default='CB2_0.30B1_0.15B_0.02k_0.1000T300N1_32n1_2880N2_0n2_0N3_0n3_0z_95NAll_92160crosslinkfinal.ord',
        help='Input .ord file'
    )
    parser.add_argument(
        'out_file', nargs='?', default=None,
        help='Output .data file (auto-generated if omitted)'
    )
    parser.add_argument(
        '-k', type=float, default=4.0,
        help='Distance scaling factor (default: 4)'
    )

    args = parser.parse_args()
    read_ord_to_data(args.ord_file, args.out_file, k=args.k)
