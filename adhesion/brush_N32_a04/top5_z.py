#!/usr/bin/env python3
"""
Parse LAMMPS data files and find the minimum z-coordinate among the
top 1% of type-1 atoms (ranked by z-coordinate, highest first).

Usage:
    python top5_z.py                  # analyze the first .data file found
    python top5_z.py <data_file>      # analyze a specific file
    python top5_z.py --all            # analyze all .data files, save to CSV
"""

import sys
import os
import csv
import math
import glob


def parse_lammps_data(filepath):
    """Parse a LAMMPS data file and return atom records.

    Each record: (atom_id, molecule_id, atom_type, charge, x, y, z, nx, ny, nz)
    """
    atoms = []
    in_atoms = False

    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()

            if line.lower().startswith('atoms'):
                in_atoms = True
                continue

            if not in_atoms:
                continue

            if not line:
                continue

            if not line[0].isdigit():
                if line and not line[0].isdigit() and not line[0] in '-.':
                    first_word = line.split()[0] if line.split() else ''
                    if first_word[0].isalpha() and first_word.lower() != 'atoms':
                        break
                    continue

            parts = line.split()
            if len(parts) >= 7:
                atom_id = int(parts[0])
                molecule_id = int(parts[1])
                atom_type = int(parts[2])
                charge = float(parts[3])
                x = float(parts[4])
                y = float(parts[5])
                z = float(parts[6])
                nx = int(parts[7]) if len(parts) > 7 else 0
                ny = int(parts[8]) if len(parts) > 8 else 0
                nz = int(parts[9]) if len(parts) > 9 else 0
                atoms.append((atom_id, molecule_id, atom_type, charge, x, y, z, nx, ny, nz))

    return atoms


def analyze_file(filepath, target_type=1, top_pct=0.01):
    """Analyze a single data file and return results dict."""
    atoms = parse_lammps_data(filepath)
    n_total = len(atoms)

    z_vals = [z for (_, _, atype, _, _, _, z, _, _, _) in atoms if atype == target_type]
    n_target = len(z_vals)

    if n_target == 0:
        return {
            'filename': os.path.basename(filepath),
            'total_atoms': n_total,
            'target_atoms': 0,
            'top_n': 0,
            'min_z_top': None,
            'max_z_top': None,
            'mean_z_top': None,
        }

    z_vals.sort()
    top_n = max(1, math.ceil(n_target * top_pct))
    top_z_values = z_vals[-top_n:]

    return {
        'filename': os.path.basename(filepath),
        'total_atoms': n_total,
        'target_atoms': n_target,
        'top_n': top_n,
        'min_z_top': top_z_values[0],
        'max_z_top': top_z_values[-1],
        'mean_z_top': sum(top_z_values) / len(top_z_values),
    }


def print_result(r):
    """Pretty-print a single result."""
    print(f"\n--- {r['filename']} ---")
    print(f"Total atoms: {r['total_atoms']}")
    print(f"Type-1 atoms: {r['target_atoms']}")
    if r['min_z_top'] is not None:
        print(f"Atoms in top 1%: {r['top_n']}")
        print(f"Minimum z in top 1%: {r['min_z_top']:.6f}")
        print(f"Maximum z in top 1%: {r['max_z_top']:.6f}")
        print(f"Mean z in top 1%:    {r['mean_z_top']:.6f}")
    else:
        print("No type-1 atoms found.")


def main():
    if '--surface' in sys.argv:
        # Output only the surface z (top 1% type-1 min z) for use with LAMMPS -v
        filepath = sys.argv[sys.argv.index('--surface') + 1] if len(sys.argv) > sys.argv.index('--surface') + 1 else None
        if not filepath:
            candidates = sorted(glob.glob('*.data'))
            if not candidates:
                print("No .data file found.")
                sys.exit(1)
            filepath = candidates[0]
        r = analyze_file(filepath)
        if r['min_z_top'] is not None:
            print(f"{r['min_z_top']:.6f}")
        else:
            print("ERROR: No type-1 atoms found.", file=sys.stderr)
            sys.exit(1)
        return

    if '--all' in sys.argv:
        # Batch mode: analyze all .data files, output CSV
        data_files = sorted(glob.glob('*.data'))
        if not data_files:
            print("No .data files found.")
            sys.exit(1)

        results = []
        for f in data_files:
            print(f"Processing {f} ...")
            r = analyze_file(f)
            print_result(r)
            results.append(r)

        # Write CSV
        csv_path = 'top1_z_results.csv'
        fieldnames = ['filename', 'total_atoms', 'target_atoms',
                      'top_n', 'min_z_top', 'max_z_top', 'mean_z_top']
        with open(csv_path, 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(results)

        print(f"\nResults saved to: {csv_path}")

    else:
        # Single-file mode
        if len(sys.argv) > 1:
            filepath = sys.argv[1]
        else:
            candidates = sorted(glob.glob('*.data'))
            if not candidates:
                print("No .data file found.")
                sys.exit(1)
            filepath = candidates[0]
            print(f"Using: {filepath}")

        r = analyze_file(filepath)
        print_result(r)


if __name__ == '__main__':
    main()
