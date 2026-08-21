import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import glob
import re
import subprocess

base_dir = os.path.dirname(os.path.abspath(__file__))

# Find all matching subdirectories and sort by the N2 value
pattern = os.path.join(base_dir, "pdms_DFS_*_adhension_eq")
dirs = sorted(glob.glob(pattern))

# Color cycle for distinct curves
colors = plt.cm.tab10.colors

fig, ax = plt.subplots(figsize=(10, 6))

for i, d in enumerate(dirs):
    csv_path = os.path.join(d, "force_10_100_1000.csv")
    log_path = os.path.join(d, "log.adhesion_adv")
    if not os.path.isfile(csv_path) or not os.path.isfile(log_path):
        continue

    # Extract label: the number before "N2" (e.g. "960" from "...n1_960N2...")
    folder_name = os.path.basename(d)
    match = re.search(r'n1_(\d+)N2', folder_name)
    label = match.group(1) if match else os.path.basename(d)

    # Extract Lx from log.adhesion_adv (first data line after "Step ... Lx" header)
    lx = None
    with open(log_path, 'r') as f:
        found_header = False
        for line in f:
            if not found_header:
                if line.strip().startswith('Step') and 'Lx' in line:
                    found_header = True
            else:
                stripped = line.strip()
                if stripped and stripped[0].isdigit():
                    parts = line.split()
                    if len(parts) >= 9:
                        lx = float(parts[8])  # 9th column (0-indexed: 8)
                    break

    if lx is None:
        print(f"WARNING: Could not extract Lx from {log_path}, skipping")
        continue

    # Read CSV
    depths, forces = [], []
    with open(csv_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= 3:
                depths.append(float(parts[1]))
                forces.append(float(parts[2]) * 10 / lx)

    color = colors[i % len(colors)]
    ax.plot(depths, forces, marker='o', markersize=2, linewidth=1.2,
            color=color, label=label)

ax.set_xlabel("Depth (Å)")
ax.set_ylabel("Force (N/m)")
ax.legend(title="N2", loc='best')
ax.set_title("Force per unit width vs Depth (All N2)")
fig.tight_layout()

out_path = os.path.join(base_dir, "force_vs_depth_all.png")
fig.savefig(out_path, dpi=150)
plt.close(fig)
print(f"Plotted {len(dirs)} curves -> {out_path}")
