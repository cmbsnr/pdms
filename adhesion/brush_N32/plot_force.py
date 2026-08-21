import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

# Paths
base_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(
    base_dir,
    "pdms_DFS_CB2_0.30B1_0.15B_0.02k_0.1000T300N1_32n1_960N2_0n2_0N3_0n3_0z_95NAll_30720crosslinkfinal_adhension_eq"
)
csv_path = os.path.join(data_dir, "force_10_100_1000.csv")
out_path = os.path.join(data_dir, "force_vs_depth.png")

# Read CSV manually (tab-separated, skip comment lines)
depths, forces = [], []
with open(csv_path, 'r') as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        parts = line.split()
        if len(parts) >= 3:
            depths.append(float(parts[1]))
            forces.append(float(parts[2]))

label = "960"

# Plot
fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(depths, forces, marker='o', markersize=3, linewidth=1.2, label=label)
ax.set_xlabel("Depth (Å)")
ax.set_ylabel("Force (nN)")
ax.legend()
ax.set_title(f"Force vs Depth ({label})")
fig.tight_layout()

# Save silently
fig.savefig(out_path, dpi=150)
plt.close(fig)
print(f"Plot saved to: {out_path}")
