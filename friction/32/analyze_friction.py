#!/usr/bin/env python3
"""
分析 in.friction 的输出 (恒力驱动法)。

用法:
  python3 analyze_friction.py <name> [--drive-start-step N] [--save PNG]
其中 <name> 是刷层data文件前缀 (同 in.friction 的 -var name)。

说明 (本体系实测, 速度弱化摩擦):
  * 恒力法测不出"终端速度"。若 Fdrive < 静摩擦(~90-100), 水珠全程钉扎,
    vx 始终≈0 -> 说明 Fdrive 不足以脱钉;
  * 若 Fdrive > 静摩擦, 水珠脱钉后持续加速, vx 单调增大 (且刷层升温),
    不存在稳态。此时 vx 首次显著非零的时间 对应脱钉。
  * 因此该脚本只画 vx/xcm/力 的时程, 供判断"钉扎 vs 滑动"和脱钉时刻,
    不再计算终端速度。
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    name = args[0] if args else "pdms_DFS_CB2_0.30B1_0.15B_0.02k_0.1000T300N1_32n1_960N2_0n2_0N3_0n3_0z_95NAll_30720crosslinkfinal_brush_out"

    f = f"{name}/{name}_friction_vel.txt"
    g = f"{name}/{name}_friction_force.txt"
    d = np.loadtxt(f, skiprows=2)
    step, vx, vy, vz, xcm, ycm, zcm = d.T
    df = np.loadtxt(g, skiprows=2)
    fstep, fwat, fric, T = df.T

    # 判断钉扎/滑动: vx 首次显著偏离0的时间 (超过 ~5e-5 A/fs 且持续)
    t0 = 0
    for i, s in enumerate(step):
        if i > 0 and np.abs(vx[i]) > 5e-5:
            t0 = i
            break
    # 末段是否持续加速 (若无终端速度, vx 末段仍在增长)
    accel = "持续加速(无终端速度)" if len(vx) > 10 and vx[-1] > vx[-5] + 1e-5 else "接近稳态/钉扎"
    # 摩擦阻力范围 (忽略初期瞬态)
    fric_mid = fric[len(fric)//2:]

    print(f"file       : {f}")
    print(f"steps      : {step[0]} .. {step[-1]}")
    print(f"depin time : step {step[t0]} (vx 首次超过5e-5 A/fs)")
    print(f"final vx   : {vx[-1]:.3e} A/fs  = {vx[-1]*1e5:.0f} m/s  ({accel})")
    print(f"friction   : 末段 {fric_mid.min():+.1f} ~ {fric_mid.max():+.1f} kcal/mol/A"
          f" (驱动开始附近≈静摩擦, 见 thermo/force 早期 fric≈-{abs(fric[0]) if len(fric)>1 else '?'})")
    print(f"note       : 若 vx 基本不动 -> 钉扎, Fdrive < 静摩擦; "
          f"若 vx 单调增长 -> 已脱钉, 无终端速度, 该法只能测静摩擦阈值")

    plt.figure(figsize=(10, 8))
    plt.subplot(3, 1, 1)
    plt.plot(step, vx, "o-", ms=3, label="vx")
    plt.axvline(step[t0], color="r", ls="--", label="depin start")
    plt.ylabel("vx (A/fs)"); plt.legend()
    plt.subplot(3, 1, 2)
    plt.plot(step, xcm, "o-", ms=3)
    plt.axvline(step[t0], color="r", ls="--")
    plt.ylabel("xcm (A)")
    plt.subplot(3, 1, 3)
    plt.plot(fstep, fric, "o-", ms=3, color="C1")
    plt.axvline(step[t0], color="r", ls="--")
    plt.ylabel("friction force (kcal/mol/A)"); plt.xlabel("time step (1 fs)")
    plt.tight_layout()
    out = f"{name}/friction_vx.png"
    plt.savefig(out, dpi=150)
    print(f"plot saved : {out}")

if __name__ == "__main__":
    main()
