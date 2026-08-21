#!/usr/bin/env python3
"""
分析 in.trap 的输出: 提取稳态牵引力(=摩擦力) 与 质心速度, 作图。

用法:
  python3 analyze_trap.py <name> [vd1 vd2 ...] [--save PNG]

读取 <name>/<name>_trap_vd<vd>_force.txt 和 _vel.txt:
  - force: step, c_fwat_x, v_fric_x, v_fspring_total, c_tcomw
  - vel:   step, v_vxw, v_xcmw, v_ycmw, v_zcmw, v_fspring_total
稳态(后段)中 v_fspring_total 的平均 = 该拖曳速度 vd 下的牵引力, 摩擦力 = -牵引力。
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def analyze(name, vd):
    f = f"{name}/{name}_trap_vd{vd}_force.txt"
    g = f"{name}/{name}_trap_vd{vd}_vel.txt"
    d = np.loadtxt(f, skiprows=2)
    v = np.loadtxt(g, skiprows=2)
    step, fwat, fric, fspring, T = d.T
    vstep, vxw, xcm, ycm, zcm, fsp2 = v.T
    # 稳态 = 后50%
    n = len(fspring)
    tail = slice(int(n*0.5), n)
    F_trap = fspring[tail].mean()
    F_fric = -F_trap          # 摩擦力 = -牵引力
    v_ss = vxw[tail].mean()
    T_ss = T[tail].mean()
    return dict(vd=float(vd), F_trap=F_trap, F_fric=F_fric, v_ss=v_ss, T=T_ss,
                step=step, fspring=fspring, vxw=vxw, xcm=xcm)

def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    name = args[0]
    vds = args[1:] or ["1.0e-4"]
    results = [analyze(name, vd) for vd in vds]

    print(f"{'vd (A/fs)':>10}  {'v_ss (A/fs)':>12}  {'牵引力':>10}  {'摩擦力':>10}  {'T':>7}")
    for r in results:
        print(f"{r['vd']:>10.1e}  {r['v_ss']:>12.2e}  {r['F_trap']:>10.1f}  {r['F_fric']:>10.1f}  {r['T']:>7.1f}")

    # 摩擦-速度曲线
    plt.figure(figsize=(8,6))
    plt.plot([r["v_ss"] for r in results], [r["F_fric"] for r in results],
             "o-", ms=8, lw=2)
    plt.xlabel("droplet COM velocity v (A/fs)  [1 A/fs = 1e5 m/s]")
    plt.ylabel("friction force (kcal/mol/A)")
    plt.title(f"PDMS brush friction vs velocity\n{name[:50]}...")
    plt.grid(alpha=0.3)
    out = f"{name}/friction_vs_v.png"
    plt.savefig(out, dpi=150)
    print(f"plot saved : {out}")

    # 最后一组的时程
    r = results[-1]
    plt.figure(figsize=(10,8))
    plt.subplot(3,1,1); plt.plot(r["step"], r["fspring"], "o-", ms=3)
    plt.axhline(r["F_trap"], color="g", ls="--", label=f"steady trap force={r['F_trap']:.1f}")
    plt.ylabel("trap force"); plt.legend()
    plt.subplot(3,1,2); plt.plot(r["step"], r["vxw"], "o-", ms=3)
    plt.axhline(r["v_ss"], color="g", ls=":", label=f"steady v={r['v_ss']:.2e}")
    plt.ylabel("vx (A/fs)"); plt.legend()
    plt.subplot(3,1,3); plt.plot(r["step"], r["xcm"], "o-", ms=3)
    plt.xlabel("time step"); plt.ylabel("xcm (A)")
    plt.tight_layout()
    out2 = f"{name}/trap_timeseries.png"
    plt.savefig(out2, dpi=150)
    print(f"plot saved : {out2}")

if __name__ == "__main__":
    main()
