#!/usr/bin/env python3
"""
表面张力数据分析脚本
=====================
分析 Martini 粗粒化水的 LAMMPS 表面张力模拟结果，
与实验值（300K 下 71.66 mN/m）进行比较。

表面张力计算公式（slab 法）：
  γ = 0.5 * Lz * [P_zz - 0.5*(P_xx + P_yy)]

单位换算：1 atm·Å = 0.0101325 mN/m

读取文件：
  - surface_tension_ts.txt : 表面张力时间序列
  - log.lammps              : LAMMPS 热力学输出
  - density_z.txt           : z 方向密度分布
"""

import sys
import os
import math
import subprocess

# 实验参考值
EXP_GAMMA = 71.66  # mN/m at 300K
EXP_TEMP = 300.0   # K

# 颜色代码（终端输出）
GREEN = '\033[92m'
YELLOW = '\033[93m'
RED = '\033[91m'
CYAN = '\033[96m'
BOLD = '\033[1m'
RESET = '\033[0m'

def read_surface_tension(filename="surface_tension_ts.txt"):
    """
    读取表面张力时间序列文件。
    格式：Step gamma(mN/m) gamma_raw(atm*A) P_lateral(atm) P_normal(atm)
    """
    data = []
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            cols = line.split()
            if len(cols) >= 5:
                try:
                    step = int(cols[0])
                    gamma = float(cols[1])
                    gamma_raw = float(cols[2])
                    p_lat = float(cols[3])
                    p_norm = float(cols[4])
                    data.append({
                        'step': step,
                        'gamma': gamma,
                        'gamma_raw': gamma_raw,
                        'p_lateral': p_lat,
                        'p_normal': p_norm
                    })
                except ValueError:
                    continue
    return data

def read_thermo(filename="log.lammps"):
    """
    读取 LAMMPS 热力学输出文件。
    """
    data = []
    with open(filename, 'r') as f:
        lines = f.readlines()

    # 找到数据起始行
    in_data = False
    headers = []
    for line in lines:
        line_stripped = line.strip()
        if 'Step' in line_stripped and 'Temp' in line_stripped:
            headers = line_stripped.split()
            in_data = True
            continue
        if in_data and line_stripped:
            cols = line_stripped.split()
            if len(cols) >= len(headers):
                row = {}
                for i, h in enumerate(headers):
                    try:
                        row[h] = float(cols[i])
                    except ValueError:
                        row[h] = cols[i]
                data.append(row)
    return data

def read_density_profile(filename="density_z.txt"):
    """
    读取沿 z 方向的密度分布。
    格式（LAMMPS ave/chunk 输出）：
      # 头部注释（以 # 开头）
      Timestep Number-of-chunks Total-count
      Chunk Coord1 Ncount density/number
    注意：Ncount 可能是浮点数！
    """
    data = []
    with open(filename, 'r') as f:
        lines = f.readlines()

    for i, line in enumerate(lines):
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        cols = line.split()
        # 跳过头部行（Timestep Number-of-chunks Total-count）
        if len(cols) == 3 and cols[0].isdigit() and not cols[1].startswith('-') and not cols[1].__contains__('.'):
            continue
        # 数据行：Chunk Coord1 Ncount density/number
        if len(cols) >= 4:
            try:
                chunk = int(cols[0])
                coord = float(cols[1])
                ncount = float(cols[2])  # Ncount 可能是浮点数
                density = float(cols[3])
                data.append({
                    'chunk': chunk,
                    'ncount': ncount,
                    'density': density,
                    'coord': coord
                })
            except (ValueError, IndexError):
                continue
    return data

def analyze(data, thermo_data, density_data, lz):
    """
    分析表面张力数据，计算统计量。
    """
    if not data:
        print(f"{RED}错误：表面张力数据为空！{RESET}")
        return

    # 使用后半段数据（平衡后的数据更可靠）
    n = len(data)
    n_equil = n // 2  # 前半段作为平衡，后半段用于分析
    equil_data = data[n_equil:]

    gammas = [d['gamma'] for d in equil_data]
    n_samples = len(gammas)

    if n_samples < 2:
        print(f"{RED}错误：平衡后样本数不足！{RESET}")
        return

    mean_gamma = sum(gammas) / n_samples

    # 样本标准差
    variance = sum((g - mean_gamma)**2 for g in gammas) / (n_samples - 1)
    std_gamma = math.sqrt(variance)

    # 标准误差（考虑到时间序列的自相关，使用 block averaging）
    # 简单估计：分成 5 块
    n_blocks = min(5, n_samples // 10)
    if n_blocks >= 2:
        block_size = n_samples // n_blocks
        block_means = []
        for i in range(n_blocks):
            start = i * block_size
            end = start + block_size if i < n_blocks - 1 else n_samples
            block_means.append(sum(gammas[start:end]) / (end - start))
        block_mean = sum(block_means) / n_blocks
        block_var = sum((b - block_mean)**2 for b in block_means) / (n_blocks - 1)
        sem = math.sqrt(block_var / n_blocks)
    else:
        sem = std_gamma / math.sqrt(n_samples)

    # 与实验值比较
    deviation = mean_gamma - EXP_GAMMA
    deviation_pct = deviation / EXP_GAMMA * 100

    # 输出分析结果
    print(f"\n{BOLD}{'='*60}{RESET}")
    print(f"{BOLD}  Martini 粗粒化水的表面张力分析结果{RESET}")
    print(f"{BOLD}{'='*60}{RESET}")

    print(f"\n{CYAN}-- 模拟参数 --{RESET}")
    print(f"  体系大小     : {n_samples} 个采样点（平衡后）")
    print(f"  z 方向盒子   : {lz:.2f} Å")
    print(f"  温度         : {EXP_TEMP} K")

    print(f"\n{CYAN}-- 表面张力结果 --{RESET}")
    print(f"  模拟值 (Martini CG) : {mean_gamma:.4f} ± {sem:.4f} mN/m")
    print(f"  实验值 (300K)       : {EXP_GAMMA:.2f} mN/m")
    print(f"  绝对偏差            : {deviation:+.4f} mN/m")
    print(f"  相对偏差            : {deviation_pct:+.2f}%")

    # 压力张量分析
    if equil_data:
        p_lat_avg = sum(d['p_lateral'] for d in equil_data) / n_samples
        p_norm_avg = sum(d['p_normal'] for d in equil_data) / n_samples
        print(f"\n{CYAN}-- 压力张量分量 --{RESET}")
        print(f"  <P_lateral> (xy平面) : {p_lat_avg:.4f} atm")
        print(f"  <P_normal>  (z方向)  : {p_norm_avg:.4f} atm")
        print(f"  ΔP = P_n - P_l       : {p_norm_avg - p_lat_avg:.4f} atm")

    # 评价
    print(f"\n{CYAN}-- 评价 --{RESET}")
    abs_dev = abs(deviation_pct)
    if abs_dev < 5:
        print(f"  {GREEN}模拟值与实验值吻合良好（偏差 < 5%）{RESET}")
    elif abs_dev < 20:
        print(f"  {YELLOW}模拟值与实验值存在中等偏差（5% < 偏差 < 20%）{RESET}")
    else:
        print(f"  {RED}模拟值与实验值偏差较大（偏差 > 20%）{RESET}")

    print(f"\n  {CYAN}说明：Martini 标准（非极化）水模型已知会低估水的表面张力。")
    print(f"  文献中 Martini 水的表面张力通常在 30-45 mN/m 范围内。")
    print(f"  这主要是因为 Martini 水珠之间缺乏方向性氢键和极化效应。")
    print(f"  如需更准确的表面张力，可使用 Martini 极化水模型（PW）。{RESET}")

    # 密度分布分析
    if density_data:
        print(f"\n{CYAN}-- 密度分布分析 --{RESET}")
        # 找出密度 > 阈值 50% 的区域作为液相
        max_dens = max(d['density'] for d in density_data)
        half_max = max_dens / 2.0
        liquid_bins = [d for d in density_data if d['density'] > half_max]
        if liquid_bins:
            z_range = liquid_bins[-1]['coord'] - liquid_bins[0]['coord']
            bulk_dens = sum(d['density'] for d in liquid_bins) / len(liquid_bins)
            print(f"  液相区域宽度 : {z_range:.1f} Å")
            print(f"  液相体密度   : {bulk_dens:.6f} atoms/Å³")
            print(f"  最大密度     : {max_dens:.6f} atoms/Å³")

    print(f"\n{BOLD}{'='*60}{RESET}\n")

    return {
        'mean': mean_gamma,
        'sem': sem,
        'std': std_gamma,
        'deviation_pct': deviation_pct,
        'n_samples': n_samples,
        'lz': lz
    }

def try_plot(data, density_data, output_prefix="surface_tension"):
    """尝试绘制图表（需要 matplotlib）"""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print(f"{YELLOW}matplotlib 未安装，跳过绘图。{RESET}")
        return

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # 1. 表面张力时间序列
    ax1 = axes[0, 0]
    steps = [d['step'] for d in data]
    gammas = [d['gamma'] for d in data]
    ax1.plot(steps, gammas, 'b-', alpha=0.7, linewidth=0.5)
    # 滑动平均
    window = max(1, len(gammas) // 50)
    if window > 1:
        smoothed = []
        for i in range(len(gammas) - window + 1):
            smoothed.append(sum(gammas[i:i+window]) / window)
        smooth_steps = steps[window//2:window//2+len(smoothed)]
        ax1.plot(smooth_steps, smoothed, 'r-', linewidth=1.5, label=f'滑动平均 (窗口={window})')

    ax1.axhline(y=EXP_GAMMA, color='green', linestyle='--', linewidth=1.5,
                label=f'实验值 = {EXP_GAMMA} mN/m')
    # 平均值
    n_equil = len(data) // 2
    mean_g = sum(d['gamma'] for d in data[n_equil:]) / max(1, len(data) - n_equil)
    ax1.axhline(y=mean_g, color='red', linestyle=':', linewidth=1.5,
                label=f'模拟均值 = {mean_g:.2f} mN/m')
    ax1.set_xlabel('Time Step')
    ax1.set_ylabel('Surface Tension (mN/m)')
    ax1.set_title('Surface Tension vs Time')
    ax1.legend(fontsize=8)
    ax1.grid(True, alpha=0.3)

    # 2. 压力张量分量
    ax2 = axes[0, 1]
    p_lat = [d['p_lateral'] for d in data]
    p_norm = [d['p_normal'] for d in data]
    ax2.plot(steps, p_lat, 'b-', alpha=0.5, linewidth=0.5, label='P_lateral (xy)')
    ax2.plot(steps, p_norm, 'r-', alpha=0.5, linewidth=0.5, label='P_normal (z)')
    ax2.set_xlabel('Time Step')
    ax2.set_ylabel('Pressure (atm)')
    ax2.set_title('Pressure Tensor Components')
    ax2.legend(fontsize=8)
    ax2.grid(True, alpha=0.3)

    # 3. 表面张力分布直方图
    ax3 = axes[1, 0]
    equil_gammas = [d['gamma'] for d in data[len(data)//2:]]
    ax3.hist(equil_gammas, bins=30, density=True, alpha=0.7, color='steelblue', edgecolor='black')
    ax3.axvline(x=mean_g, color='red', linestyle='--', linewidth=2,
                label=f'Mean = {mean_g:.2f}')
    ax3.axvline(x=EXP_GAMMA, color='green', linestyle='--', linewidth=2,
                label=f'Exp = {EXP_GAMMA}')
    ax3.set_xlabel('Surface Tension (mN/m)')
    ax3.set_ylabel('Probability Density')
    ax3.set_title('Distribution of Surface Tension (Equilibrated)')
    ax3.legend(fontsize=8)

    # 4. 密度分布
    ax4 = axes[1, 1]
    if density_data:
        coords = [d['coord'] for d in density_data]
        dens = [d['density'] for d in density_data]
        ax4.plot(coords, dens, 'b-', linewidth=1.5)
        ax4.fill_between(coords, 0, dens, alpha=0.3, color='blue')
        ax4.set_xlabel('z Coordinate (Å)')
        ax4.set_ylabel('Number Density (atoms/Å³)')
        ax4.set_title('Density Profile along z (Slab Normal)')
        ax4.grid(True, alpha=0.3)
    else:
        ax4.text(0.5, 0.5, 'No density data', ha='center', va='center')
        ax4.set_title('Density Profile (No Data)')

    plt.tight_layout()
    outfile = f"{output_prefix}_analysis.png"
    plt.savefig(outfile, dpi=150, bbox_inches='tight')
    print(f"图表已保存至: {outfile}")
    plt.close()

def main():
    print(f"{BOLD}{CYAN}")
    print("╔══════════════════════════════════════════════════════╗")
    print("║   Martini 粗粒化水表面张力分析工具                 ║")
    print("║   参考实验值: 71.66 mN/m (300K)                    ║")
    print("╚══════════════════════════════════════════════════════╝")
    print(f"{RESET}")

    # 检查文件
    files = {
        "surface_tension_ts.txt": "表面张力时间序列",
        "log.lammps": "LAMMPS 热力学日志",
    }

    for fname, fdesc in files.items():
        if not os.path.exists(fname):
            print(f"{RED}错误：找不到文件 {fname} ({fdesc})！{RESET}")
            print("请先运行 LAMMPS 表面张力模拟。")
            sys.exit(1)

    # 读取数据
    print("正在读取表面张力时间序列...")
    data = read_surface_tension("surface_tension_ts.txt")
    print(f"  读取到 {len(data)} 个时间点")

    print("正在读取热力学日志...")
    thermo_data = read_thermo("log.lammps")
    print(f"  读取到 {len(thermo_data)} 个热力学记录")

    # 获取 Lz（盒子 z 方向长度）
    lz = None
    if thermo_data:
        lz = thermo_data[-1].get('Lz', None)
    if lz is None and data:
        # 从 raw gamma 反推 Lz
        lz = 393.5  # 默认值 3 * 131.16 ≈ 393.5
    print(f"  z 方向盒子长度 Lz = {lz:.2f} Å")

    # 读取密度分布
    density_data = []
    if os.path.exists("density_z.txt"):
        print("正在读取密度分布...")
        density_data = read_density_profile("density_z.txt")
        print(f"  读取到 {len(density_data)} 个密度区间")

    # 分析
    result = analyze(data, thermo_data, density_data, lz)

    # 绘图
    try_plot(data, density_data)

    # 输出汇总
    if result:
        with open("surface_tension_summary.txt", 'w') as f:
            f.write("="*60 + "\n")
            f.write("Martini 粗粒化水表面张力分析汇总\n")
            f.write("="*60 + "\n\n")
            f.write(f"模拟值 (Martini CG) : {result['mean']:.4f} ± {result['sem']:.4f} mN/m\n")
            f.write(f"实验值 (300K)       : {EXP_GAMMA:.2f} mN/m\n")
            f.write(f"相对偏差            : {result['deviation_pct']:+.2f}%\n")
            f.write(f"采样点数            : {result['n_samples']}\n")
            f.write(f"z 盒子长度          : {result['lz']:.2f} Å\n")
        print(f"汇总文件已保存至: surface_tension_summary.txt")

if __name__ == "__main__":
    main()
