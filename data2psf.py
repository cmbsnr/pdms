# -*- coding: utf-8 -*-
"""
Created on Thu Mar 13 09:31:53 2025

@author: cmbsn
"""

#!/usr/bin/env python3
import sys

def parse_lammps_data(filename):
    """
    解析 LAMMPS data 文件，提取 Masses、Atoms 和 Bonds 信息
    """
    atoms = []
    bonds = []
    masses = {}
    section = None
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            # 忽略空行或注释行
            if not line or line.startswith("#"):
                continue
            # 检查段落标识
            if line.startswith("Masses"):
                section = "Masses"
                continue
            elif line.startswith("Atoms"):
                section = "Atoms"
                continue
            elif line.startswith("Bonds"):
                section = "Bonds"
                continue
            
            # 根据当前段落进行处理
            if section == "Masses":
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        atom_type = int(parts[0])
                        mass = float(parts[1])
                        masses[atom_type] = mass
                    except ValueError:
                        continue
            elif section == "Atoms":
                parts = line.split()
                # 假设 LAMMPS Atoms 段格式为：
                # atom-ID molecule-ID atom-type charge x y z
                if len(parts) >= 7:
                    try:
                        atom_id = int(parts[0])
                        molecule_id = int(parts[1])
                        atom_type = int(parts[2])
                        charge = float(parts[3])
                        x = float(parts[4])
                        y = float(parts[5])
                        z = float(parts[6])
                        atoms.append({
                            "id": atom_id,
                            "molecule": molecule_id,
                            "type": atom_type,
                            "charge": charge,
                            "x": x,
                            "y": y,
                            "z": z
                        })
                    except ValueError:
                        continue
            elif section == "Bonds":
                parts = line.split()
                # 假设 Bonds 段格式为：
                # bond-ID bond-type atom1 atom2
                if len(parts) >= 4:
                    try:
                        # 忽略 bond-ID 和 bond-type，只保留原子序号
                        atom1 = int(parts[2])
                        atom2 = int(parts[3])
                        bonds.append((atom1, atom2))
                    except ValueError:
                        continue
    return atoms, bonds, masses

def write_psf(atoms, bonds, masses, output_filename):
    """
    将解析得到的原子和键信息写入 PSF 文件
    PSF 文件格式示例：
      PSF
     
             1 !NTITLE
       REMARKS Converted from LAMMPS data file
     
            X !NATOM
       [atom lines...]
     
            Y !NBOND: bonds
       [键信息，每行最多 4 对]
    """
    with open(output_filename, 'w') as f:
        # 写入 PSF 文件头部信息
        f.write("PSF\n")
        
        # 写入原子数及原子信息
        f.write("{:>8} !NATOM\n".format(len(atoms)))
        for atom in atoms:
            # PSF 原子行通常包含：原子序号、分段名、残基号、残基名、原子名、原子类型、charge、坐标(x,y,z)
            # 这里统一使用默认值：分段名 "SYS"，残基号 1，残基名 "RES"，原子名按原子序号生成
            atom_id = atom["id"]
            segname = "1"
            resid = "H"
            resname = "H"
            charge = atom["charge"]
            # 注意字段宽度及对齐方式可根据实际需求调整
            f.write("{:>8}{:>7}{:>10}{:>5}      {:>8.6f}        1.0079\n".format(
                atom_id, segname, resid, resname, charge
            ))
        
        f.write("\n")
        # 写入键（Bonds）信息
        f.write("{:>8} !NBOND: bonds\n".format(len(bonds)))
        # 每行输出最多 4 对键（共8个数字）
        line_buf = ""
        count = 0
        for (a1, a2) in bonds:
            line_buf += "{:>8}{:>8}".format(a1, a2)
            count += 1
            if count % 4 == 0:
                f.write(line_buf + "\n")
                line_buf = ""
        if line_buf:
            f.write(line_buf + "\n")
        # 如果需要添加 ANGLES、DIHEDRALS 等信息，可在此处扩展

def convert_lammps_to_psf(lammps_file, psf_file):
    atoms, bonds, masses = parse_lammps_data(lammps_file)
    write_psf(atoms, bonds, masses, psf_file)
    
if __name__ == "__main__":
    lammps_file = 'PDMS_50.data' if len(sys.argv) < 2 else sys.argv[1]
    psf_file = lammps_file[:-4] + 'psf'
    convert_lammps_to_psf(lammps_file, psf_file)
    print("转换完成！")
