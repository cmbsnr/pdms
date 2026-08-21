#!/bin/bash
# 用法: ./run.sh <input_name> [np]
# 自动计算 PDMS 表面位置、步数，传给 LAMMPS
set -e

NAME="${1:?Usage: $0 <input_name> [np]}"
NP="${2:-8}"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
DATA_FILE="${NAME}.data"

[ -f "$DATA_FILE" ] || { echo "错误: 找不到 $DATA_FILE"; exit 1; }

# 参数
INDENT_SPEED=0.25      # Å/ps
DT=0.002               # ps
INIT_GAP=20.0          # 压头上方初始间距 Å
INDENT_DEPTH=20.0      # 压入深度 Å

# 计算 PDMS 表面 z（top 1% type-1 min z）
SURFACE_Z=$(python3 "${SCRIPT_DIR}/top5_z.py" --surface "$DATA_FILE")
echo "PDMS surface (top1% type-1 min z) = $SURFACE_Z"

# 计算步数: total_travel / (speed * dt) = (init_gap + indent_depth) / (0.25 * 0.002)
PER_STEP=$(python3 -c "print(${INDENT_SPEED} * ${DT})")
TOTAL_TRAVEL=$(python3 -c "print(${INIT_GAP} + ${INDENT_DEPTH})")
NSTEPS_INDENT=$(python3 -c "print(int((${INIT_GAP} + ${INDENT_DEPTH}) / (${INDENT_SPEED} * ${DT})))")
NSTEPS_RETRACT=$((NSTEPS_INDENT * 2))

echo "total_travel = ${TOTAL_TRAVEL} Å, per_step = ${PER_STEP} Å"
echo "nsteps_indent = ${NSTEPS_INDENT}, nsteps_retract = ${NSTEPS_RETRACT}"

mpirun -np "$NP" lmp_mpi \
    -in "${SCRIPT_DIR}/in.adhesion_adv" \
    -var input_name "$NAME" \
    -var surface_z "$SURFACE_Z" \
    -var total_travel "$TOTAL_TRAVEL" \
    -var nsteps_indent "$NSTEPS_INDENT" \
    -var nsteps_retract "$NSTEPS_RETRACT"
