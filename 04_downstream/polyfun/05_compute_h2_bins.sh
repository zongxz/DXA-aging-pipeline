#!/bin/bash
# ---------------------------------------------------------------------------
# DXA-aging-pipeline / 04_downstream
# 第 5 步：分箱 h2 回归，输出 SNPVAR
# 注意: 文中 /path/to/... 为占位符，运行前请按本机实际路径替换（见 config/）。
# ---------------------------------------------------------------------------
cd /path/to/project/GWAS/PolyFun

mkdir -p h2 logs

for chr in {1..22}; do
    ln -sf \
        ../ldscores/double_right.${chr}.bins.parquet \
        h2/double_right.${chr}.bins.parquet

    ln -sf \
        ../ldscores/double_right.${chr}.l2.ldscore.parquet \
        h2/double_right.${chr}.l2.ldscore.parquet
        
    ln -sf \
        ../ldscores/double_right.${chr}.l2.M \
        h2/double_right.${chr}.l2.M
done


set -euo pipefail

cd /path/to/project/GWAS/PolyFun

PYTHON="python"
POLYFUN="software/polyfun/polyfun.py"

SUMSTATS="munge/double_right.parquet"
OUTPREFIX="h2/double_right"
W_LD_CHR="data/baselineLF2.2.UKB/weights.UKB."

LOGDIR="logs"
LOG="${LOGDIR}/double_right_compute_h2_bins.log"

mkdir -p h2 "$LOGDIR"

echo
echo "============================================================"
echo "PolyFun Stage 4: compute-h2-bins"
echo "============================================================"
echo "[INFO] PolyFun   : $POLYFUN"
echo "[INFO] Sumstats  : $SUMSTATS"
echo "[INFO] Output    : $OUTPREFIX"
echo "[INFO] W-LD-CHR  : $W_LD_CHR"
echo "[INFO] Log       : $LOG"
echo

# ------------------------------------------------------------
# 1. Check summary statistics
# ------------------------------------------------------------

if [ ! -s "$SUMSTATS" ]; then
    echo "[ERROR] Summary statistics not found:"
    echo "        $SUMSTATS"
    exit 1
fi

echo "[CHECK] Summary statistics: OK"

# ------------------------------------------------------------
# 2. Check PolyFun LD scores and bins
# ------------------------------------------------------------

echo
echo "[CHECK] Checking PolyFun LD scores and bins..."

for chr in {1..22}; do

    BINS="h2/double_right.${chr}.bins.parquet"
    LDSCORE="h2/double_right.${chr}.l2.ldscore.parquet"

    if [ ! -e "$BINS" ]; then
        echo "[ERROR] Missing bins:"
        echo "        $BINS"
        exit 1
    fi

    if [ ! -e "$LDSCORE" ]; then
        echo "[ERROR] Missing LD score:"
        echo "        $LDSCORE"
        exit 1
    fi

    echo "CHR${chr}: OK"

done

echo "[CHECK] PolyFun LD scores and bins CHR1-22: OK"

# ------------------------------------------------------------
# 3. Check S-LDSC weights
# ------------------------------------------------------------

echo
echo "[CHECK] Checking S-LDSC weights..."

for chr in {1..22}; do

    FILE="${W_LD_CHR}${chr}.l2.ldscore.parquet"

    if [ ! -s "$FILE" ]; then
        echo "[ERROR] Missing S-LDSC weight:"
        echo "        $FILE"
        exit 1
    fi

    echo "CHR${chr}: OK"

done

echo "[CHECK] S-LDSC weights CHR1-22: OK"

# ------------------------------------------------------------
# 4. Run PolyFun Stage 4
# ------------------------------------------------------------

echo
echo "============================================================"
echo "[RUN] PolyFun --compute-h2-bins"
echo "============================================================"
echo

$PYTHON "$POLYFUN" \
    --compute-h2-bins \
    --allow-missing \
    --output-prefix "$OUTPREFIX" \
    --sumstats "$SUMSTATS" \
    --w-ld-chr "$W_LD_CHR" \
    2>&1 | tee "$LOG"

# ------------------------------------------------------------
# 5. Finish
# ------------------------------------------------------------

echo
echo "============================================================"
echo "[DONE] PolyFun Stage 4 completed"
echo "============================================================"

echo
echo "[OUTPUT FILES]"

find h2 -maxdepth 1 -type f \
    -name "double_right*" \
    -printf '%f\n' | sort

echo
echo "[LOG]"
echo "$LOG"