#!/bin/bash
# ---------------------------------------------------------------------------
# DXA-aging-pipeline / 04_downstream
# 第 7 步：按 3 Mb 窗口生成精细定位任务清单
# 注意: 文中 /path/to/... 为占位符，运行前请按本机实际路径替换（见 config/）。
# ---------------------------------------------------------------------------
set -euo pipefail

cd /path/to/project/GWAS/PolyFun

PYTHON="python"
SUMSTATS="finemap/double_right_snpvar.tsv.gz"
REGIONS="software/polyfun/ukb_regions.tsv.gz"

OUTDIR="finemap/jobs"

mkdir -p "$OUTDIR"

echo "============================================================"
echo "PolyFun Stage 6: Create local-LD fine-mapping jobs"
echo "============================================================"

if [ ! -s "$SUMSTATS" ]; then
    echo "[ERROR] Missing:"
    echo "$SUMSTATS"
    exit 1
fi

if [ ! -s "$REGIONS" ]; then
    echo "[ERROR] Missing:"
    echo "$REGIONS"
    exit 1
fi

for chr in {1..22}; do

    OUT="${OUTDIR}/double_right_chr${chr}.jobs"

    echo
    echo "------------------------------------------------------------"
    echo "CHR${chr}"
    echo "------------------------------------------------------------"

    rm -f "$OUT"

    "$PYTHON" software/polyfun/create_finemapper_jobs.py \
        --sumstats "$SUMSTATS" \
        --n 48850 \
        --method susie \
        --max-num-causal 2 \
        --chr "$chr" \
        --out-prefix finemap/double_right \
        --jobs-file "$OUT" \
        --regions-file "$REGIONS" \
        --python3 "$PYTHON"

    echo "[INFO] Created:"
    echo "       $OUT"

    echo "[INFO] Number of jobs:"
    wc -l "$OUT"

done

echo
echo "============================================================"
echo "Stage 6 completed"
echo "============================================================"

echo
echo "[JOB SUMMARY]"

for chr in {1..22}; do
    FILE="${OUTDIR}/double_right_chr${chr}.jobs"
    printf "CHR%-2s : " "$chr"
    wc -l < "$FILE"
done