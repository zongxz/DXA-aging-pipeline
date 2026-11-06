#!/bin/bash
# ---------------------------------------------------------------------------
# DXA-aging-pipeline / 04_downstream
# 第 4 步：逐染色体计算 LD score
# 注意: 文中 /path/to/... 为占位符，运行前请按本机实际路径替换（见 config/）。
# ---------------------------------------------------------------------------

set -euo pipefail

cd /path/to/project/GWAS/PolyFun

SUMSTATS_PREFIX="ldscores/double_right"
BFILE_CHR="/path/to/project/GWAS/PolyFun/genotype_cm"
LOGDIR="logs"

mkdir -p "$LOGDIR"

for chr in {1..22}; do

    OUT="${SUMSTATS_PREFIX}.${chr}.l2.ldscore.parquet"
    LOG="${LOGDIR}/double_right_ldscore_chr${chr}_genetic_map.log"

    echo
    echo "============================================================"
    echo "CHR${chr}"
    echo "============================================================"

    # 已经成功完成的染色体直接跳过
    if [ -s "$OUT" ]; then
        echo "[SKIP] CHR${chr} already completed:"
        echo "       $OUT"
        continue
    fi

    echo "[RUN] CHR${chr}"
    echo "[LOG] $LOG"

    python software/polyfun/polyfun.py \
        --compute-ldscores \
        --allow-missing \
        --output-prefix "$SUMSTATS_PREFIX" \
        --bfile-chr "$BFILE_CHR/" \
        --chr "$chr" \
        2>&1 | tee "$LOG"

    # 如果 python 返回成功，检查输出是否存在
    if [ ! -s "$OUT" ]; then
        echo
        echo "[ERROR] CHR${chr} finished but output file was not found:"
        echo "        $OUT"
        exit 1
    fi

    echo "[DONE] CHR${chr}"
done

echo
echo "============================================================"
echo "ALL CHROMOSOMES COMPLETED"
echo "============================================================"