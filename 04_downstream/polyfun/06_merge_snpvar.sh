#!/bin/bash
# ---------------------------------------------------------------------------
# DXA-aging-pipeline / 04_downstream
# 第 6 步：合并 22 条染色体的 SNPVAR
# 注意: 文中 /path/to/... 为占位符，运行前请按本机实际路径替换（见 config/）。
# ---------------------------------------------------------------------------
set -euo pipefail

cd /path/to/project/GWAS/PolyFun

OUT="finemap/double_right_snpvar.tsv.gz"
TMP="finemap/double_right_snpvar.tmp.tsv"

mkdir -p finemap

echo "============================================================"
echo "PolyFun Stage 5: Prepare SNPVAR summary statistics"
echo "============================================================"

rm -f "$TMP" "$OUT"

first=1

for chr in {1..22}; do

    FILE="h2/double_right.${chr}.snpvar_constrained.gz"

    if [ ! -s "$FILE" ]; then
        echo "[ERROR] Missing SNPVAR:"
        echo "        $FILE"
        exit 1
    fi

    echo "[CHR${chr}] $FILE"

    if [ "$first" -eq 1 ]; then
        zcat "$FILE" > "$TMP"
        first=0
    else
        zcat "$FILE" | tail -n +2 >> "$TMP"
    fi

done

echo
echo "[CHECK] SNPVAR rows:"
wc -l "$TMP"

echo
echo "[CHECK] Header:"
head -1 "$TMP"

echo
echo "[CHECK] First 3 rows:"
head -4 "$TMP"

echo
echo "[CHECK] Compressing..."

gzip -c "$TMP" > "$OUT"

rm -f "$TMP"

echo
echo "[DONE]"
ls -lh "$OUT"

echo
echo "============================================================"
echo "Stage 5 completed"
echo "============================================================"