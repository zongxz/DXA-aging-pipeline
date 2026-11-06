#!/bin/bash
# ---------------------------------------------------------------------------
# DXA-aging-pipeline / 04_downstream
# 第 3 步：构建带真实遗传距离(cM)的去重基因型（关键：--cm-map）
# 注意: 文中 /path/to/... 为占位符，运行前请按本机实际路径替换（见 config/）。
# ---------------------------------------------------------------------------

set -euo pipefail

BASE=/path/to/project/GWAS
GENO=$BASE/info_pass
POLYFUN=$BASE/PolyFun
MAP=$POLYFUN/genetic_map/genetic_map_b37

TMP=$POLYFUN/genotype_tmp
OUT=$POLYFUN/genotype_cm

mkdir -p "$TMP" "$OUT"

for chr in {1..22}; do

    echo
    echo "=========================================="
    echo "CHR${chr}"
    echo "=========================================="

    # --------------------------------------------------
    # 1. 创建临时 BIM，并给每一行唯一 ID
    # --------------------------------------------------
    awk '{
        print $1, "ROW" NR, $3, $4, $5, $6
    }' \
    "$GENO/${chr}.bim" \
    > "$TMP/${chr}.bim"

    ln -sf "$GENO/${chr}.bed" "$TMP/${chr}.bed"
    ln -sf "$GENO/${chr}.fam" "$TMP/${chr}.fam"

    # --------------------------------------------------
    # 2. 按 CHR + BP + 无序 allele pair 找重复
    #    每组保留第一条，后面的加入 exclude
    # --------------------------------------------------
    awk '
    {
        a=$5
        b=$6

        if (a > b) {
            tmp=a
            a=b
            b=tmp
        }

        key=$1 "\t" $4 "\t" a "\t" b

        if (seen[key]++) {
            print $2
        }
    }
    ' \
    "$TMP/${chr}.bim" \
    > "$TMP/${chr}.exclude"

    n_remove=$(wc -l < "$TMP/${chr}.exclude")

    echo "Duplicate variants to remove: $n_remove"

    # --------------------------------------------------
    # 3. 重新生成完整同步的 BED/BIM/FAM
    #    同时写入真实 genetic map CM
    # --------------------------------------------------
    plink \
        --bfile "$TMP/${chr}" \
        --exclude "$TMP/${chr}.exclude" \
        --cm-map \
        "$MAP/genetic_map_chr${chr}_combined_b37.txt" "$chr" \
        --make-bed \
        --out "$OUT/${chr}" \
        --allow-no-sex

    # --------------------------------------------------
    # 4. QC：variant 数量
    # --------------------------------------------------
    n_final=$(wc -l < "$OUT/${chr}.bim")
    n_original=$(wc -l < "$GENO/${chr}.bim")

    echo "Original variants : $n_original"
    echo "Removed variants  : $n_remove"
    echo "Final variants    : $n_final"

    if [ "$n_final" -ne "$((n_original - n_remove))" ]; then
        echo "ERROR: variant count mismatch on CHR${chr}"
        exit 1
    fi

    # --------------------------------------------------
    # 5. QC：检查最终 BIM 是否还有重复
    # --------------------------------------------------
    awk '
    {
        a=$5
        b=$6

        if (a > b) {
            tmp=a
            a=b
            b=tmp
        }

        key=$1 "\t" $4 "\t" a "\t" b
        count[key]++
    }

    END {
        n=0
        rec=0

        for (k in count) {
            if (count[k] > 1) {
                n++
                rec += count[k]
            }
        }

        printf "Final duplicate groups = %d\n", n
        printf "Final duplicate records = %d\n", rec

        if (n != 0) {
            exit 1
        }
    }
    ' "$OUT/${chr}.bim"

    # --------------------------------------------------
    # 6. QC：检查 CM
    # --------------------------------------------------
    awk '
    {
        cm=$3

        if (cm < min || NR==1)
            min=cm

        if (cm > max || NR==1)
            max=cm

        if (cm == 0)
            zero++

        if (cm < 0)
            negative++
    }

    END {
        printf "CM_min = %.8f\n", min
        printf "CM_max = %.8f\n", max
        printf "CM_zero = %d\n", zero
        printf "CM_negative = %d\n", negative
    }
    ' "$OUT/${chr}.bim"

done

echo
echo "=========================================="
echo "ALL CHROMOSOMES COMPLETED"
echo "=========================================="