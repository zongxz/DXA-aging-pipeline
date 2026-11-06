#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# DXA-aging-pipeline / 04_downstream
# QC：诊断 SNPVAR 中的重复 rsID 并分类
# 注意: 文中 /path/to/... 为占位符，运行前请按本机实际路径替换（见 config/）。
# ---------------------------------------------------------------------------

import gzip
import os
from collections import defaultdict

SNPVAR = "finemap/single_right_knee_snpvar.tsv.gz"
GENO_DIR = "genotype_cm"
OUT = "finemap/duplicate_qc"

os.makedirs(OUT, exist_ok=True)

# ============================================================
# 1. 读取 genotype_cm
#
# exact:
#   CHR BP A1 A2
#
# swap:
#   CHR BP A2 A1
# ============================================================

geno_exact = set()
geno_swap = set()

print("Reading genotype_cm...")

for chrom in range(1, 23):

    bim = f"{GENO_DIR}/{chrom}.bim"

    with open(bim) as f:
        for line in f:

            x = line.rstrip("\n").split()

            if len(x) < 6:
                continue

            chr_, snp, cm, bp, a1, a2 = x[:6]

            geno_exact.add(
                (chr_, bp, a1, a2)
            )

            geno_swap.add(
                (chr_, bp, a2, a1)
            )

    print(f"  chr{chrom} done")

print()
print("genotype exact keys:", len(geno_exact))
print("genotype swap keys :", len(geno_swap))


# ============================================================
# 2. 读取 SNPVAR
#
# 只保存 duplicate SNP ID 的记录
# ============================================================

print()
print("Scanning SNPVAR...")

counts = defaultdict(int)

with gzip.open(SNPVAR, "rt") as f:

    header = f.readline().rstrip("\n").split("\t")
    idx = {x:i for i,x in enumerate(header)}

    for line in f:

        x = line.rstrip("\n").split("\t")

        snp = x[idx["SNP"]]

        counts[snp] += 1


dup_ids = {
    snp for snp, n in counts.items()
    if n > 1
}

print("Total SNP IDs :", len(counts))
print("Duplicate IDs :", len(dup_ids))


# ============================================================
# 3. 保存 duplicate records
# ============================================================

dup_records = defaultdict(list)

with gzip.open(SNPVAR, "rt") as f:

    header_line = f.readline()
    header = header_line.rstrip("\n").split("\t")
    idx = {x:i for i,x in enumerate(header)}

    for line in f:

        x = line.rstrip("\n").split("\t")

        snp = x[idx["SNP"]]

        if snp not in dup_ids:
            continue

        chr_ = x[idx["CHR"]]
        bp = x[idx["BP"]]
        a1 = x[idx["A1"]]
        a2 = x[idx["A2"]]

        dup_records[snp].append({
            "chr": chr_,
            "bp": bp,
            "a1": a1,
            "a2": a2,
            "line": line.rstrip("\n")
        })


# ============================================================
# 4. 对每个 duplicate SNP 分类
# ============================================================

classification = []

for snp, rows in dup_records.items():

    # 每一条 summary-stat record 的 genotype_cm 匹配状态
    statuses = []

    for r in rows:

        exact_key = (
            r["chr"],
            r["bp"],
            r["a1"],
            r["a2"]
        )

        swap_key = (
            r["chr"],
            r["bp"],
            r["a2"],
            r["a1"]
        )

        if exact_key in geno_exact:

            status = "EXACT"

        elif swap_key in geno_exact:

            status = "SWAP"

        else:

            status = "NO_MATCH"

        statuses.append(status)

    n_exact = statuses.count("EXACT")
    n_swap = statuses.count("SWAP")
    n_no = statuses.count("NO_MATCH")

    # --------------------------------------------------------
    # 分类
    # --------------------------------------------------------

    if n_exact == len(rows):

        category = "ALL_EXACT"

    elif n_exact + n_swap == len(rows):

        category = "ALL_MATCHED"

    elif n_exact + n_swap > 0:

        category = "PARTIAL_MATCH"

    else:

        category = "ALL_NO_MATCH"

    # --------------------------------------------------------
    # 特殊情况：
    #
    # 相同 SNP ID + 相同 allele pair 出现多次
    #
    # 例如 rs56158851
    # C CTGTTT
    # CTGTTT C
    #
    # 这两条实际上属于同一个 allele pair。
    # --------------------------------------------------------

    allele_groups = defaultdict(list)

    for r, status in zip(rows, statuses):

        key = (
            r["chr"],
            r["bp"],
            tuple(sorted([r["a1"], r["a2"]]))
        )

        allele_groups[key].append(
            (r, status)
        )

    same_allele_duplicate = any(
        len(v) > 1
        for v in allele_groups.values()
    )

    if same_allele_duplicate:

        category = "SAME_ALLELE_DUPLICATE"

    classification.append({
        "SNP": snp,
        "N": len(rows),
        "N_EXACT": n_exact,
        "N_SWAP": n_swap,
        "N_NO_MATCH": n_no,
        "CATEGORY": category,
        "rows": rows,
        "statuses": statuses
    })


# ============================================================
# 5. 输出 summary
# ============================================================

summary_file = f"{OUT}/duplicate_summary.tsv"

with open(summary_file, "w") as out:

    out.write(
        "SNP\tN\tN_EXACT\tN_SWAP\tN_NO_MATCH\tCATEGORY\n"
    )

    for z in sorted(
        classification,
        key=lambda x: x["SNP"]
    ):

        out.write(
            f"{z['SNP']}\t"
            f"{z['N']}\t"
            f"{z['N_EXACT']}\t"
            f"{z['N_SWAP']}\t"
            f"{z['N_NO_MATCH']}\t"
            f"{z['CATEGORY']}\n"
        )


# ============================================================
# 6. 输出详细结果
# ============================================================

detail_file = f"{OUT}/duplicate_detail.tsv"

with open(detail_file, "w") as out:

    out.write(
        "SNP\tCHR\tBP\tA1\tA2\tSTATUS\tCATEGORY\n"
    )

    for z in classification:

        for r, status in zip(
            z["rows"],
            z["statuses"]
        ):

            out.write(
                f"{z['SNP']}\t"
                f"{r['chr']}\t"
                f"{r['bp']}\t"
                f"{r['a1']}\t"
                f"{r['a2']}\t"
                f"{status}\t"
                f"{z['CATEGORY']}\n"
            )


# ============================================================
# 7. 统计
# ============================================================

cat_count = defaultdict(int)

for z in classification:
    cat_count[z["CATEGORY"]] += 1

print()
print("============================================================")
print("DUPLICATE QC SUMMARY")
print("============================================================")

for k in [
    "ALL_EXACT",
    "ALL_MATCHED",
    "PARTIAL_MATCH",
    "SAME_ALLELE_DUPLICATE",
    "ALL_NO_MATCH"
]:

    print(
        f"{k:25s}: {cat_count[k]}"
    )

print()
print("Summary:")
print(summary_file)

print()
print("Detail:")
print(detail_file)


# ============================================================
# 8. 特别输出 SAME_ALLELE_DUPLICATE
# ============================================================

special = [
    z for z in classification
    if z["CATEGORY"] == "SAME_ALLELE_DUPLICATE"
]

print()
print("============================================================")
print("SAME ALLELE DUPLICATES")
print("============================================================")

for z in special[:100]:

    print()
    print(
        z["SNP"],
        "N =", z["N"]
    )

    for r, status in zip(
        z["rows"],
        z["statuses"]
    ):

        print(
            "   ",
            r["chr"],
            r["bp"],
            r["a1"],
            r["a2"],
            status
        )

print()
print(
    "Showing first",
    min(100, len(special)),
    "of",
    len(special)
)
