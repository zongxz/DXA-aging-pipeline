#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# DXA-aging-pipeline / 04_downstream
# 第 8 步：核心——修重复位点、切区域汇总统计、并发跑 SuSiE finemapper
# 注意: 文中 /path/to/... 为占位符，运行前请按本机实际路径替换（见 config/）。
# ---------------------------------------------------------------------------
set -u
set -o pipefail

###############################################################################
# PolyFun local fine-mapping
# Trait: double_right
#
# Approach 3
#
# Strategy
# ---------------------------------------------------------------------------
# 1. Collect all local fine-mapping jobs.
#
# 2. Read the full SNPVAR file once.
#
# 3. Read duplicate QC:
#
#       finemap/duplicate_qc/duplicate_detail.tsv
#
#    Format:
#
#       SNP CHR BP A1 A2 STATUS CATEGORY
#
# 4. Only repair:
#
#       CATEGORY = SAME_ALLELE_DUPLICATE
#
#    For these SNPs:
#
#       STATUS = EXACT  -> KEEP
#       STATUS = SWAP   -> REMOVE
#
#    ALL_EXACT and ALL_MATCHED are NOT modified.
#
# 5. Generate region-specific sumstats.
#
#    Existing unaffected regions are kept.
#    Existing affected regions are backed up and regenerated.
#
# 6. Run local PLINK LD:
#
#       --geno genotype_cm/<chr>
#
# 7. Use:
#
#       --allow-swapped-indel-alleles
#
#    but NEVER:
#
#       --allow-missing
#
# 8. Existing successful fine-mapping results are skipped.
#
# 9. Failed jobs are recorded in:
#
#       finemap/failed_jobs.txt
###############################################################################

BASE="/path/to/project/GWAS/PolyFun"

POLYFUN="$BASE/software/polyfun"
FINEMAPPER="$POLYFUN/finemapper.py"

SUMSTATS="$BASE/finemap/double_right_snpvar.tsv.gz"

JOBS_DIR="$BASE/finemap/jobs_local"

REGION_SUMSTATS_DIR="$BASE/finemap/region_sumstats"

LOG_DIR="$BASE/finemap/logs"

RESULT_DIR="$BASE/finemap"

DUP_QC_DIR="$BASE/finemap/duplicate_qc"

ALL_JOBS="$BASE/finemap/all_local_jobs.txt"

RUN_JOBS="$BASE/finemap/run_local_region_jobs.txt"

FAILED_JOBS="$BASE/finemap/failed_jobs.txt"

MAX_JOBS=8

N=48850
MEMORY=1
MAX_NUM_CAUSAL=2

cd "$BASE" || exit 1

mkdir -p "$REGION_SUMSTATS_DIR"
mkdir -p "$LOG_DIR"
mkdir -p "$DUP_QC_DIR"

###############################################################################
# 1. Check input
###############################################################################

echo
echo "============================================================"
echo "[INFO] Checking input files"
echo "============================================================"
echo

if [[ ! -f "$SUMSTATS" ]]; then
    echo "[ERROR] SNPVAR not found:"
    echo "        $SUMSTATS"
    exit 1
fi

if [[ ! -f "$FINEMAPPER" ]]; then
    echo "[ERROR] finemapper.py not found:"
    echo "        $FINEMAPPER"
    exit 1
fi

if [[ ! -f "$DUP_QC_DIR/duplicate_detail.tsv" ]]; then
    echo "[ERROR] duplicate_detail.tsv not found:"
    echo "        $DUP_QC_DIR/duplicate_detail.tsv"
    echo
    echo "Run first:"
    echo "  python scripts/check_duplicate_snps.py"
    exit 1
fi

###############################################################################
# 2. Check genotype_cm
###############################################################################

echo "[INFO] Checking genotype_cm..."

for chr in {1..22}; do

    for ext in bed bim fam; do

        FILE="$BASE/genotype_cm/${chr}.${ext}"

        if [[ ! -f "$FILE" ]]; then
            echo "[ERROR] Missing:"
            echo "        $FILE"
            exit 1
        fi

    done

done

echo "[INFO] genotype_cm check passed."

###############################################################################
# 3. Collect all jobs
###############################################################################

echo
echo "============================================================"
echo "[INFO] Collecting local jobs"
echo "============================================================"
echo

: > "$ALL_JOBS"

for chr in {1..22}; do

    JOB_FILE="$JOBS_DIR/double_right_chr${chr}.jobs"

    if [[ ! -f "$JOB_FILE" ]]; then
        echo "[ERROR] Missing job file:"
        echo "        $JOB_FILE"
        exit 1
    fi

    cat "$JOB_FILE" >> "$ALL_JOBS"

done

TOTAL_JOBS=$(wc -l < "$ALL_JOBS")

echo "[INFO] Total jobs: $TOTAL_JOBS"

if [[ "$TOTAL_JOBS" -ne 2748 ]]; then
    echo "[WARNING] Expected 2748 jobs, found $TOTAL_JOBS"
fi

###############################################################################
# 4. Validate original jobs
###############################################################################

echo
echo "[INFO] Validating original jobs..."

if grep -q -- '--ld ' "$ALL_JOBS"; then

    echo "[ERROR] --ld found in original jobs:"
    grep -- '--ld ' "$ALL_JOBS" | head

    exit 1

fi

if ! grep -q -- '--geno genotype_cm/' "$ALL_JOBS"; then

    echo "[ERROR] No local --geno jobs found."

    exit 1

fi

echo "[INFO] Original job validation passed."

###############################################################################
# 5. Generate / repair region-specific sumstats
###############################################################################

echo
echo "============================================================"
echo "[INFO] Preparing region-specific sumstats"
echo "============================================================"
echo

python - <<'PY'

import os
import shlex
import pandas as pd

BASE = "/path/to/project/GWAS/PolyFun"

SUMSTATS = os.path.join(
    BASE,
    "finemap",
    "double_right_snpvar.tsv.gz"
)

JOBS = os.path.join(
    BASE,
    "finemap",
    "all_local_jobs.txt"
)

OUTDIR = os.path.join(
    BASE,
    "finemap",
    "region_sumstats"
)

DETAIL = os.path.join(
    BASE,
    "finemap",
    "duplicate_qc",
    "duplicate_detail.tsv"
)

os.makedirs(OUTDIR, exist_ok=True)

###############################################################################
# Parse jobs
###############################################################################

regions = []

with open(JOBS) as f:

    for line in f:

        line = line.strip()

        if not line:
            continue

        args = shlex.split(line)

        def get_arg(name):

            try:
                i = args.index(name)
                return args[i + 1]

            except (ValueError, IndexError):

                raise RuntimeError(
                    f"Cannot parse {name} from job:\n{line}"
                )

        chr_num = int(get_arg("--chr"))
        start = int(get_arg("--start"))
        end = int(get_arg("--end"))
        out = get_arg("--out")

        regions.append(
            {
                "chr": chr_num,
                "start": start,
                "end": end,
                "out": out
            }
        )

print(
    f"[INFO] Parsed regions: {len(regions)}"
)

###############################################################################
# Read duplicate QC
###############################################################################

dup = pd.read_csv(
    DETAIL,
    sep="\t",
    dtype=str
)

required_dup = [
    "SNP",
    "CHR",
    "BP",
    "A1",
    "A2",
    "STATUS",
    "CATEGORY"
]

missing_dup = [
    x for x in required_dup
    if x not in dup.columns
]

if missing_dup:

    raise RuntimeError(
        "Missing columns in duplicate_detail.tsv: "
        + str(missing_dup)
    )

###############################################################################
# Identify SAME_ALLELE_DUPLICATE
###############################################################################

same = dup[
    dup["CATEGORY"] == "SAME_ALLELE_DUPLICATE"
].copy()

same_ids = set(
    same["SNP"].astype(str)
)

print(
    f"[INFO] SAME_ALLELE_DUPLICATE SNPs: "
    f"{len(same_ids)}"
)

if len(same_ids) != 66:

    raise RuntimeError(
        f"Expected 66 SAME_ALLELE_DUPLICATE SNPs, "
        f"found {len(same_ids)}"
    )

###############################################################################
# Validate status
###############################################################################

status_count = (
    same
    .groupby(
        ["SNP", "STATUS"]
    )
    .size()
    .unstack(
        fill_value=0
    )
)

if "EXACT" not in status_count.columns:
    raise RuntimeError(
        "No EXACT records found."
    )

if "SWAP" not in status_count.columns:
    raise RuntimeError(
        "No SWAP records found."
    )

bad = status_count[
    (status_count["EXACT"] != 1) |
    (status_count["SWAP"] != 1)
]

if len(bad) > 0:

    print(
        "[ERROR] Unexpected SAME_ALLELE_DUPLICATE structure:"
    )

    print(
        bad.to_string()
    )

    raise RuntimeError(
        "SAME_ALLELE_DUPLICATE QC failed."
    )

print(
    "[INFO] All 66 SNPs have exactly "
    "1 EXACT + 1 SWAP."
)

###############################################################################
# Build exact / swap keys
###############################################################################

same_exact = same[
    same["STATUS"] == "EXACT"
].copy()

same_swap = same[
    same["STATUS"] == "SWAP"
].copy()

exact_keys = set(
    zip(
        same_exact["SNP"],
        same_exact["CHR"].astype(int),
        same_exact["BP"].astype(int),
        same_exact["A1"],
        same_exact["A2"]
    )
)

swap_keys = set(
    zip(
        same_swap["SNP"],
        same_swap["CHR"].astype(int),
        same_swap["BP"].astype(int),
        same_swap["A1"],
        same_swap["A2"]
    )
)

###############################################################################
# Read full SNPVAR
###############################################################################

print()
print("[INFO] Reading full SNPVAR...")

required = [
    "CHR",
    "SNP",
    "BP",
    "A1",
    "A2",
    "SNPVAR",
    "N",
    "Z"
]

df = pd.read_csv(
    SUMSTATS,
    sep="\t"
)

missing = [
    x for x in required
    if x not in df.columns
]

if missing:

    raise RuntimeError(
        "Missing SNPVAR columns: "
        + str(missing)
    )

df["CHR"] = df["CHR"].astype(int)
df["BP"] = df["BP"].astype(int)
df["SNP"] = df["SNP"].astype(str)
df["A1"] = df["A1"].astype(str)
df["A2"] = df["A2"].astype(str)

print(
    f"[INFO] Full SNPVAR rows: "
    f"{len(df):,}"
)

###############################################################################
# Check the 66 duplicate SNPs
###############################################################################

same_mask = df["SNP"].isin(
    same_ids
)

same_df = df.loc[
    same_mask,
    [
        "SNP",
        "CHR",
        "BP",
        "A1",
        "A2"
    ]
].copy()

print(
    f"[INFO] Rows belonging to the 66 SNPs: "
    f"{len(same_df)}"
)

if len(same_df) != 132:

    raise RuntimeError(
        f"Expected 132 rows "
        f"(66 SNPs x 2), found {len(same_df)}"
    )

###############################################################################
# Remove only SWAP rows
###############################################################################

remove_index = []

for idx, row in same_df.iterrows():

    key = (
        row["SNP"],
        int(row["CHR"]),
        int(row["BP"]),
        row["A1"],
        row["A2"]
    )

    if key in swap_keys:

        remove_index.append(idx)

    elif key in exact_keys:

        pass

    else:

        raise RuntimeError(
            "Unexpected SAME_ALLELE record:\n"
            + str(key)
        )

print(
    f"[INFO] SWAP rows to remove: "
    f"{len(remove_index)}"
)

if len(remove_index) != 66:

    raise RuntimeError(
        f"Expected 66 SWAP rows, "
        f"found {len(remove_index)}"
    )

###############################################################################
# Create repaired SNPVAR in memory
###############################################################################

df_fixed = df.drop(
    index=remove_index
).copy()

print(
    f"[INFO] SNPVAR rows after repair: "
    f"{len(df_fixed):,}"
)

###############################################################################
# Final check: each SAME_ALLELE SNP exactly once
###############################################################################

remaining_count = (
    df_fixed[
        df_fixed["SNP"].isin(same_ids)
    ]
    .groupby("SNP")
    .size()
)

bad_remaining = remaining_count[
    remaining_count != 1
]

if len(bad_remaining) > 0:

    print(
        "[ERROR] Remaining duplicate counts:"
    )

    print(
        bad_remaining.to_string()
    )

    raise RuntimeError(
        "Post-repair duplicate QC failed."
    )

print(
    "[INFO] All 66 SAME_ALLELE_DUPLICATE SNPs "
    "now occur exactly once."
)

###############################################################################
# Generate region-specific sumstats
###############################################################################

n_existing = 0
n_created = 0
n_repaired = 0

affected_regions = []

for i, r in enumerate(
    regions,
    start=1
):

    chr_num = r["chr"]
    start = r["start"]
    end = r["end"]

    filename = (
        f"double_right."
        f"chr{chr_num}."
        f"{start}_{end}.tsv.gz"
    )

    outfile = os.path.join(
        OUTDIR,
        filename
    )

    ###########################################################################
    # Original region
    ###########################################################################

    region_original = df[
        (df["CHR"] == chr_num) &
        (df["BP"] >= start) &
        (df["BP"] <= end)
    ][required].copy()

    ###########################################################################
    # Check whether this region contains one of the 66 problematic SNPs
    ###########################################################################

    affected = region_original[
        "SNP"
    ].isin(
        same_ids
    ).any()

    ###########################################################################
    # Unaffected existing region:
    # keep exactly as it is.
    ###########################################################################

    if (
        os.path.exists(outfile)
        and os.path.getsize(outfile) > 0
        and not affected
    ):

        n_existing += 1
        continue

    ###########################################################################
    # Generate repaired region
    ###########################################################################

    region_fixed = df_fixed[
        (df_fixed["CHR"] == chr_num) &
        (df_fixed["BP"] >= start) &
        (df_fixed["BP"] <= end)
    ][required].copy()

    if len(region_fixed) == 0:

        print(
            f"[WARNING] Empty region: "
            f"chr{chr_num}:{start}-{end}"
        )

        continue

    ###########################################################################
    # Backup existing affected region
    ###########################################################################

    if (
        affected
        and os.path.exists(outfile)
        and os.path.getsize(outfile) > 0
    ):

        backup = (
            outfile +
            ".before_duplicate_fix"
        )

        if not os.path.exists(backup):

            os.rename(
                outfile,
                backup
            )

            print(
                f"[BACKUP] {os.path.basename(outfile)}"
            )

    ###########################################################################
    # Write region sumstats
    ###########################################################################

    region_fixed.to_csv(
        outfile,
        sep="\t",
        index=False,
        compression="gzip"
    )

    if affected:

        n_repaired += 1

        affected_regions.append(
            (
                chr_num,
                start,
                end
            )
        )

    else:

        n_created += 1

    if (
        i == 1
        or i % 100 == 0
        or i == len(regions)
    ):

        print(
            f"[INFO] {i}/{len(regions)} regions processed"
        )

###############################################################################
# Save affected region list
###############################################################################

affected_file = os.path.join(
    BASE,
    "finemap",
    "duplicate_qc",
    "affected_regions.tsv"
)

with open(
    affected_file,
    "w"
) as f:

    f.write(
        "CHR\tSTART\tEND\n"
    )

    for x in affected_regions:

        f.write(
            f"{x[0]}\t{x[1]}\t{x[2]}\n"
        )

###############################################################################
# Summary
###############################################################################

print()
print(
    "============================================================"
)

print(
    "[INFO] Region sumstats summary"
)

print(
    "============================================================"
)

print(
    f"Total regions:             {len(regions)}"
)

print(
    f"Existing unaffected:       {n_existing}"
)

print(
    f"Newly created:             {n_created}"
)

print(
    f"Repaired affected:         {n_repaired}"
)

print(
    f"Affected region list:      {affected_file}"
)

print(
    "============================================================"
)

###############################################################################
# Final QC of affected regions
###############################################################################

print()
print(
    "[INFO] Running final region duplicate QC..."
)

bad_regions = []

for chr_num, start, end in affected_regions:

    filename = (
        f"double_right."
        f"chr{chr_num}."
        f"{start}_{end}.tsv.gz"
    )

    outfile = os.path.join(
        OUTDIR,
        filename
    )

    x = pd.read_csv(
        outfile,
        sep="\t",
        usecols=["SNP"]
    )

    counts = x["SNP"].value_counts()

    duplicated = counts[
        counts > 1
    ]

    duplicated_same = [
        snp
        for snp in duplicated.index
        if snp in same_ids
    ]

    if duplicated_same:

        bad_regions.append(
            (
                chr_num,
                start,
                end,
                duplicated_same
            )
        )

if bad_regions:

    print(
        "[ERROR] SAME_ALLELE_DUPLICATE "
        "still duplicated:"
    )

    for x in bad_regions:
        print(x)

    raise RuntimeError(
        "Final region duplicate QC failed."
    )

print(
    "[INFO] Final region duplicate QC passed."
)

print(
    "[INFO] Region sumstats preparation completed."
)

PY

###############################################################################
# 6. Validate region sumstats
###############################################################################

echo
echo "============================================================"
echo "[INFO] Validating region-specific sumstats"
echo "============================================================"
echo

N_REGION_SUMSTATS=$(find \
    "$REGION_SUMSTATS_DIR" \
    -type f \
    -name '*.tsv.gz' \
    | wc -l)

echo "[INFO] Region sumstats files: $N_REGION_SUMSTATS"

if [[ "$N_REGION_SUMSTATS" -eq 0 ]]; then

    echo "[ERROR] No region-specific sumstats found."

    exit 1

fi

###############################################################################
# 7. Create execution commands
###############################################################################

echo
echo "============================================================"
echo "[INFO] Creating execution commands"
echo "============================================================"
echo

python - <<'PY'

import os
import shlex

BASE = "/path/to/project/GWAS/PolyFun"

JOBS = os.path.join(
    BASE,
    "finemap",
    "all_local_jobs.txt"
)

OUTFILE = os.path.join(
    BASE,
    "finemap",
    "run_local_region_jobs.txt"
)

REGION_DIR = os.path.join(
    BASE,
    "finemap",
    "region_sumstats"
)

commands = []

with open(JOBS) as f:

    for line in f:

        line = line.strip()

        if not line:
            continue

        args = shlex.split(line)

        def get_arg(name):

            try:

                i = args.index(name)

                return args[i + 1]

            except (
                ValueError,
                IndexError
            ):

                raise RuntimeError(
                    f"Cannot parse {name}:\n{line}"
                )

        chr_num = int(
            get_arg("--chr")
        )

        start = int(
            get_arg("--start")
        )

        end = int(
            get_arg("--end")
        )

        out = get_arg("--out")

        region_sumstats = os.path.join(
            REGION_DIR,
            f"double_right."
            f"chr{chr_num}."
            f"{start}_{end}.tsv.gz"
        )

        if not os.path.exists(
            region_sumstats
        ):

            raise RuntimeError(
                "Missing region sumstats:\n"
                + region_sumstats
            )

        cmd = [
            "python",
            "software/polyfun/finemapper.py",

            "--chr",
            str(chr_num),

            "--start",
            str(start),

            "--end",
            str(end),

            "--out",
            out,

            "--method",
            "susie",

            "--sumstats",
            region_sumstats,

            "--n",
            "48850",

            "--memory",
            "1",

            "--max-num-causal",
            "2",

            "--geno",
            f"genotype_cm/{chr_num}",

            "--allow-swapped-indel-alleles"
        ]

        commands.append(
            " ".join(
                shlex.quote(x)
                for x in cmd
            )
        )

with open(
    OUTFILE,
    "w"
) as f:

    for cmd in commands:

        f.write(
            cmd + "\n"
        )

print(
    f"[INFO] Generated "
    f"{len(commands)} execution commands"
)

print(
    f"[INFO] Output: {OUTFILE}"
)

PY

TOTAL_RUN_JOBS=$(wc -l < "$RUN_JOBS")

echo "[INFO] Execution jobs: $TOTAL_RUN_JOBS"

if [[ "$TOTAL_RUN_JOBS" -ne "$TOTAL_JOBS" ]]; then

    echo "[ERROR] Number of execution jobs does not match."

    exit 1

fi

###############################################################################
# 8. Validate execution commands
###############################################################################

echo
echo "============================================================"
echo "[INFO] Validating execution commands"
echo "============================================================"
echo

if grep -q -- '--ld ' "$RUN_JOBS"; then

    echo "[ERROR] --ld found."

    exit 1

fi

if grep -q -- '--allow-missing' "$RUN_JOBS"; then

    echo "[ERROR] --allow-missing found."

    exit 1

fi

N_LOCAL=$(grep -c \
    -- '--geno genotype_cm/' \
    "$RUN_JOBS")

N_SWAPPED=$(grep -c \
    -- '--allow-swapped-indel-alleles' \
    "$RUN_JOBS")

echo "[INFO] Local LD jobs:"
echo "       $N_LOCAL / $TOTAL_RUN_JOBS"

echo "[INFO] allow-swapped-indel-alleles:"
echo "       $N_SWAPPED / $TOTAL_RUN_JOBS"

if [[ "$N_LOCAL" -ne "$TOTAL_RUN_JOBS" ]]; then

    echo "[ERROR] Not all jobs use local genotype_cm."

    exit 1

fi

if [[ "$N_SWAPPED" -ne "$TOTAL_RUN_JOBS" ]]; then

    echo "[ERROR] Not all jobs contain "
    echo "        --allow-swapped-indel-alleles"

    exit 1

fi

###############################################################################
# 9. Run fine-mapping
###############################################################################

echo
echo "============================================================"
echo "[INFO] Starting fine-mapping"
echo "============================================================"

echo "[INFO] MAX_JOBS = $MAX_JOBS"
echo

RUNNING=0
STARTED=0
COMPLETED=0
SKIPPED=0

PIDS=()
PID_JOB_IDS=()

: > "$FAILED_JOBS"

while IFS= read -r CMD; do

    [[ -z "$CMD" ]] && continue

    STARTED=$((STARTED + 1))

    ###########################################################################
    # Extract output file
    ###########################################################################

    OUT=$(
        echo "$CMD" |
        sed -n \
        's/.*--out \([^ ]*\).*/\1/p'
    )

    ###########################################################################
    # Skip existing successful/nonempty result
    ###########################################################################

    if [[
        -n "$OUT" &&
        -f "$OUT" &&
        -s "$OUT"
    ]]; then

        SKIPPED=$((SKIPPED + 1))

        echo "[SKIP] $OUT"

        continue

    fi

    ###########################################################################
    # Start job
    ###########################################################################

    JOB_ID="$STARTED"

    LOG="$LOG_DIR/job_${JOB_ID}.log"

    echo
    echo "[START] Job $JOB_ID"
    echo "        Log: $LOG"
    echo "        $CMD"

    (
        echo "============================================================"
        echo "JOB ID: $JOB_ID"
        echo "START: $(date)"
        echo "============================================================"
        echo
        echo "$CMD"
        echo

        eval "$CMD"

        STATUS=$?

        echo
        echo "============================================================"
        echo "END: $(date)"
        echo "STATUS: $STATUS"
        echo "============================================================"

        exit $STATUS

    ) > "$LOG" 2>&1 &

    PID=$!

    PIDS+=(
        "$PID"
    )

    PID_JOB_IDS+=(
        "$JOB_ID"
    )

    RUNNING=$((RUNNING + 1))

    ###########################################################################
    # Maintain MAX_JOBS
    ###########################################################################

    if [[ "$RUNNING" -ge "$MAX_JOBS" ]]; then

        wait "${PIDS[0]}"

        STATUS=$?

        FINISHED_JOB="${PID_JOB_IDS[0]}"

        if [[ "$STATUS" -eq 0 ]]; then

            COMPLETED=$((COMPLETED + 1))

            echo "[DONE] Job $FINISHED_JOB"

        else

            echo "$FINISHED_JOB" \
                >> "$FAILED_JOBS"

            echo "[FAIL] Job $FINISHED_JOB"

        fi

        PIDS=(
            "${PIDS[@]:1}"
        )

        PID_JOB_IDS=(
            "${PID_JOB_IDS[@]:1}"
        )

        RUNNING=$((RUNNING - 1))

    fi

done < "$RUN_JOBS"

###############################################################################
# 10. Wait for remaining jobs
###############################################################################

echo
echo "[INFO] Waiting for remaining jobs..."

while [[ "$RUNNING" -gt 0 ]]; do

    wait "${PIDS[0]}"

    STATUS=$?

    FINISHED_JOB="${PID_JOB_IDS[0]}"

    if [[ "$STATUS" -eq 0 ]]; then

        COMPLETED=$((COMPLETED + 1))

        echo "[DONE] Job $FINISHED_JOB"

    else

        echo "$FINISHED_JOB" \
            >> "$FAILED_JOBS"

        echo "[FAIL] Job $FINISHED_JOB"

    fi

    PIDS=(
        "${PIDS[@]:1}"
    )

    PID_JOB_IDS=(
        "${PID_JOB_IDS[@]:1}"
    )

    RUNNING=$((RUNNING - 1))

done

###############################################################################
# 11. Final summary
###############################################################################

RESULT_COUNT=$(find \
    "$RESULT_DIR" \
    -maxdepth 1 \
    -type f \
    -name 'double_right.chr*.gz' \
    | wc -l)

FAILED_COUNT=0

if [[ -s "$FAILED_JOBS" ]]; then

    FAILED_COUNT=$(
        wc -l < "$FAILED_JOBS"
    )

fi

echo
echo "============================================================"
echo "Fine-mapping summary"
echo "============================================================"

echo "Total jobs:       $TOTAL_RUN_JOBS"
echo "Started:          $STARTED"
echo "Completed:        $COMPLETED"
echo "Skipped:          $SKIPPED"
echo "Failed:           $FAILED_COUNT"
echo "Result files:     $RESULT_COUNT"
echo "MAX_JOBS:         $MAX_JOBS"

echo
echo "Region sumstats:"
echo "  $REGION_SUMSTATS_DIR"

echo
echo "Affected regions:"
echo "  $DUP_QC_DIR/affected_regions.tsv"

echo
echo "Logs:"
echo "  $LOG_DIR"

echo
echo "Failed jobs:"
echo "  $FAILED_JOBS"

echo "============================================================"

###############################################################################
# 12. Exit status
###############################################################################

if [[ "$FAILED_COUNT" -gt 0 ]]; then

    exit 1

fi

exit 0