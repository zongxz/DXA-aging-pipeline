# ---------------------------------------------------------------------------
# DXA-aging-pipeline / 04_downstream
# QC：精细定位输出的一致性检查（SNP/位置映射、PIP 累加）
# 注意: 文中 /path/to/... 为占位符，运行前请按本机实际路径替换（见 config/）。
# ---------------------------------------------------------------------------

python - <<'PY'
import pandas as pd
from pathlib import Path

for trait in ["double_right", "single_right_hip", "single_right_knee"]:
    candidates = [
        Path(f"{trait}_finemap_snp_summary.tsv"),
        Path(f"{trait}_finemap_snp_summary.tsv.gz"),
    ]
    f = next((p for p in candidates if p.exists()), None)

    if f is None:
        print(f"\n{trait}: summary file not found")
        continue

    d = pd.read_csv(f, sep="\t")
    print(f"\n========== {trait} ==========")
    print("File:", f.name)
    print("Total rows:", len(d))
    print("Unique SNP IDs:", d["SNP"].nunique())
    print("Unique CHR:BP:", d[["CHR", "BP"]].drop_duplicates().shape[0])

    # 同一 SNP ID 对应多个染色体位置
    snp_pos = d.groupby("SNP")[["CHR", "BP"]].apply(
        lambda x: len(x.drop_duplicates())
    )
    print("SNP IDs mapping to multiple positions:", int((snp_pos > 1).sum()))

    # 同一染色体位置对应多个 SNP ID
    pos_snp = d.groupby(["CHR", "BP"])["SNP"].nunique()
    print("Positions mapping to multiple SNP IDs:", int((pos_snp > 1).sum()))

    # 展示重复记录
    dup_snp = d[d.duplicated("SNP", keep=False)].sort_values(
        ["SNP", "CHR", "BP"]
    )
    print("\nDuplicate SNP examples:")
    print(
        dup_snp[["CHR", "SNP", "BP", "A1", "A2", "PIP_MAX", "N_REGIONS"]]
        .head(15).to_string(index=False)
    )

    dup_pos = d[d.duplicated(["CHR", "BP"], keep=False)].sort_values(
        ["CHR", "BP", "SNP"]
    )
    print("\nDuplicate position examples:")
    print(
        dup_pos[["CHR", "SNP", "BP", "A1", "A2", "PIP_MAX", "N_REGIONS"]]
        .head(15).to_string(index=False)
    )
PY


python - <<'PY'
import pandas as pd

for trait in ["double_right", "single_right_hip", "single_right_knee"]:
    f = f"{trait}_credible_sets.tsv"
    d = pd.read_csv(f, sep="\t")

    cs = (
        d.groupby(["REGION", "CREDIBLE_SET"])
         .agg(
             N_SNP=("SNP", "nunique"),
             PIP_SUM=("PIP", "sum"),
             PIP_MIN=("PIP", "min"),
             PIP_MAX=("PIP", "max"),
         )
         .reset_index()
    )

    print(f"\n========== {trait} ==========")
    print("Number of region-CS groups:", len(cs))
    print("\nPIP sum distribution:")
    print(cs["PIP_SUM"].describe().to_string())
    print("\nGroups with PIP sum < 0.90:", (cs["PIP_SUM"] < 0.90).sum())
    print("Groups with PIP sum >= 0.90:", (cs["PIP_SUM"] >= 0.90).sum())
    print("\nLowest PIP sums:")
    print(cs.nsmallest(10, "PIP_SUM").to_string(index=False))
PY