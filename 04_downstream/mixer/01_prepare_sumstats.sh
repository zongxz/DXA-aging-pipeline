# ---------------------------------------------------------------------------
# DXA-aging-pipeline / 04_downstream
# 整理 GWAS 为 MiXeR 输入（SNP/CHR/BP/A1/A2/N/Z）
# 注意: 文中 /path/to/... 为占位符，运行前请按本机实际路径替换（见 config/）。
# ---------------------------------------------------------------------------
REF=/path/to/project/GWAS/MiXeR/ref
GWAS=/path/to/project/GWAS/gwas_result
OUT=/path/to/project/GWAS/MiXeR/sumstats

mkdir -p "$OUT"

python - <<'PY'
import pandas as pd
import numpy as np
import gzip
import os

GWAS = os.path.expanduser(
    "/path/to/project/GWAS/gwas_result"
)
OUT = os.path.expanduser(
    "/path/to/project/GWAS/MiXeR/sumstats"
)

os.makedirs(OUT, exist_ok=True)

traits = {
   # "single_right_knee": "single_right_knee.txt",
   # "single_right_hip":  "single_right_hip.txt",
   "double_right":  "double_right.txt",
}

for trait, fn in traits.items():

    infile = os.path.join(GWAS, fn)
    outfile = os.path.join(OUT, trait + ".sumstats.gz")

    print(f"\n===== {trait} =====")
    print("Reading:", infile)

    df = pd.read_csv(infile, sep=r"\s+", dtype={
        "SNP": str,
        "CHR": str,
        "A1": str,
        "A2": str
    })

    required = ["CHR", "SNP", "POS", "A1", "A2", "N", "BETA", "SE"]

    missing = [x for x in required if x not in df.columns]
    if missing:
        raise ValueError(f"{trait}: missing columns: {missing}")

    # MiXeR Z score
    df["Z"] = df["BETA"] / df["SE"]

    # 保留 MiXeR 需要的字段
    out = df[
        ["SNP", "CHR", "POS", "A1", "A2", "N", "Z"]
    ].copy()

    out = out.rename(columns={"POS": "BP"})

    # 数值化
    out["N"] = pd.to_numeric(out["N"], errors="coerce")
    out["Z"] = pd.to_numeric(out["Z"], errors="coerce")

    # allele 大写
    out["A1"] = out["A1"].str.upper()
    out["A2"] = out["A2"].str.upper()

    # QC
    n0 = len(out)

    out = out.dropna(
        subset=["SNP", "CHR", "BP", "A1", "A2", "N", "Z"]
    )

    # 去除非法 allele
    out = out[
        out["A1"].isin(["A", "C", "G", "T"]) &
        out["A2"].isin(["A", "C", "G", "T"])
    ]

    # 去除 A1 == A2
    out = out[out["A1"] != out["A2"]]

    # 去除 N <= 0
    out = out[out["N"] > 0]

    # 去除 Z 非有限值
    out = out[np.isfinite(out["Z"])]

    # 去除重复 SNP
    out = out.drop_duplicates("SNP", keep="first")

    # 排序
    out["CHR_NUM"] = pd.to_numeric(out["CHR"], errors="coerce")
    out["BP"] = pd.to_numeric(out["BP"], errors="coerce")

    out = out.dropna(subset=["CHR_NUM", "BP"])

    out = out[out["CHR_NUM"].between(1, 22)]

    out = out.sort_values(["CHR_NUM", "BP"])

    out = out.drop(columns=["CHR_NUM"])

    # 输出
    out.to_csv(
        outfile,
        sep="\t",
        index=False,
        compression="gzip"
    )

    print("Original SNPs :", n0)
    print("Final SNPs    :", len(out))
    print("Median N      :", out["N"].median())
    print("Min N         :", out["N"].min())
    print("Max N         :", out["N"].max())
    print("Median |Z|    :", out["Z"].abs().median())
    print("Output        :", outfile)
PY