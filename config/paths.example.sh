# ---------------------------------------------------------------------------
# DXA-aging-pipeline / config
# 路径配置模板（bash 侧）。复制为 config/paths.sh 后按需修改。
#   cp config/paths.example.sh config/paths.sh
#   source config/paths.sh
# 注意: 本仓库代码里写的是 /path/to/... 占位符，可用
#       bash scripts/instantiate_paths.sh 一次性替换。
# ---------------------------------------------------------------------------

# ---- 项目根 ----
export PROJECT_ROOT="/path/to/project"      # Python 侧：元数据 / GWAS / 下游结果
export EPI_ROOT="/path/to/epi"              # R 侧：F1–F6 的 data/ 与 result/

# ---- 共享参考数据 ----
export REF_ROOT="/path/to/ref"              # 1000G EUR、HapMap3、PRS-CS LD 块、LDSC 权重、GTEx
export RAW_GENO_DIR="/path/to/raw_genotype" # QC 前原始插补基因型（per-chromosome bed/bim/fam）

# ---- 上游共享表 ----
export SHARED_ROOT="/path/to/shared_project" # ID 对照表、蛋白组矩阵、scPagwas 辅助数据

# ---- 第三方工具 ----
export TOOLS_ROOT="/path/to/tools"
export PLINK2="${TOOLS_ROOT}/plink2"
export PLINK19="${TOOLS_ROOT}/plink"
export GCTA="${TOOLS_ROOT}/gcta64"
export LDSC="${TOOLS_ROOT}/ldsc/ldsc.py"
export PRSCS="${TOOLS_ROOT}/PRScs/PRScs.py"
export POLYFUN="${TOOLS_ROOT}/polyfun"
export SPREIDXCAN="${TOOLS_ROOT}/MetaXcan/software/SPrediXcan.py"

# ---- MiXeR 容器 ----
export MIXER_IMAGE="ghcr.io/precimed/gsa-mixer:latest"
export MIXER_DIR="${PROJECT_ROOT}/GWAS/MiXeR"
export MIXER_THREADS=24
export MIXER_N_REP=20

# ---- 通用 ----
export N_CORES=20          # LAVA / doParallel 并发数
export R_BIN="Rscript"
export PY_BIN="python"
