# ---------------------------------------------------------------------------
# DXA-aging-pipeline / config
# 路径配置模板（R 侧）。复制为 config/paths.R 后按需修改。
#   cp config/paths.example.R config/paths.R
# 用法: 在 R 脚本开头 source("config/paths.R")，之后用 PROJECT_ROOT / REF_ROOT 等变量。
# 注意: 本仓库代码里写的是 /path/to/... 占位符，可用
#       bash scripts/instantiate_paths.sh 一次性替换。
# ---------------------------------------------------------------------------

PROJECT_ROOT  <- Sys.getenv("PROJECT_ROOT",  "/path/to/project")
EPI_ROOT      <- Sys.getenv("EPI_ROOT",      "/path/to/epi")
REF_ROOT      <- Sys.getenv("REF_ROOT",      "/path/to/ref")
SHARED_ROOT   <- Sys.getenv("SHARED_ROOT",   "/path/to/shared_project")
TOOLS_ROOT    <- Sys.getenv("TOOLS_ROOT",    "/path/to/tools")
RAW_GENO_DIR  <- Sys.getenv("RAW_GENO_DIR",  "/path/to/raw_genotype")

# ---- 常用子目录 ----
GWAS_DIR      <- file.path(PROJECT_ROOT, "GWAS")
GWAS_RESULT   <- file.path(GWAS_DIR, "gwas_result")
LDSC_DIR      <- file.path(GWAS_DIR, "LDSC")
LAVA_DIR      <- file.path(GWAS_DIR, "LAVA")
NGWAMA_DIR    <- file.path(GWAS_DIR, "NGWAMA")
PRSCS_DIR     <- file.path(GWAS_DIR, "PRS_CS")
POLYFUN_DIR   <- file.path(GWAS_DIR, "PolyFun")
EPI_DATA      <- file.path(EPI_ROOT, "data")
EPI_RESULT    <- file.path(EPI_ROOT, "result")

# ---- 参考数据 ----
G1000_EUR     <- file.path(REF_ROOT, "g1000_eur", "g1000_eur")
HAPMAP3_SNPS  <- file.path(REF_ROOT, "hapmap", "w_hm3.snplist")
LDSC_EUR_LD   <- file.path(REF_ROOT, "eur_w_ld_chr")

# ---- 并发 ----
N_CORES <- as.integer(Sys.getenv("N_CORES", 20))
