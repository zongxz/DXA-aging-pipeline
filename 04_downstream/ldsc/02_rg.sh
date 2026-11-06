# ---------------------------------------------------------------------------
# DXA-aging-pipeline / 04_downstream
# LDSC 第 2 步：成对遗传相关 rg 与跨性状截距
# 注意: 文中 /path/to/... 为占位符，运行前请按本机实际路径替换（见 config/）。
# ---------------------------------------------------------------------------
/path/to/shared_project/gwas/LDreg/ldsc/ldsc.py \
    --rg \
    /path/to/project/GWAS/LDSC/double_right.sumstats.gz,/path/to/project/GWAS/LDSC/single_right_knee.sumstats.gz \
    --ref-ld-chr /path/to/shared_project/gwas/LDreg/data/eur_w_ld_chr/ \
    --w-ld-chr /path/to/shared_project/gwas/LDreg/data/eur_w_ld_chr/ \
    --out /path/to/project/GWAS/LDSC/double_right_vs_single_right_knee

/path/to/shared_project/gwas/LDreg/ldsc/ldsc.py \
    --rg \
    /path/to/project/GWAS/LDSC/double_right.sumstats.gz,/path/to/project/GWAS/LDSC/single_right_hip.sumstats.gz \
    --ref-ld-chr /path/to/shared_project/gwas/LDreg/data/eur_w_ld_chr/ \
    --w-ld-chr /path/to/shared_project/gwas/LDreg/data/eur_w_ld_chr/ \
    --out /path/to/project/GWAS/LDSC/double_right_vs_single_right_hip

/path/to/shared_project/gwas/LDreg/ldsc/ldsc.py \
    --rg \
    /path/to/project/GWAS/LDSC/single_right_knee.sumstats.gz,/path/to/project/GWAS/LDSC/single_right_hip.sumstats.gz \
    --ref-ld-chr /path/to/shared_project/gwas/LDreg/data/eur_w_ld_chr/ \
    --w-ld-chr /path/to/shared_project/gwas/LDreg/data/eur_w_ld_chr/ \
    --out /path/to/project/GWAS/LDSC/single_right_knee_vs_single_right_hip


/path/to/shared_project/gwas/LDreg/ldsc/ldsc.py \
    --rg \
    /path/to/project/GWAS/LDSC/double_right.sumstats.gz,/path/to/project/GWAS/LDSC/double_right.sumstats.gz \
    --ref-ld-chr /path/to/shared_project/gwas/LDreg/data/eur_w_ld_chr/ \
    --w-ld-chr /path/to/shared_project/gwas/LDreg/data/eur_w_ld_chr/ \
    --out /path/to/project/GWAS/LDSC/double_right_rg

/path/to/shared_project/gwas/LDreg/ldsc/ldsc.py \
    --rg \
    /path/to/project/GWAS/LDSC/single_right_knee.sumstats.gz,/path/to/project/GWAS/LDSC/single_right_knee.sumstats.gz \
    --ref-ld-chr /path/to/shared_project/gwas/LDreg/data/eur_w_ld_chr/ \
    --w-ld-chr /path/to/shared_project/gwas/LDreg/data/eur_w_ld_chr/ \
    --out /path/to/project/GWAS/LDSC/single_right_knee_rg
    
/path/to/shared_project/gwas/LDreg/ldsc/ldsc.py \
    --rg \
    /path/to/project/GWAS/LDSC/single_right_hip.sumstats.gz,/path/to/project/GWAS/LDSC/single_right_hip.sumstats.gz \
    --ref-ld-chr /path/to/shared_project/gwas/LDreg/data/eur_w_ld_chr/ \
    --w-ld-chr /path/to/shared_project/gwas/LDreg/data/eur_w_ld_chr/ \
    --out /path/to/project/GWAS/LDSC/single_right_hip_rg

