# ---------------------------------------------------------------------------
# DXA-aging-pipeline / 02_epi_outcomes
# 单部位（膝或髋）年龄加速的 scPagwas 分析
# 注意: 文中 /path/to/... 为占位符，运行前请按本机实际路径替换（见 config/）。
# ---------------------------------------------------------------------------
library(scPagwas)
library(Seurat)
library(data.table)
library(dplyr)
rm(list = ls())
setwd("/path/to/epi/result/F4_scPagwas")
source('/path/to/shared_project/singlecell/code/scPagwas_main.R')


start_time <- Sys.time()
load('/path/to/shared_project/singlecell/data/scPagwas_val/genes.by.pathway_kegg.RData')
load('/path/to/shared_project/singlecell/data/scPagwas_val/chrom_ld_new.RData')
load('/path/to/shared_project/singlecell/data/scPagwas_val/block_annotation_hg37_new.RData')


Pagwas <- scPagwas_main(Pagwas =NULL,
                        gwas_data='/path/to/epi/data/GWAS_scPagwas/single_right_hip_scpagwas.txt', # GWAS 数据    
                        Single_data ='/path/to/epi/data/GSE317068/GSE317068_merge_cell.rds', # Seurat 格式的单细胞数据 
                        #Single_data =sce.all, # Seurat 格式的单细胞数据    
                        output.prefix="single_right_hip", #输出结果文件的前缀    
                        output.dirs="single_right_hip", 
                        Pathway_list=Genes_by_pathway_kegg, # 通路数据    
                        n.cores=50, # 多核计算,根据自己的实际情况修改    
                        assay="RNA",    
                        singlecell=T,    
                        iters_singlecell = 100, #背景校正单细胞 p 值计算的迭代次数    
                        celltype=T,    
                        block_annotation = block_annotation,    
                        chrom_ld = chrom_ld)

save(Pagwas,file="single_right_hip/single_right_hip_scPagwas.RData")
end_time <- Sys.time()
cat("Total running time:\n")
print(end_time - start_time)


rm(list = ls())
start_time <- Sys.time()
load('/path/to/shared_project/singlecell/data/scPagwas_val/genes.by.pathway_kegg.RData')
load('/path/to/shared_project/singlecell/data/scPagwas_val/chrom_ld_new.RData')
load('/path/to/shared_project/singlecell/data/scPagwas_val/block_annotation_hg37_new.RData')


Pagwas <- scPagwas_main(Pagwas =NULL,
                        gwas_data='/path/to/epi/data/GWAS_scPagwas/single_right_knee_scpagwas.txt', # GWAS 数据    
                        Single_data ='/path/to/epi/data/GSE317068/GSE317068_merge_cell.rds', # Seurat 格式的单细胞数据 
                        #Single_data =sce.all, # Seurat 格式的单细胞数据    
                        output.prefix="single_right_knee", #输出结果文件的前缀    
                        output.dirs="single_right_knee", 
                        Pathway_list=Genes_by_pathway_kegg, # 通路数据    
                        n.cores=50, # 多核计算,根据自己的实际情况修改    
                        assay="RNA",    
                        singlecell=T,    
                        iters_singlecell = 100, #背景校正单细胞 p 值计算的迭代次数    
                        celltype=T,    
                        block_annotation = block_annotation,    
                        chrom_ld = chrom_ld)

save(Pagwas,file="single_right_knee/single_right_knee_scPagwas.RData")
end_time <- Sys.time()
cat("Total running time:\n")
print(end_time - start_time)