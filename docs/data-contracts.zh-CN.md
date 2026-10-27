# 输入与输出约定

以下说明只针对当前已整理的可移植脚本，不代表缺失的上游步骤已经实现。除非脚本另有说明，表格使用 UTF-8、带列名的 TSV。受控数据和个体级信息不得提交到 Git。

## 01 · ieQTL 遗传调控

### 共享表整理

`01_ieQTL_cell_type_regulation/01_ieQTL_discovery_and_sharing/00_prepare_sharing_tables.py` 接受 pair 汇总 TSV（列：`tissue`、`cell_type`、`pair`）和合并后的显著 ieQTL gzip TSV（列至少包括 `tissue`、`cell_type`、`pair`、`pheno_id`、`variant_id`），输出变异层面及 ieGene 层面的共享统计和共享定义。

`01_ieQTL_cell_type_regulation/01_ieQTL_discovery_and_sharing/01_plot_detectability.py` 需要 `tissue`、`cell_type`、`pair`、`effective_fraction_sample_size`、`iegene_detection_proportion`（0–1）。这些是下游结果；DWLS 细胞比例估计和 OmiGA 拟合尚未纳入。

### 共定位

`supporting_analyses_not_in_requested_figures/01_GWAS_ieQTL_colocalization/01_run_colocalization.R` 需要每个 trait 目录下的 `chr1.tsv` 至 `chr26.tsv`。GWAS 列为 `SNP`、`Freq`、`b`、`se`、`p`；ieQTL 列为 `pheno_id`、`variant_id`、`af`、`beta_g2`、`beta_se_g2`、`pval_g2`。

`02_plot_candidate_loci.py` 需要候选清单列 `trait`、`tissue`、`cell_type`、`gene`、`lead_snp`、`pph3`、`pph4`，以及区域 GWAS/ieQTL 文件。`03_plot_ear_length_case.py` 还需要 PLINK 基因型前缀和可执行程序计算 LD。使用前确认基因组版本、染色体命名、等位基因和变异 ID 一致。

bulk eQTL–ieQTL 比较是另一个分支，当前尚无完整可移植脚本序列，不要把 GWAS–ieQTL 结果当作 bulk-QTL 比较。

## 02 · 复杂性状细胞定位

结果图脚本需要如下结构的结果目录：`<tissue>/<trait>/Merged_celltype_summary.csv`，以及包含 `cell_type`、`cell_lineage` 的谱系映射 TSV。汇总表需要 `celltype`、`mean_TRS`、`merged_BG_adjp`、`n_cells`。

PCC 背景排序脚本需要 `pcc_top1000_<context>.tsv` 格式的文件，列为 `gene` 和 `PCC_clean`。上游数据清单还要说明性状列表、表型单位、基因组版本、scPagwas 版本、阈值和映射配置；这些上游步骤当前未打包。

## 03 · 跨物种疾病分析

表达桥接和哮喘案例脚本需要一对一同源基因记录，包含 `Human_homology_type`、`Sheep_gene_name`、`Human_gene_name_hg19`，以及可选的 `Human_gene_name_hg38`。运行时用 `--human-ortholog-column` 显式选择与人类对象基因符号相符的列，避免静默使用错误版本。Seurat 导出脚本按用户指定的分组列和 RNA layer 导出均值表达、组内细胞数及会话信息；这是细胞均值，不是按生物学样本汇总的原始 counts。人类 `.h5ad` 需具有指定分组列并明确使用的表达层。细胞标签需事先完成协调。

GWAS 位点窗口预处理脚本只接受已整理的性状 GWAS/COJO 表，输出合并位点周围窗口及处理清单；不负责原始 summary statistics 协调、组装版本转换或等位基因 harmonization。运行前需确认输入的坐标组装版本和坐标约定一致。

银屑病汇总需要 `candidate_pcc_with_stable_ortholog_flag.tsv.gz` 和 `candidate_stable_fraction_enrichment_summary.tsv`。皮炎 GO 绘图脚本需要预先生成的保守差异基因表和 GO 富集表。哮喘案例需要整理后的人类疾病关联结果、羊肺 pseudobulk、人类气道对象和同源基因表。

每种疾病都应记录 GWAS 来源/版本、组装版本、坐标转换、表达层、pseudobulk 分组、同源基因版本、细胞标签对应表、富集背景和阈值。疾病 GWAS 到细胞映射及富集结果生成尚未完整打包。
