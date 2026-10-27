# 脚本补充、编号与优化记录（更新至 2026-10-27）

本次核对了 `02｜SheepGTEx｜eQTL/ieQTL 图件与机制` 和 `03｜西农 scPagwas｜分析与可视化` 的近期任务记录，以及两个本地源项目的实际代码。任务记录用于定位；下面的新增模块以实际读到的源脚本为依据。

## 已直接补入或修改的代码

1. **ieGene 共享图与 E5 匹配背景**：加入采纳的共享图脚本 `01_ieQTL_cell_type_regulation/01_ieQTL_discovery_and_sharing/02_plot_ieGene_sharing.py`，以及 `supporting_analyses_not_in_requested_figures/02_regulatory_annotation/00_e5_matched_background_analysis.py`。后者以同一基因的 nominal ieQTL 变异池为背景，按等位基因频率及 cis 距离匹配；不应解读为因果检验。需核对 E5 文件的组装版本、区间坐标约定，以及 lead/背景字段定义。
2. **独立 LD 计算**：`supporting_analyses_not_in_requested_figures/01_GWAS_ieQTL_colocalization/00_compute_locus_ld.py`。从源项目 `plot_gwas_ieqtl_ld_batch.py` 的 PLINK 1.9 调用逻辑拆出，只做 LD，保留日志、工具版本、命令及 `ld.tsv`，不需要先画图。用户明确提供基因型前缀、lead ID、染色体和区域；校验 lead 在 BIM 中唯一且位于区域内，要求新输出目录。未计算的 LD 保留缺失，不能当作零。现有案例绘图仍有自己的 LD 调用，独立 LD 表尚未全部接入这些绘图模块。
3. **Seurat 分组平均表达导出**：`03_cross_species_disease_analysis/02_conserved_expression_and_enrichment/00_export_seurat_group_mean.R`。基于源项目 `export_sheep_seurat_pseudobulk.R`，显式指定分组列和 RNA layer，导出 gene × group 均值矩阵、各组细胞数和 R 会话信息。这里的均值矩阵用于表达相关，不能当作按生物学重复聚合的原始 counts 用于差异表达。不会自动把缺失 data layer 改成另一种归一化。
4. **GWAS 位点窗口预处理**：`03_cross_species_disease_analysis/01_human_GWAS_to_sheep_cells/00_prepare_500kb_gwas_loci.py` 从源项目独立出已合并位点的 ±500 kb 裁切和清单记录。它只覆盖这个明确的预处理步骤，不包含原始 GWAS 协调、坐标转换、同源映射或服务器任务编排；不能把整个 `human_icpag4_final` 目录描述成已移植。
5. **PCC 排名通用化**：`02_complex_trait_cell_mapping/03_candidate_genes_and_pathways/00_rank_candidate_contexts.py` 新增 `--contexts`，通过清单指定任意组织/性状/细胞和文件，自动调整面板数量；同时保存独立的 `ranked_candidates.tsv`。默认仍可读原论文十个背景。新增 `--include-all-genes`，允许显式关闭原论文按基因名称进行的过滤。示例清单见 `config/pcc_contexts.example.tsv`。
6. **共定位错误可追溯**：`01_run_colocalization.R` 现在拒绝覆盖已有输出，检查最小共享变异数和 P 值范围，把 `coloc.abf` 失败写入 `coloc_errors.tsv` 并发出警告。
7. **跨物种表达桥接**：组织/细胞类型桥接和哮喘案例改为显式选择 hg19/hg38 同源基因列，避免注释表带有两列时误选造成静默错配；哮喘表达汇总保留稀疏矩阵，并按基因逐列执行 rank-sum，避免把完整单细胞矩阵稠密化。

## 用法示例

```bash
python supporting_analyses_not_in_requested_figures/01_GWAS_ieQTL_colocalization/00_compute_locus_ld.py \
  --plink /path/to/plink --bfile /path/to/genotypes \
  --lead VARIANT_ID --chrom 3 --start 7000000 --end 8000000 \
  --out-dir results/new_locus_ld

Rscript 03_cross_species_disease_analysis/02_conserved_expression_and_enrichment/00_export_seurat_group_mean.R \
  /path/to/atlas.rds cell_type data results/group_mean.tsv

python 02_complex_trait_cell_mapping/03_candidate_genes_and_pathways/00_rank_candidate_contexts.py \
  --input-dir /path/to/pcc_tables --contexts config/pcc_contexts.example.tsv \
  --out-dir results/new_pcc_contexts --include-all-genes
```

路径和示例条目需替换为自己的数据。PLINK 输入要求同一组装版本，区域坐标是 1-based inclusive。

## 找到了，但暂时不直接纳入的代码

- `build_human_sample_pseudobulk.py`：可处理 10X 输入，但当前任意 H5 回退读法、重复基因名称后缀策略，以及失败后仍返回成功退出码，都应先改正再作为公开入口。
- `human_icpag4_final` 下除已移植的 locus-window 预处理外，还有坐标/对象转换和任务生成代码。这些脚本依赖特定服务器、原始 runner 和数据布局，需逐一追溯再改成独立输入，不能仅搬过来声称已复现。
- 单细胞疾病富集代码在不同分析目录里有多个分支。任务记录明确提示：保守疾病基因的富集、另一套 GSEA 和通路着色使用的统计量不同，不能把有显著结果的分支替换进论文的另一分支。

## 仍影响公共复用的主要问题

- `01_run_colocalization.R` 当前把 GWAS 按连续性状处理，并用标准误和频率推算样本量。这个推算有表型标准化假设；自己的数据应改用准确样本量。病例对照性状需另设类型、病例比例等配置。
- 当前共定位实际按每个基因的输入变异集合运行，R 脚本未独立实施 ±500-kb 过滤；窗口由上游输入决定。旧文档所称 ±500 kb 不能理解为此脚本自动裁切。
- 哮喘案例脚本当前的细胞级 `ranksums` 会把单个细胞作为检验观测值；如果细胞来自少量生物学样本，这会产生伪重复风险。公开分析中应以样本为统计单位，或明确把该结果限定为描述性探索；仅把矩阵稀疏化不会解决推断单位问题。
- 人羊表达桥接使用按细胞组求均值的表达矩阵，适合比较表达模式，但不能提供样本间变异或差异表达推断。若需要稳健的疾病差异表达，仍需恢复按供体/样本聚合的 count pseudobulk 及对应模型。
- 多个绘图脚本仍有论文案例默认基因、PP.H3/PP.H4 或固定性状集合。后续需改成显式配置，避免换了数据却沿用论文数字。
- DWLS、OmiGA、完整 scPagwas 和虚拟敲除的准确原始运行配置仍未恢复。这些缺口影响端到端复现，不能用软件通用示例冒充论文代码。
