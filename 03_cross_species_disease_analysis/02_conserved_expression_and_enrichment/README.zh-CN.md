# 保守表达与功能富集

## 预期顺序

1. 运行 `00_export_seurat_group_mean.R <seurat_rds> <group_column> <layer> <out_tsv>`，从已注释 Seurat 对象按指定分组（例如 tissue）导出平均表达、每组细胞数和会话信息。
2. 准备人类组织和细胞类型表达输入，记录表达层和汇总方法。
3. 按指定的一对一同源基因列筛选两物种表达并统一基因 ID。
4. 运行 `02_tissue_expression_bridge.py` 和 `03_celltype_expression_bridge.py`，生成跨物种表达比较结果。显式指定 `--human-ortholog-column`，使其与 H5AD 中的人类基因符号版本一致。
5. 定义保守疾病相关基因，使用记录明确的背景全集进行富集分析（上游富集计算目前未包含）。
6. 使用 `05_plot_ad_go_enrichment.R` 和 `06_plot_psoriasis_orthologs.py` 绘制提供的结果表。

导出的分组均值是细胞层面的平均表达，用于表达相关等桥接分析；这不是按生物学样本聚合的原始 counts，也不提供生物学重复估计。桥接脚本不会自动确定同源基因、细胞标签映射或完整疾病差异表达/富集结果。请保存每一步中间表和映射文件，便于审计。
