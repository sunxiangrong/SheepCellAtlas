# 03 · 人类疾病的跨物种分析

主线：**人类疾病 GWAS 定位羊细胞背景 → 检查相关表达程序 → 比较人羊疾病基因基线表达 → 解释保守基因功能**。三个证据层分别整理，不合并成单一结论。

| 子目录 | 对应内容 | 当前代码状态 |
|---|---|---|
| `01_human_GWAS_to_sheep_cells/` | Fig. 8a 全局疾病关联；Extended Data Fig. 8 哮喘–肺和银屑病–淋巴结 | 现有 ±500 kb 位点准备脚本和哮喘案例后处理不能替代完整疾病 GWAS 到羊细胞映射。 |
| `02_conserved_expression_and_enrichment/` | Fig. 8d–g 银屑病 PCC/稳定同源基因、银屑病与 AD 角质形成细胞表达比较、保守疾病 DEG 的 GO；扩展背景的 PCC/表达与过敏原响应富集 | 有通用表达均值导出、组织/细胞类型桥接及部分结果图；指定银屑病皮肤/软角 TRS 小提琴与 UMAP、病例对应散点图、保守 DEG 和分别计算的 GO 生成链未完整纳入。 |

## 三类证据

1. **TRS／细胞关联**：疾病 GWAS 定位到组织和细胞背景（Fig. 8a–c）。
2. **稳定表达同源基因的 PCC 富集**：检查疾病相关表达程序是否偏向稳定表达同源基因（Fig. 8d）。
3. **保守表达疾病 DEG 与 GO**：用人羊疾病相关差异表达基因及功能富集解释保守基因（Fig. 8e–g）。

AD 在当前图文范围内只承担角质形成细胞表达比较和功能富集，不补写 AD–TRS。功能图按用户确认的面板称 GO，不写成 KEGG。跨物种基线表达保守不能单独证明绵羊模型再现人类疾病。

## 当前可复用脚本

- `01_human_GWAS_to_sheep_cells/00_prepare_500kb_gwas_loci.py`：从已经协调并合并的位点表准备窗口，不做 GWAS 协调或跨物种映射。
- `01_human_GWAS_to_sheep_cells/04_asthma_lung_airway_case.py`：基于已有输入整理哮喘肺/气道扩展案例，不是全局疾病关联流程。
- `02_conserved_expression_and_enrichment/00_export_seurat_group_mean.R`：导出指定 Seurat 分组的平均表达。
- `02_conserved_expression_and_enrichment/02_tissue_expression_bridge.py` 与 `03_celltype_expression_bridge.py`：计算组织或细胞类型层面的人羊表达桥接。
- `02_conserved_expression_and_enrichment/05_plot_ad_go_enrichment.R`、`06_plot_psoriasis_orthologs.py`：从既有结果表绘制部分案例结果；不生成全部上游富集统计。

复现仍需明确人类 GWAS/表达数据版本、同源基因和标签映射、细胞背景、疾病 DEG 定义及 GO 背景集。未打包的步骤在清单中如实标为缺失。
