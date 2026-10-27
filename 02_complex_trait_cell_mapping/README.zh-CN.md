# 02 · 绵羊复杂性状的细胞定位

主线：**羊 GWAS → 细胞关联 → 代表性背景 → 候选基因与功能 → 预测网络响应**。子目录按分析顺序编号，图号用于说明对应面板，不作为目录结构。

| 子目录 | 对应内容 | 应保留的分析 | 当前代码状态 |
|---|---|---|---|
| `01_global_trait_association/` | Fig. 6a | 15 个性状 × 8 个谱系关联及多重校正 | 指定面板的完整 scPagwas 运行和汇总流程未纳入。现有谱系图脚本仅是部分结果绘图，不能据此认定 Fig. 6a 已复现。 |
| `02_representative_contexts/` | Fig. 6b–e | ALLGE3–概念胚、胸深–肺的 TRS 分布与嵌入图 | 这两个指定背景的完整结果/作图链未纳入。 |
| `03_candidate_genes_and_pathways/` | Fig. 6f | 对应背景 PCC 候选；分别报告通路活性和基因集富集（KEGG） | 现有脚本可整理选定背景 PCC 排序，不能替代候选基因与 KEGG 两类完整分析。 |
| `04_tissue_level_extensions/` | Extended Data Fig. 6a–f | TRS 全局图、ALLGE3–输卵管、OHW–PBMC landscape/violin/PCC | landscape 绘图脚本仅处理已有汇总表，且不是完整面板链；OHW–PBMC PCC 是 tissue–trait 运行层面的排序，不是非经典单核细胞内排名。 |
| `05_virtual_perturbation/` | Extended Data Fig. 6g | 耳长–直肠 C7⁺ 基质细胞，三个单基因扰动及 KCNN3 预测响应 | 原始对象、软件参数和结果生成脚本未纳入。该案例是独立探索性背景，不是胸深–肺的下游验证。 |

### 现有可用脚本

- `01_global_trait_association/01_plot_lineage_associations.py`：从汇总表绘制谱系关联；不是 15 性状 × 8 谱系完整分析。
- `03_candidate_genes_and_pathways/00_rank_candidate_contexts.py`：整理预先选定背景的 PCC 候选排序；不生成 scPagwas 结果，也不完成 KEGG 活性或富集分析。
- `04_tissue_level_extensions/00_plot_trait_landscape.py`：绘制预先计算的性状/细胞背景汇总结果；不覆盖指定扩展图的所有数据处理步骤。

主流程缺少的上游代码继续留作缺口，不用通用绘图脚本补位。数据格式和统计单位以 `../docs/data-contracts.zh-CN.md` 为准。
