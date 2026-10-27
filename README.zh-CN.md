# SheepCA 遗传调控与性状分析

代码按论文中的三个分析问题组织，目录编号表示分析模块和步骤顺序，不按图号建文件夹。每个模块先说明分析目的，再按数据处理、统计分析和结果整理排列脚本；没有纳入或尚未恢复的步骤会明确标注。

```text
01_ieQTL_cell_type_regulation/       ieQTL 发现及 bulk eQTL 对照
02_complex_trait_cell_mapping/       绵羊复杂性状细胞定位
03_cross_species_disease_analysis/   人类疾病的跨物种分析
supporting_analyses_not_in_requested_figures/  本次三条主线之外的分析
```

目录组织参考 [FarmGTEx Cattle Cell Atlas pipeline](https://github.com/FarmGTEx/CattleCellAtlas_pipeline_V0)：按主题分模块、模块内按步骤编号，并用 README 说明用途和输入输出。这里仅借鉴组织方式，不照搬其分析内容。

## 主线与当前状态

| 模块 | 论文主线 | 代码现状 |
|---|---|---|
| 01 | bulk 表达、基因型和推断细胞比例 → genotype × cell-fraction interaction ieQTL → 与 bulk eQTL 比较、展示代表案例（Fig. 5） | 含 ieGene 共享表准备和部分结果图；DWLS、OmiGA 实际运行及完整 bulk 对照仍待补。不是 direct single-cell eQTL。 |
| 02 | 羊 GWAS → 谱系/组织/细胞定位 → PCC 候选与 KEGG → 虚拟扰动（Fig. 6 及 Extended Data Fig. 6） | 含部分结果整理脚本；指定面板的 15×8 全局关联、代表性背景、KEGG 和虚拟扰动流程尚未完整纳入。 |
| 03 | 人类疾病 GWAS → 羊细胞背景 → 稳定表达同源基因 → 保守疾病 DEG 与 GO（Fig. 8 及 Extended Data Fig. 8） | 含表达桥接和部分案例后处理；疾病映射、指定银屑病定位、疾病 DEG/GO 完整生成流程仍有缺口。 |

脚本只声明其实际输入和输出，不用绘图代码代替缺失的上游统计分析。Fig. 6 的 OHW–PBMC PCC 是 tissue–trait 分析层面的排序；耳长–直肠 C7⁺ 基质细胞虚拟扰动是独立的探索性背景。Fig. 8 中 TRS/细胞关联、稳定表达同源基因 PCC 富集、保守疾病 DEG 与 GO 是三类不同证据；AD 仅按现有证据用于表达比较和功能分析，不补写 AD–TRS。

本仓库是公开代码候选稿，不包含受限基因型、表达矩阵、单细胞对象或个体元数据。输入格式与尚缺步骤见 `docs/data-contracts.zh-CN.md` 和 `PUBLICATION_MANIFEST.zh-CN.tsv`。环境示例见 `environment.yml`；公开前按 `docs/release-checklist.md` 检查。
