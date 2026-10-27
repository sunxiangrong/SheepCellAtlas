# 01 · 细胞类型交互 eQTL（ieQTL）

本模块回答：遗传效应是否随 bulk 样本中的细胞类型比例变化，以及这些交互信号与同组织 bulk eQTL 有何关系。ieQTL 来自 bulk 表达、基因型和推断细胞比例的交互模型，不称为 direct single-cell eQTL。

## 分析顺序

| 步骤 | 内容 | 代码状态 |
|---|---|---|
| 01 | 用 SheepCA 参考从 bulk RNA-seq 推断细胞比例（DWLS） | 论文实际运行脚本和参数未纳入。 |
| 02 | 对齐基因型、表达、细胞比例和协变量，拟合 genotype × cell fraction 模型（OmiGA） | 实际命令、版本、协变量和阈值待补。 |
| 03 | 汇总 ieQTL/ieGene 检出及组织共享 | `01_ieQTL_discovery_and_sharing/` 含共享表准备和结果图；不负责模型拟合。 |
| 04 | 比较 bulk eQTL 与 ieQTL 的共享和效应关系 | `02_bulk_eQTL_ieQTL_comparison/` 当前记录流程边界；完整可复用脚本仍缺。 |
| 05 | 展示与主图不重复的代表性交互案例 | 尚未确认并纳入符合主图/补图分工的完整案例脚本。 |

## 已核实的结果口径

目前整理记录为 980 条组织–细胞类型–基因 ieQTL 记录、879 个唯一 ieGene；其中 336 条对应组织的 bulk eGene 未达到 `qval_g1 < 0.05`。应表述为“未达到 bulk eGene 显著标准”，不能推断 bulk 中没有该基因或任何调控信号。

主图和补图案例需互补：补图只加入未在主图出现的 PPH3-high/PPH4-low 与 PPH4-high 案例；不重复 GAPVD1、ECD、WASHC2C，也不把 STXBP1 GWAS 案例或 GWAS–ieQTL 共定位、E5 调控注释混入本模块。相应分析文件暂放在仓库根目录的 `supporting_analyses_not_in_requested_figures/`。

## 当前脚本

- `01_ieQTL_discovery_and_sharing/00_prepare_sharing_tables.py`：从已有显著 ieQTL 与组织配对汇总表整理共享统计。
- `01_ieQTL_discovery_and_sharing/01_plot_detectability.py`：从预先汇总表绘制检出情况。
- `01_ieQTL_discovery_and_sharing/02_plot_ieGene_sharing.py`：绘制 ieGene 组织共享关系。
- `02_bulk_eQTL_ieQTL_comparison/README.zh-CN.md`：记录比较目标、应使用的统计单位及待补脚本；不是已完成分析。

公开复现前仍需补上 DWLS、OmiGA 和 bulk eQTL–ieQTL 完整比较链，并提供每一步的数据字段、筛选规则和版本信息。
