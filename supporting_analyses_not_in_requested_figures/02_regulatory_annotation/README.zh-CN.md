# 匹配背景的调控注释

`00_e5_matched_background_analysis.py` 比较 PPH3-high 与 PPH4-high ieQTL lead 变异落在组织特异 FAANG E5 增强子中的比例，并与匹配的 cis 检验背景进行比较。

## 运行顺序

1. 准备 lead 表，至少包含 `group`（`PPH3-high` 或 `PPH4-high`）、`tissue`、`cell_type`、`phenotype2` 和 `variant_id`。
2. 提供各组织/细胞类型的 nominal ieQTL 文件，以及与变异坐标版本一致的 E5 BED。
3. 先运行 `python 00_e5_matched_background_analysis.py --help`，再指定新的 `--out-dir`；记录随机种子和置换次数。
4. 一并保存 lead 注释表、背景抽样分布、汇总表和图。

背景从同一基因的 nominal ieQTL 变异中抽取，并按等位基因频率和 cis 距离匹配，具体匹配分数以脚本实现为准。这是该匹配设计下的富集比较，不是全基因组因果检验。变异坐标版本、BED 版本和区间坐标约定必须匹配。
