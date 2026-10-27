# ieQTL 发现结果与共享

## 预期顺序

1. 以 SheepCA 单细胞参考估计 bulk 样本的细胞比例。论文使用 DWLS，但准确的运行命令和配置还未纳入本仓库。
2. 对齐基因型、表达、细胞比例、样本协变量和组织/细胞标签；使用论文所用 OmiGA 版本拟合交互模型。该上游分析尚未纳入。
3. 按研究阈值导出显著 ieQTL 合并表和 tissue–cell-type pair 汇总表。
上面完整流程的第 1–3 步仍未打包。下列三个脚本从 `00` 起仅在本子目录内编号，输入是这些上游步骤生成的结果：

1. `00_prepare_sharing_tables.py --summary ... --merged ... --out-dir ...`：生成变异层面和 ieGene 层面的组织共享表。
2. `01_plot_detectability.py --help`：绘制预先汇总的 ieQTL 检出情况。
3. `02_plot_ieGene_sharing.py --help`：绘制 ieGene 在组织间及组织内细胞类型间的共享关系；不会重新拟合 ieQTL。

## 表格字段

合并 ieQTL 表为 gzip 压缩 TSV，至少包含 `tissue`、`cell_type`、`pair`、`pheno_id`、`variant_id`。pair 汇总表需要 `tissue`、`cell_type`、`pair`。检出率绘图表需要 `tissue`、`cell_type`、`pair`、`effective_fraction_sample_size` 和 `iegene_detection_proportion`（0–1）。

共享脚本只负责下游汇总，不负责估计细胞比例、拟合 OmiGA 或决定显著性阈值。请在数据清单中记录这些上游参数。
