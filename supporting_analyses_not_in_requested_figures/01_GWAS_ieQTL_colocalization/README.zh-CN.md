# GWAS–ieQTL 共定位

## 预期顺序

1. 统一 GWAS 与 ieQTL 的变异 ID、等位基因、基因组版本、等位基因频率和 cis 窗口。
2. 如需预先计算 LD，先运行 `00_compute_locus_ld.py --help`；该步骤独立输出 LD 和运行记录。再运行 `01_run_colocalization.R`，计算每个性状/细胞背景的共定位。
3. 准备候选位点清单，运行 `02_plot_candidate_loci.py` 绘制选定位点。
4. 运行 `03_plot_ear_length_case.py` 复现耳长 `STXBP1` 案例；需要提供用于计算 LD 的 PLINK 基因型文件。

R 脚本输入字段见 [`../../docs/data-contracts.zh-CN.md`](../../docs/data-contracts.zh-CN.md)。共定位区域由上游输入决定，R 脚本不会自动裁切 ±500 kb；最小共享变异数由命令参数指定。输出和错误日志保存在新目录。绘图脚本只负责展示结果，不能替代输入协调或共定位统计。
