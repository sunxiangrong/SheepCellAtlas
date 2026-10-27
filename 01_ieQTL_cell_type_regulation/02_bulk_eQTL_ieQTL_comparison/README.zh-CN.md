# bulk eQTL–ieQTL 比较

该分支比较 bulk 组织与细胞类型交互分析中的分子 QTL 证据，**不是** GWAS–ieQTL 共定位流程。

## 应分别保留的结果

1. 基因层面重叠/共享表，并明确统计单位（例如 gene × tissue 或 gene × tissue–cell-type pair）。
2. 用于共定位或区域图的匹配 bulk-QTL 和 ieQTL 位点数据。
3. 效应值相关表，明确比较的是哪些系数（例如 bulk `g1` 与 ieQTL `g1` 或 `g2`）。
4. 调控注释表，记录基因组版本、注释来源/版本、坐标规则和区间重叠条件。

源项目中有案例脚本和图件，但当前公开候选包尚未形成从源表到上述结果的完整可移植流程。最终相关系数定义、信号单位、阈值和完整案例输入需要核实后再加入。

输入要求见 [`../../docs/data-contracts.zh-CN.md`](../../docs/data-contracts.zh-CN.md)。匹配的 bulk-QTL 区域轨迹不能称作 GWAS、bulk-QTL、ieQTL 三方共定位。
