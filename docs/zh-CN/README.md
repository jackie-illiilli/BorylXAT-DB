# 硼自由基 C-Cl 活化数据库

[English README](../../README.md) | [English Database Structure](../../Database_Structure.md) | [数据库结构说明](Database_Structure.md)

本仓库是硼自由基在路易斯碱活化下进行 C-Cl 键活化研究的数据生产、整理、分析与建模代码库，也是对应论文项目的配套仓库。项目目标是通过高通量量子化学计算构建过渡态数据库，并进一步开展机理分析与机器学习预测。

本项目主要包含以下几部分：

- 硼烷、路易斯碱和氯代底物反应空间的枚举
- 基于 RDKit 与 xTB 的构象生成与筛选
- 基态、限制优化、TS 搜索、单点能校正和 IRC 验证的 Gaussian 输入生成
- ASE SQLite 与 Parquet 格式的结构化数据库构建
- 描述符提取、统计分析、方法学 benchmark 和机器学习建模

## 项目范围

根据论文大纲，本研究关注的组合空间由以下三类反应物构成：

- 55 种硼自由基
- 386 种路易斯碱
- 179 种氯代底物

这些组合会经过热力学筛选后再进入过渡态计算。当前仓库内已经包含最终数据库文件：

- `boron_ccl2.db`：ASE SQLite 数据库，共 `50056` 条结构
- `boron_ccl_dataset2.parquet`：扁平化 Parquet 数据集

当前数据库中的结构类别分布为：

| 类别 | 数量 |
| --- | ---: |
| `B` | 55 |
| `LB` | 387 |
| `Cl` | 179 |
| `complex_r` | 20010 |
| `complex_p` | 20010 |
| `c_radical` | 179 |
| `ts` | 9236 |

论文叙述中的数量与当前仓库产物之间存在少量差异，通常对应不同版本的数据筛选或导出结果。

## 科学问题

这个代码库主要围绕三个问题展开：

1. 哪些硼自由基 / 路易斯碱 / 底物组合在热力学上可行？
2. 通过验证的卤原子转移过渡态具有什么样的几何与能量特征？
3. 能否把 DFT 数据转化为描述符和机器学习模型，用于预测 `ΔG‡` 及相关反应趋势？

## 工作流概览

仓库的实际工作流以 notebook 为主：

1. 从反应物表中枚举反应位点和组合。
2. 用 RDKit 生成初始三维结构，并用 xTB/CREST 做构象搜索。
3. 构建硼自由基、路易斯碱、氯代底物及 B-LB 复合物的基态结构。
4. 通过产物样式的 B-LB-Cl 几何与底物结构组合，生成 TS 初猜。
5. 进行限制优化、TS 优化、单点能校正和 IRC 验证。
6. 解析 Gaussian 输出，提取能量、几何、电荷和自旋密度等信息。
7. 组装最终的 ASE / Parquet 数据库。
8. 构建描述符，进行 benchmark 和机器学习建模。

## Notebook 说明

根目录下的主要 notebook 如下：

| Notebook | 作用 |
| --- | --- |
| `1_Calc_Reactant.ipynb` | 反应物整理、反应位点枚举、xTB 构象采样，以及硼自由基、路易斯碱、氯代底物和 B-LB 复合物的 DFT 输入生成 |
| `2_Calc_TS.ipynb` | TS 初猜构建、限制优化、TS 搜索、单点能校正、IRC 分析和 TS 汇总 |
| `3_Build_DataBase.ipynb` | 将解析后的结构和能量信息整合为 `boron_ccl2.db` 与 `boron_ccl_dataset2.parquet`，并生成分析图 |
| `4_Benchmark.ipynb` | 小规模参考数据集上的 DFT 方法 benchmark |
| `5_draw_molecule.ipynb` | 论文和 SI 所需分子图绘制 |
| `6_Modeling.ipynb` | 描述符生成、CatBoost 建模以及 OOD 分析 |

`main.py` 目前只是占位脚本，真正的工作流在 notebooks 和 `DFTStructureGenerator` 包中。

## 核心代码结构

`DFTStructureGenerator/` 中包含主要功能模块：

- `B_N_Cl.py`：反应位点识别、组合生成、DFT 作业准备、TS 结构生成、IRC 任务处理、TS 汇总导出与反应 AAM SMILES 标注
- `Build_DataBase.py`：把解析后的结构字典写入 ASE 与 Parquet 数据库，并计算 `ΔG‡` 和 `ΔG_rxn`
- `descriptor.py`：构建 B-LB 与氯代底物描述符映射，并生成 ML 可用特征矩阵
- `logfile_process.py`：Gaussian log 解析，包括能量、结构、虚频、IRC 轨迹、电荷/多重度以及错误识别
- `xtb_process.py`：xTB/CREST 构象流程与 PBS 脚本生成
- `mol_manipulation.py`：几何变换、TS/IRC 输入构建和 Gaussian 错误修复辅助
- `FormatConverter.py`：RDKit 分子、XYZ、Gaussian 输入和电荷文件之间的转换
- `draw.py`：散点图、分布图、性能图和相关性图的绘图工具
- `Tool.py`：一些小型几何与通用工具函数

## 关键数据文件

除了主数据库文件外，仓库还包含很多中间产物：

- `Data/csvs/reactants_B.csv`：硼自由基条目，共 `55` 条
- `Data/csvs/reactants_N.csv`：路易斯碱及反应位点条目，共 `415` 行
- `Data/csvs/reactants_Cl.csv`：氯代底物条目，共 `179` 条
- `Data/csvs/reactants_B_N.csv`：筛选后的 B-LB 组合，共 `20010` 条
- `Data/csvs/reactants_B_N_full.csv`：更大的组合表，共 `22825` 条
- `Data/csvs/Benchmark_Result.csv`：benchmark 结果表
- `Data/descriptor/*.pkl`：建模所用描述符映射
- `Figure/`：分析与论文配图输出目录

数据库字段细节请见 [Database_Structure.md](Database_Structure.md)。

## 电子结构计算设置

根据当前 notebook，核心 Gaussian 设置为：

- 几何优化与频率：`B3LYP/6-31G(d) + D3BJ + SMD(toluene)`
- 单点能校正：`wB97X-D/6-311+G(d,p) + SMD(toluene)`
- TS 工作流：限制优化 -> TS 优化 -> IRC 验证

代码中还区分了开壳层与闭壳层物种的自旋多重度，并支持导出带波函数输出的单点能任务，用于后续电荷分析。

## 环境依赖

当前 `pyproject.toml` 只列出了非常少的依赖。结合代码实际使用情况，建议环境至少包含：

- Python 3.12+
- `numpy`, `pandas`, `scipy`, `matplotlib`, `tqdm`
- `rdkit`
- `ase`
- `catboost`
- `scikit-learn`
- `seaborn`
- `morfeus`
- `ipykernel`
- `openpyxl` 与 `pyarrow`

此外，这个工作流还依赖外部程序：

- Gaussian
- xTB / CREST
- 支持 PBS 批处理脚本的 HPC 环境

## 使用顺序

如果想按完整流程复现，建议按以下顺序运行：

```text
1_Calc_Reactant.ipynb
2_Calc_TS.ipynb
3_Build_DataBase.ipynb
4_Benchmark.ipynb
6_Modeling.ipynb
```

如果只想使用最终数据，可直接从这些文件开始：

- `boron_ccl2.db`
- `boron_ccl_dataset2.parquet`
- `Data/descriptor/*.pkl`

## 复用说明

- `B_00001`、`LB_00001`、`Cl_00001_r`、`B_00001_LB_00001_Cl_00001` 这类命名规则贯穿整个仓库。
- 不少脚本默认沿用历史 HPC 目录结构，重新部署时通常需要先调整路径。
- 仓库中已经包含许多生成后的结果文件；若要从零完整复现，仍需要原始 Gaussian/xTB 计算环境。
- 论文叙述与当前仓库实现大体一致，但个别统计量可能因为版本更新存在轻微差异。
