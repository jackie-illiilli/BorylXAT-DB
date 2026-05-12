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

这些组合会经过热力学筛选后再进入过渡态计算。当前生成的最终数据库文件为：

- `boron_ccl.db`：ASE SQLite 数据库，共 `50057` 条结构
- `boron_ccl_dataset.parquet`：扁平化 Parquet 数据集

当前 ASE 数据库中的结构类别分布为：

| 类别 | 数量 |
| --- | ---: |
| `B` | 55 |
| `LB` | 387 |
| `Cl` | 179 |
| `complex_r` | 20010 |
| `complex_p` | 20010 |
| `c_radical` | 179 |
| `ts` | 9237 |

正文中的反应空间按 386 种 Lewis 碱统计，因为其中一个单体 Lewis 碱条目（`LB_00623`）与任意硼自由基形成 B-LB 复合物时在热力学上都不够稳定。数据库为了保留来源信息，仍在 `LB` 类别中保留该分子；但它不会出现在筛选后的 B-LB 复合物集合或 TS 反应条目中。

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
| `3_Build_DataBase.ipynb` | 将解析后的结构和能量信息整合为 `boron_ccl.db` 与 `boron_ccl_dataset.parquet`；同时保留可由审稿人直接运行的数据库检查与 Figure 1 类汇总图 |
| `4_Benchmark.ipynb` | DFT 方法 benchmark；已区分依赖原始 Gaussian 结果的部分和可直接基于已整理数据运行的结果分析部分 |
| `5_Modeling.ipynb` | 描述符生成、已提取描述符加载、CatBoost 建模、验证图与 OOD 分析 |
| `6_Draw_Figures.ipynb` | 集中的论文/SI 作图 notebook；整合了原先分散在数据库构建与分子作图流程中的绘图代码 |

`main.py` 目前只是占位脚本，真正的工作流在 notebooks 和 `DFTStructureGenerator` 包中。

## Notebook 润色与审稿运行标记

当前 notebooks 已按论文审稿和复用场景做过整理，同时保留了完整计算来源追溯所需的生产流程。

| Notebook | 已完成整理 | 审稿可用性 |
| --- | --- | --- |
| `1_Calc_Reactant.ipynb` | 重写不合时宜的 markdown 批注，收紧反应物枚举/优化相关注释，删除未使用的 import 和遗留变量。 | 记录生产级反应物优化流程；完整重跑需要 Gaussian、xTB/CREST 和历史计算目录。 |
| `2_Calc_TS.ipynb` | 清理未使用变量和 import，将 `complete_target()` 移入 `DFTStructureGenerator/borane_xat_workflow.py`，并把 target CSV 补全部分移动到末尾作为可选补充/导出步骤。 | 记录 TS 初猜、限制优化、TS 优化、SPE、IRC 和汇总导出流程；完整重跑需要原始 Gaussian/xTB 文件夹。 |
| `3_Build_DataBase.ipynb` | 添加审稿运行标记，区分原始 Gaussian 解析和已发布数据库检查；论文作图相关功能集中到 `6_Draw_Figures.ipynb`。 | 审稿人可直接检查已发布的 `boron_ccl.db` 和 `boron_ccl_dataset.parquet`。 |
| `4_Benchmark.ipynb` | 清理 benchmark 注释，重新组织 setup、方法定义和结果分析，并把原始结构/输入收集部分和已提交结果分析部分分开。 | 审稿人可直接加载 `Data/csvs/Benchmark_Result.csv` 并重新生成 benchmark 汇总图。 |
| `5_Modeling.ipynb` | 润色 import 和注释，拆分描述符生成与描述符加载，把描述符生成标记为可选，并补上 Figure 6B 的 `plt.savefig(...)`。 | 审稿人可使用 `Data/descriptor/` 中预提取的描述符，跳过耗时描述符生成后直接进行模型验证和绘图。 |
| `6_Draw_Figures.ipynb` | 把原先分散在数据库构建和分子作图流程中的论文/SI 作图代码整合到一个 notebook；原 Figure 4 标记和图片文件名已改为 Figure 5。 | 当前推荐的审稿作图入口，可从仓库内已提交的数据产品重新生成论文和 SI 图片。 |

历史上的 `B_N_Cl` 命名更适合改称为 `borane_xat_workflow`：当前模块不只是静态的 B/N/Cl 组合表，而是覆盖硼自由基、Lewis 碱复合物、氯代底物、XAT 过渡态生成、汇总导出和反应标注的完整工作流。

Notebook 中使用的标记含义如下：

| 标记 | 含义 |
| --- | --- |
| `[REVIEWER-RUNNABLE]` | 基于仓库中已提交的数据文件即可运行的审稿分析单元 |
| `[RAW-GAUSSIAN/E:/work]` | 依赖历史 Gaussian/xTB 输出目录的来源追溯或完整复现单元，尤其是 `E:/work` 路径 |
| `[OPTIONAL-DESCRIPTOR-GENERATION]` | 计算量较大的描述符生成单元；仓库已提供预计算描述符文件，便于快速审稿 |

## 核心代码结构

`DFTStructureGenerator/` 中包含主要功能模块：

- `borane_xat_workflow.py`：反应位点识别、组合生成、DFT 作业准备、TS 结构生成、IRC 任务处理、TS 汇总导出与反应 AAM SMILES 标注
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

## 图片输出

当前推荐从 `6_Draw_Figures.ipynb` 统一导出论文和 SI 图片。清理过程中，原先 Figure 4 相关的 notebook 标记和图片文件名已经统一改为 Figure 5。已提交和 notebook 可重新导出的关键输出包括：

- `Figure/` 中的 Figure 1 与 Figure 3 汇总图
- Figure 5 系列：`Figure5_reaction_landscape.png`、`Figure5A_TS_type_distribution.png`、`Figure5B_BCl_BDE_by_B_type.png`、`Figure5B_BCl_BDFE_by_LB_type.png`、`Figure5B_CCl_BDFE_by_hybridization.png`、`Figure5C_BEP_residual.png` 和 `Figure5D_TS_geometry.png`
- Benchmark 输出：`FigureS17_Benchmark_MAE_R2_combined.png`
- 重新运行 `5_Modeling.ipynb` 时导出的模型图：`Figure6B_model_validation.png`、`Figure6C_model_feature_importance.png` 和 `FigureS21_descriptor_correlation_map.png`
- `Figure/FigureInSI/` 下的 SI 分子网格图与分布图

图片清理时也同步把已经生成的 Figure 4 图片文件名改为 Figure 5，因此当前 `Figure/` 目录和 `6_Draw_Figures.ipynb` 中的论文编号保持一致。

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
5_Modeling.ipynb
6_Draw_Figures.ipynb
```

如果审稿人不需要重新连接原始 Gaussian 计算目录，建议使用以下路径：

1. 在 `3_Build_DataBase.ipynb` 中只运行基于已提交数据库的检查和汇总分析单元。
2. 在 `4_Benchmark.ipynb` 中从 benchmark 结果表加载之后的分析部分开始运行。
3. 在 `5_Modeling.ipynb` 中从预计算描述符加载部分开始运行。
4. 使用 `6_Draw_Figures.ipynb` 从仓库内的数据产品重新导出论文和 SI 图片。

如果只想使用最终数据，可直接从这些文件开始：

- `boron_ccl.db`
- `boron_ccl_dataset.parquet`
- `Data/descriptor/*.pkl`

## 审稿代码可用性

仓库已经按“可直接审稿分析”和“完整计算来源复现”两种场景组织。审稿人无需访问原始私有计算目录，也可以检查发布数据集并重新生成主要分析图。

审稿人可直接运行的入口：

- `3_Build_DataBase.ipynb`：基于已提交 ASE/Parquet 数据库的检查与汇总分析
- `4_Benchmark.ipynb`：基于 `Data/csvs/Benchmark_Result.csv` 的 benchmark 方法分析
- `5_Modeling.ipynb`：基于已提交描述符文件的描述符加载、CatBoost 模型验证、特征重要性、OOD 分析和 Figure 6 导出
- `6_Draw_Figures.ipynb`：基于已提交 CSV/数据库产品的论文和 SI 图片重新生成

完整来源复现入口：

- `1_Calc_Reactant.ipynb`：反应物枚举、xTB 构象搜索输入、Gaussian 优化/SPE 输入和反应物表解析
- `2_Calc_TS.ipynb`：TS 初猜、限制优化、TS 搜索、SPE 校正、IRC 验证和 TS CSV 补全
- `3_Build_DataBase.ipynb` 与 `4_Benchmark.ipynb` 中的 raw section：解析原始 Gaussian log、收集 benchmark 结构和结果

完整来源复现需要 Gaussian、xTB/CREST、PBS 风格 HPC 环境，并需要根据本地环境调整历史 `E:/work` 计算路径。标准审稿建议运行 `[REVIEWER-RUNNABLE]` 单元，跳过 `[RAW-GAUSSIAN/E:/work]` 单元。在 `5_Modeling.ipynb` 中，除非需要重新生成描述符，否则建议跳过 `[OPTIONAL-DESCRIPTOR-GENERATION]`，直接使用 `Data/descriptor/` 中已提取好的描述符以节省时间。

## 复用说明

- `B_00001`、`LB_00001`、`Cl_00001_r`、`B_00001_LB_00001_Cl_00001` 这类命名规则贯穿整个仓库。
- 不少脚本默认沿用历史 HPC 目录结构，重新部署时通常需要先调整路径。
- 仓库中已经包含许多生成后的结果文件；若要从零完整复现，仍需要原始 Gaussian/xTB 计算环境。
- 论文审稿时建议优先运行带 `[REVIEWER-RUNNABLE]` 标记的单元，并优先使用仓库内已提交的数据库、benchmark 表和描述符 pickle 文件；只有需要完整计算来源追溯时才运行 `[RAW-GAUSSIAN/E:/work]` 单元。
- 正文中的 Lewis 碱数量为 386，而数据库保留 387 个单体 `LB` 条目；额外的 `LB_00623` 仅作为来源记录保留，不参与稳定 B-LB 复合物和 TS 条目。
