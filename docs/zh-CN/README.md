# BorylXAT-DB: 硼自由基介导 C-Cl 原子转移过渡态数据库

[English README](../../README.md) | [English Database Structure](../../Database_Structure.md) | [数据库结构说明](Database_Structure.md)

本仓库包含 **BorylXAT-DB** 的工作流代码、已发布数据接口以及论文和审稿修订阶段使用的分析 notebook。项目聚焦于路易斯碱配位硼自由基介导的 C-Cl 原子转移反应，并提供从结构生成、过渡态搜索到描述符建模和统计分析的完整记录。

公开代码仓库：<https://github.com/jackie-illiilli/BorylXAT-DB>  
Zenodo 数据归档：<https://doi.org/10.5281/zenodo.20134535>

## 仓库覆盖内容

本项目主要包括：

- 硼自由基、Lewis 碱和氯代底物反应空间的枚举
- 基于 RDKit 与 xTB 的构象生成与筛选
- Gaussian 基态优化、限制优化、TS 搜索、单点能校正和 IRC 验证输入生成
- ASE SQLite 与 Parquet 结构化数据库构建
- 描述符提取、benchmark、BEP 分析和机器学习建模
- 审稿人可直接运行的主文与修订分析 notebook

论文设计空间包含：

- `55` 个硼自由基
- `386` 个 Lewis 碱
- `179` 个氯代底物

当前已发布的数据库文件为：

- `BorylXAT-DB.db`：ASE SQLite 结构数据库
- `BorylXAT-DB.parquet`：扁平化 Parquet 导出
- `BorylXAT-DB_qh_update.db`：QHARM 热化学更新数据库，供 QHARM 修订 notebook 使用

当前 ASE 数据库中的类别分布如下：

| 类别 | 数量 |
| --- | ---: |
| `B` | 55 |
| `LB` | 387 |
| `Cl` | 179 |
| `complex_r` | 20010 |
| `complex_p` | 20010 |
| `c_radical` | 179 |
| `ts` | 8980 |

数据库中的单体 `LB` 数量是 `387`，比正文中的 `386` 多一个，是因为 `LB_00623` 仅作为来源记录保留，没有进入正文筛选后稳定的 B-LB 复合物和 TS 数据集中。

## 依赖文件

为回应审稿阶段对可复现性的要求，仓库现已提供三种机器可读依赖文件：

- `pyproject.toml`：包元数据与核心 Python 依赖
- `requirements.txt`：适合 `pip` / `venv` 的审稿与分析环境
- `environment.yml`：适合 `conda` 的完整环境，尤其适用于安装 `rdkit`

公开 notebook 和模块使用到的主要 Python 包包括：

- `ase`
- `catboost`
- `ipykernel`
- `jupyterlab`
- `matplotlib`
- `morfeus-ml`
- `numpy`
- `openpyxl`
- `pandas`
- `pyarrow`
- `rdkit`
- `scikit-learn`
- `scipy`
- `seaborn`
- `tqdm`
- `xgboost`，用于部分修订版 baseline notebook

只有在完整重跑原始量化计算时才需要以下外部软件：

- Gaussian
- xTB / CREST
- 支持 PBS 批处理脚本的 HPC 环境

## 环境安装

### 方案 1：pip / venv

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Windows PowerShell：

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

### 方案 2：conda

```bash
conda env create -f environment.yml
conda activate borylxat-db
```

### 方案 3：uv

```bash
uv sync --extra revision
```

如果需要运行带 `xgboost` 的修订版 baseline notebook，建议带上 `revision` extra。

## 关键数据文件及用途

当前工作区已经包含主要发布数据文件。如果你是从干净仓库重新克隆且本地缺失这些大文件，可以从 Zenodo 下载后放到仓库根目录。

| 文件 | 用途 | 主要对应 notebook |
| --- | --- | --- |
| `BorylXAT-DB.db` | 已发布 ASE 结构数据库 | `3_Build_DataBase.ipynb`、`5_Modeling.ipynb`、作图 notebook |
| `BorylXAT-DB.parquet` | 扁平化结构/反应表 | 表格检查和下游分析 |
| `BorylXAT-DB_qh_update.db` | QHARM 更新数据库 | QHARM 修订 notebook |
| `Data/TS/Borane_all.csv` | 公开反应主表 | 建模与 BEP 分析 |
| `Data/TS/result_filter_extended.csv` | 扩展筛选结果表 | 热力学筛选与覆盖率修订 notebook |
| `Data/csvs/Benchmark.csv` | 当前 benchmark 汇总表 | `4_Benchmark.ipynb` |
| `Data/csvs/Experiment.csv` | 实验验证输入表 | `revision_JACS_Experiment.ipynb` |
| `Data/descriptor/*.pkl` | 预计算描述符文件 | `5_Modeling.ipynb`、`5_Modeling_xTB.ipynb`、修订建模 notebook |

如需下载主数据库文件，可使用：

```bash
curl -L -o BorylXAT-DB.db "https://zenodo.org/records/20134535/files/BorylXAT-DB.db?download=1"
curl -L -o BorylXAT-DB.parquet "https://zenodo.org/records/20134535/files/BorylXAT-DB.parquet?download=1"
```

Windows PowerShell：

```powershell
Invoke-WebRequest -Uri "https://zenodo.org/records/20134535/files/BorylXAT-DB.db?download=1" -OutFile "BorylXAT-DB.db"
Invoke-WebRequest -Uri "https://zenodo.org/records/20134535/files/BorylXAT-DB.parquet?download=1" -OutFile "BorylXAT-DB.parquet"
```

## Notebook 说明

### 主 notebook

| Notebook | 作用 | 是否可直接基于已发布文件运行 | 主要输入 |
| --- | --- | --- | --- |
| `1_Calc_Reactant.ipynb` | 反应物准备、枚举、xTB 构象与 Gaussian 输入生成 | 否，主要用于来源追溯 | 原始计算目录、Gaussian、xTB/CREST |
| `2_Calc_TS.ipynb` | TS 初猜、限制优化、TS 搜索、IRC 工作流 | 否，主要用于来源追溯 | 原始计算目录、Gaussian |
| `3_Build_DataBase.ipynb` | 数据库构建和已发布数据库检查 | 检查部分可直接运行 | `BorylXAT-DB.db`、`BorylXAT-DB.parquet` |
| `4_Benchmark.ipynb` | benchmark 分析 | 是 | `Data/csvs/Benchmark.csv` |
| `5_Modeling.ipynb` | DFT 描述符 CatBoost 建模 | 是 | `Data/TS/Borane_all.csv` 与 `Data/descriptor/` |
| `5_Modeling_xTB.ipynb` | xTB 描述符建模 | 是 | `Data/TS/Borane_all.csv`、`BNdes_xtb.pkl`、`Cldes_xtb.pkl` |
| `6_Draw_Figures.ipynb` | 论文和 SI 图片统一导出 | 是 | 已发布 CSV、descriptor、database 和输出文件 |

### 审稿修订新增 notebook

| Notebook | 用途 | 备注 |
| --- | --- | --- |
| `revision_ml_baseline_ablation_bep_dependence.ipynb` | ML baseline、ablation 与 BEP 依赖分析 | `xgboost` 为可选但推荐 |
| `revision_thermodynamic_filter_BEP_ML.ipynb` | 热力学筛选辅助集分析 | 使用已发布和修订 CSV 输出 |
| `revision_JACS_Experiment.ipynb` | 外部实验趋势验证 | 使用 `Data/csvs/Experiment.csv` 和 `output/jacs_experiment/` |
| `revision_qharm_barrier_analysis.ipynb` | RRHO 与 QHARM 能垒比较 | 需要 `BorylXAT-DB_qh_update.db` |
| `revision_qharm_modeling.ipynb` | QHARM 更新后的建模分析 | 优先使用 QHARM 描述符 |
| `revision_basis_set_bond_length_comparison.ipynb` | 6-31G(d) 与 6-31+G(d) 几何敏感性比较 | 使用已发布修订输出表 |
| `7_QHARM_Thermodynamic_Filter_Coverage.ipynb` | QHARM favorable-domain 覆盖率分析 | 使用修订 CSV 输出 |

Notebook 中的审稿标签含义如下：

| 标签 | 含义 |
| --- | --- |
| `[REVIEWER-RUNNABLE]` | 可直接基于仓库和 Zenodo 发布文件运行 |
| `[RAW-GAUSSIAN/E:/work]` | 依赖历史 Gaussian/xTB 原始工作目录 |
| `[OPTIONAL-DESCRIPTOR-GENERATION]` | 计算量较大，通常可跳过，直接使用仓库提供的预计算描述符 |

## 仓库结构

```text
.
|-- DFTStructureGenerator/     可复用 Python 模块，负责解析、描述符、数据库和绘图
|-- Data/
|   |-- ChemDraw/              反应物库 ChemDraw 源文件
|   |-- csvs/                  中间和已发布 CSV 汇总表
|   |-- descriptor/            已发布描述符 pickle 文件
|   `-- TS/                    已发布 TS 级 CSV 表和坐标
|-- Figure/                    论文与 SI 图片导出结果
|-- output/                    修订分析输出与已执行 notebook
|-- 1_Calc_Reactant.ipynb
|-- 2_Calc_TS.ipynb
|-- 3_Build_DataBase.ipynb
|-- 4_Benchmark.ipynb
|-- 5_Modeling.ipynb
|-- 5_Modeling_xTB.ipynb
|-- 6_Draw_Figures.ipynb
|-- revision_*.ipynb
|-- 7_QHARM_Thermodynamic_Filter_Coverage.ipynb
`-- Database_Structure.md
```

## 核心 Python 模块

`DFTStructureGenerator/` 中的主要模块包括：

- `borane_xat_workflow.py`：反应位点识别、组合生成、DFT 输入准备、TS 结构生成和汇总导出
- `Build_DataBase.py`：构建 `BorylXAT-DB.db` 与 `BorylXAT-DB.parquet`
- `descriptor.py`：生成描述符映射并组装 ML 特征
- `thermochemistry.py`：RRHO / QHARM 数据库读取工具，供修订 notebook 使用
- `logfile_process.py`：Gaussian log 解析
- `xtb_process.py`：xTB / CREST 构象工作流辅助函数
- `mol_manipulation.py`：几何变换与 Gaussian 输入辅助
- `FormatConverter.py`：RDKit、XYZ 与 Gaussian 相关格式转换
- `draw.py`：绘图辅助函数
- `project_paths.py`：仓库相对路径与原始计算目录环境变量接口

## 原始计算路径说明

历史上的原始计算工作目录位于 `E:/work/...`。仓库中的审稿可运行部分已经尽量避免强依赖该路径，但如果要完整重跑来源工作流，仍需要提供相应的本地归档目录。

路径相关设置如下：

- 环境变量：`BORYLXAT_RAW_CALC_ROOT`
- 默认回退路径：`E:/work/B_Cl_Nu`

例如在 PowerShell 中：

```powershell
$env:BORYLXAT_RAW_CALC_ROOT = "D:\\path\\to\\raw\\workspace"
```

## 推荐使用顺序

完整来源复现顺序：

```text
1_Calc_Reactant.ipynb
2_Calc_TS.ipynb
3_Build_DataBase.ipynb
4_Benchmark.ipynb
5_Modeling.ipynb
5_Modeling_xTB.ipynb
6_Draw_Figures.ipynb
```

如果只做审稿式复用，建议：

1. 在 `3_Build_DataBase.ipynb` 中仅运行数据库检查部分。
2. 在 `4_Benchmark.ipynb` 中从 `Benchmark.csv` 读取部分开始运行。
3. 在 `5_Modeling.ipynb` 或 `5_Modeling_xTB.ipynb` 中从描述符加载部分开始运行。
4. 需要修订分析时，直接运行对应的 `revision_*.ipynb`。
5. 使用 `6_Draw_Figures.ipynb` 统一导出论文和 SI 图片。

## 程序化读取示例

```python
from ase.db import connect

db = connect("BorylXAT-DB.db")
ts_rows = list(db.select(category="ts"))
print(len(ts_rows))
print(ts_rows[0].key_value_pairs)
```

## 补充说明

- `main.py` 只是占位脚本，真正的工作流在 notebook 和 `DFTStructureGenerator/` 模块中。
- 审稿运行路径与完整 Gaussian/xTB 来源路径是有意分开的。
- 原来的 `Data/csvs/Benchmark_Result.csv` 已退役，当前 benchmark 汇总表是 `Data/csvs/Benchmark.csv`。
- QHARM 修订分析依赖 `BorylXAT-DB_qh_update.db` 以及 `Data/descriptor/` 下的 QHARM 描述符文件。
