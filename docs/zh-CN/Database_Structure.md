# 硼自由基催化 C-Cl 活化反应数据库结构

[English Version](../../Database_Structure.md) | [项目 README](README.md)

基于 `Build_DataBase.py` 的实现逻辑，数据库构建流程会生成两个核心文件，用于结构化存储量子化学计算得到的三维结构和热力学数据。两种格式面向不同使用场景：

1. **`boron_ccl.db`**：ASE SQLite 数据库，适合通过 ASE 高效存储 `Atoms` 对象、检索键值对元数据并追踪反应映射关系。
2. **`boron_ccl_dataset.parquet`**：Parquet 数据集，以扁平化二维表的形式存储坐标、原子序数及描述符，更适合 Pandas、PyArrow 和机器学习任务。

---

## 包含的化学物种类别

数据库中的所有结构会根据 `key` 的命名规则被划分为以下七类：

| Category | 描述 | 开/闭壳层 | 特殊信息 |
| --- | --- | --- | --- |
| `B` | 硼自由基催化剂 | 开壳层 | 自旋密度 |
| `LB` | 路易斯碱 / 亲核试剂 | 闭壳层 | — |
| `Cl` | 氯代底物（`Cl_xxxxx_r`） | 闭壳层 | — |
| `complex_r` | 反应物复合物（B-LB reactant） | 开壳层 | 自旋密度 |
| `complex_p` | 产物复合物（B-LB product） | 闭壳层 | — |
| `ts` | 过渡态 | — | 虚频 / IRC / 能垒 |
| `c_radical` | 碳自由基产物（`Cl_xxxxx_p`） | 开壳层 | 自旋密度 |

---

## 1. ASE SQLite 数据库（`boron_ccl.db`）

除了 ASE 原生的 `Atoms` 对象外，数据库还包含：

- **键值对（`key_value_pairs`）**：便于筛选和查询
- **附加数据（`data`）**：用于存储列表、字典等复杂字段

### 键值对（`key_value_pairs` / `kvp`）

**所有物种共有字段**

| 字段 | 类型 | 说明 |
| --- | --- | --- |
| `category` | str | 结构类别 |
| `B_id` | float | 硼自由基编号 |
| `LB_id` | float | 路易斯碱编号 |
| `Cl_id` | float | 氯代底物编号 |
| `gibbs_hartree` | float | 绝对吉布斯自由能（Hartree） |
| `charge` | int | 体系总电荷 |
| `temperature_K` | float | 热力学温度（`298.15`） |
| `solvent` | str | 隐式溶剂（`toluene`） |
| `smiles` | str | SMILES 或反应 AAM 字符串 |
| `source_key` | str | 原始结构标识符 |

**非 TS 物种额外字段**

| 字段 | 类型 | 说明 |
| --- | --- | --- |
| `dipole_moment_debye` | float | 偶极矩（Debye） |
| `homo_energy_kcal` | float | HOMO 轨道能量（kcal/mol） |
| `lumo_energy_kcal` | float | LUMO 轨道能量（kcal/mol） |

**TS 专有字段**

| 字段 | 类型 | 说明 |
| --- | --- | --- |
| `is_transition_state` | bool | 对 TS 条目恒为 `True` |
| `expected_imaginary_freqs` | int | 预期虚频个数（`1`） |
| `imaginary_frequency_cm_1` | float | 虚频频率（cm^-1） |
| `barrier_kcal` | float | 活化自由能 `ΔG‡`（kcal/mol） |
| `delta_g_rxn_kcal` | float | 反应自由能 `ΔG_rxn`（kcal/mol） |
| `reactant_complex_key` | str | 反应物复合物的 key |
| `reactant_cl_key` | str | 氯代底物的 key |
| `product_complex_key` | str | 产物复合物的 key |
| `product_c_radical_key` | str | 碳自由基产物的 key |

### 附加数据块（`data`）

**所有物种共有字段**

| 字段 | 类型 | 说明 |
| --- | --- | --- |
| `formula` | str | 化学分子式 |

**除 B、LB 和 TS 外的物种（`Cl_r`、`c_radical`、`complex_r`、`complex_p`）**

| 字段 | 类型 | 说明 |
| --- | --- | --- |
| `hirshfeld_charges` | List[float] | 全部原子的 Hirshfeld 电荷 |

**开壳层物种：B / `complex_r` / `c_radical`**

| 字段 | 类型 | 说明 |
| --- | --- | --- |
| `spin_densities` | List[float] | 全部原子的 Mulliken 自旋密度 |

**TS 条目**

| 字段 | 类型 | 说明 |
| --- | --- | --- |
| `imaginary_freq_displacement` | List[List[float]] | 虚频模式位移向量（`N x 3`） |
| `irc_forward_positions` | List[List[float]] | IRC 正向末端坐标（`N x 3`） |
| `irc_reverse_positions` | List[List[float]] | IRC 逆向末端坐标（`N x 3`） |

**反应物（`complex_r` 与 `Cl_r`）**

| 字段 | 类型 | 说明 |
| --- | --- | --- |
| `associated_ts_keys` | List[str] | 相关过渡态 key 列表 |

---

## 2. Parquet 数据集（`boron_ccl_dataset.parquet`）

Parquet 文件保存为统一的二维表 `DataFrame`，主要包含以下几类列：

### 基础信息列

| 列名 | 类型 | 说明 |
| --- | --- | --- |
| `key` | str | 唯一结构标识符 |
| `category` | str | 结构类别 |
| `B_id`, `LB_id`, `Cl_id` | float | 组分编号 |
| `smiles` | str | SMILES 或反应 AAM 字符串 |
| `natoms` | int | 原子总数 |
| `formula` | str | 化学分子式 |
| `charge` | int | 体系总电荷 |
| `gibbs_hartree` | float | 吉布斯自由能（Hartree） |

### 三维结构列

| 列名 | 类型 | 说明 |
| --- | --- | --- |
| `numbers` | List[int] | 原子序数列表 |
| `positions` | List[List[float]] | 三维坐标（`N x 3`） |

### 电子结构描述符列

| 列名 | 类型 | 适用类别 | 说明 |
| --- | --- | --- | --- |
| `hirshfeld_charges` | List[float] | 除 B、LB、TS 外的物种 | Hirshfeld 电荷 |
| `dipole_moment_debye` | float | 所有非 TS 物种 | 偶极矩（Debye） |
| `spin_densities` | List[float] | B / `complex_r` / `c_radical` | Mulliken 自旋密度 |
| `homo_energy_kcal` | float | 所有非 TS 物种 | HOMO 能量（kcal/mol） |
| `lumo_energy_kcal` | float | 所有非 TS 物种 | LUMO 能量（kcal/mol） |

### TS 专属列

| 列名 | 类型 | 说明 |
| --- | --- | --- |
| `imaginary_frequency_cm_1` | float | 虚频频率（cm^-1） |
| `imaginary_freq_displacement` | List[List[float]] | 虚频模式位移向量（`N x 3`） |
| `irc_forward_positions` | List[List[float]] | IRC 正向末端坐标 |
| `irc_reverse_positions` | List[List[float]] | IRC 逆向末端坐标 |
| `barrier_kcal` | float | 活化自由能 `ΔG‡`（kcal/mol） |
| `delta_g_rxn_kcal` | float | 反应自由能 `ΔG_rxn`（kcal/mol） |
| `reactant_complex` | str | 反应物复合物 key |
| `reactant_cl` | str | 氯代底物 key |
| `product_complex` | str | 产物复合物 key |
| `product_c_radical` | str | 碳自由基产物 key |
