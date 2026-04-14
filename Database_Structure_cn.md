# 硼自由基催化C-Cl活化反应数据库结构 (Database Structure)

基于 `Build_DataBase_v2.py` 脚本逻辑，该构建过程会生成两个核心数据库文件，用于结构化存储量子化学计算（DFT）得到的三维结构和热力学数据。两者各自用途不同：

1. **`boron_ccl.db`** (ASE SQLite 数据库): 主要配合 ASE (Atomic Simulation Environment) 使用，支持高效的结构对象存储、键值对查询、能垒和映射关系追踪。
2. **`boron_ccl_dataset.parquet`** (Parquet 数据集): 以二维表格格式存储完整的扁平化数据（包含坐标和原子序数），可以直接被 Pandas、PyArrow 加载，非常适合机器学习与深度学习任务（如搭建图神经网络）。

---

## 包含的化学物种体系 (Categories)

数据库中的所有化学结构通过命名规则 (`key`) 的正则匹配被归类为以下 7 种 (`category`):

| Category      | 描述                               | 开/闭壳层 | 特殊描述符    |
| ------------- | ---------------------------------- | --------- | ------------- |
| `B`         | 硼烷催化剂 (Borane)                | 开壳层    | 自旋密度      |
| `LB`        | 路易斯碱/亲核试剂 (Lewis Base)     | 闭壳层    | —            |
| `Cl`        | 反应底物 (氯代底物 `Cl_xxxxx_r`) | 闭壳层    | —            |
| `complex_r` | 反应物复合物 (B-LB reactant)       | 开壳层    | 自旋密度      |
| `complex_p` | 产物复合物 (B-LB product)          | 闭壳层    | —            |
| `ts`        | 过渡态 (Transition State)          | —        | 虚频/IRC/能垒 |
| `c_radical` | 碳自由基产物 (`Cl_xxxxx_p`)      | 开壳层    | 自旋密度      |

---

## 1. ASE SQLite 数据库 (`boron_ccl.db`)

ASE 数据库除了原生存储 `Atoms` 对象外，包含**键值对 (`key_value_pairs`)**（支持 SQL 式过滤）和**附加数据 (`data`)**（存储复杂的列表/字典）。

### 键值对 (`key_value_pairs` / `kvp`)

**所有物种共有字段：**

| 字段              | 类型  | 说明                       |
| ----------------- | ----- | -------------------------- |
| `category`      | str   | 物种类别（7种之一）        |
| `B_id`          | float | 硼烷编号                   |
| `LB_id`         | float | 路易斯碱编号               |
| `Cl_id`         | float | 氯代底物编号               |
| `gibbs_hartree` | float | 绝对吉布斯自由能 (Hartree) |
| `charge`        | int   | 体系净电荷                 |
| `temperature_K` | float | 热力学温度 (298.15)        |
| `solvent`       | str   | 隐式溶剂 ("toluene")       |
| `smiles`        | str   | SMILES / 反应AAM表达式     |
| `source_key`    | str   | 原始结构标识键             |

**非 TS 物种额外字段：**

| 字段                    | 类型  | 说明                     |
| ----------------------- | ----- | ------------------------ |
| `dipole_moment_debye` | float | 分子偶极矩 (Debye)       |
| `homo_energy_kcal`    | float | HOMO 轨道能量 (kcal/mol) |
| `lumo_energy_kcal`    | float | LUMO 轨道能量 (kcal/mol) |

**TS 过渡态专有字段：**

| 字段                         | 类型  | 说明                          |
| ---------------------------- | ----- | ----------------------------- |
| `is_transition_state`      | bool  | True                          |
| `expected_imaginary_freqs` | int   | 预期虚频个数 (1)              |
| `imaginary_frequency_cm_1` | float | 虚频频率 (cm⁻¹, 负值)       |
| `barrier_kcal`             | float | ΔG‡ 活化能 (kcal/mol)       |
| `delta_g_rxn_kcal`         | float | ΔG_rxn 反应自由能 (kcal/mol) |
| `reactant_complex_key`     | str   | 反应物复合物 key              |
| `reactant_cl_key`          | str   | 氯代底物 key                  |
| `product_complex_key`      | str   | 产物复合物 key                |
| `product_c_radical_key`    | str   | 碳自由基产物 key              |

### 附加数据块 (`data`)

**所有物种：**

| 字段        | 类型 | 说明       |
| ----------- | ---- | ---------- |
| `formula` | str  | 化学分子式 |

**除 B、LB 和 TS 之外的物种（Cl_r, c_radical, complex_r/p）：**

| 字段                  | 类型        | 说明                      |
| --------------------- | ----------- | ------------------------- |
| `hirshfeld_charges` | List[float] | 所有原子的 Hirshfeld 电荷 |

**B / complex_r / c_radical（开壳层）：**

| 字段               | 类型        | 说明                         |
| ------------------ | ----------- | ---------------------------- |
| `spin_densities` | List[float] | 所有原子的 Mulliken 自旋密度 |

**TS 过渡态：**

| 字段                            | 类型              | 说明                        |
| ------------------------------- | ----------------- | --------------------------- |
| `imaginary_freq_displacement` | List[List[float]] | 虚频振动位移矢量 (N×3)     |
| `irc_forward_positions`       | List[List[float]] | IRC 正向末端原子坐标 (N×3) |
| `irc_reverse_positions`       | List[List[float]] | IRC 逆向末端原子坐标 (N×3) |

**反应物 (complex_r & Cl_r)：**

| 字段                   | 类型      | 说明                      |
| ---------------------- | --------- | ------------------------- |
| `associated_ts_keys` | List[str] | 关联的所有过渡态 key 列表 |

---

## 2. Parquet 数据集 (`boron_ccl_dataset.parquet`)

统一的结构化二维表（`DataFrame`），包含以下列：

### 基础信息列

| 列名                           | 类型  | 说明                   |
| ------------------------------ | ----- | ---------------------- |
| `key`                        | str   | 结构唯一标识主键       |
| `category`                   | str   | 类别                   |
| `B_id`, `LB_id`, `Cl_id` | float | 组分编号               |
| `smiles`                     | str   | SMILES / AAM           |
| `natoms`                     | int   | 原子总数               |
| `formula`                    | str   | 化学分子式             |
| `charge`                     | int   | 系统总电荷             |
| `gibbs_hartree`              | float | 吉布斯自由能 (Hartree) |

### 3D空间特征列

| 列名          | 类型              | 说明                |
| ------------- | ----------------- | ------------------- |
| `numbers`   | List[int]         | 原子序数列表        |
| `positions` | List[List[float]] | 原子 3D 坐标 (N×3) |

### 电子结构描述符列

| 列名                    | 类型        | 适用类别                  | 说明                 |
| ----------------------- | ----------- | ------------------------- | -------------------- |
| `hirshfeld_charges`   | List[float] | 除 B, LB, TS 外的其他物种 | Hirshfeld 电荷       |
| `dipole_moment_debye` | float       | 非 TS 全部                | 偶极矩 (Debye)       |
| `spin_densities`      | List[float] | B / complex_r / c_radical | Mulliken 自旋密度    |
| `homo_energy_kcal`    | float       | 非 TS 全部                | HOMO 能量 (kcal/mol) |
| `lumo_energy_kcal`    | float       | 非 TS 全部                | LUMO 能量 (kcal/mol) |

### TS 过渡态专属列

| 列名                            | 类型              | 说明                          |
| ------------------------------- | ----------------- | ----------------------------- |
| `imaginary_frequency_cm_1`    | float             | 虚频频率 (cm⁻¹)             |
| `imaginary_freq_displacement` | List[List[float]] | 虚频位移矢量 (N×3)           |
| `irc_forward_positions`       | List[List[float]] | IRC 正向末端坐标              |
| `irc_reverse_positions`       | List[List[float]] | IRC 逆向末端坐标              |
| `barrier_kcal`                | float             | ΔG‡ 活化能 (kcal/mol)       |
| `delta_g_rxn_kcal`            | float             | ΔG_rxn 反应自由能 (kcal/mol) |
| `reactant_complex`            | str               | 反应物复合物 key              |
| `reactant_cl`                 | str               | 氯代底物 key                  |
| `product_complex`             | str               | 产物复合物 key                |
| `product_c_radical`           | str               | 碳自由基产物 key              |
