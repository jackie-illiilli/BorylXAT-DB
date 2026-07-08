from __future__ import annotations

from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd


matplotlib.use("Agg")

PROJECT_ROOT = Path(__file__).resolve().parents[2]
ATTEMPTED_CSV = PROJECT_ROOT / "Data" / "csvs" / "First_12000.csv"
FINAL_TS_CSV = PROJECT_ROOT / "Data" / "TS" / "Borane_all.csv"

OUTPUT_DIR = PROJECT_ROOT / "output" / "jupyter-notebook"
FIGURE_PATH = OUTPUT_DIR / "ts_missing_type_distribution.png"
COUNTS_PATH = OUTPUT_DIR / "ts_missing_type_counts.csv"
COMBO_PATH = OUTPUT_DIR / "ts_missing_type_combinations.csv"
ANNOTATED_PATH = OUTPUT_DIR / "ts_missing_typed_rows.csv"

B_TYPE_ORDER = ["BR3", "R2BH", "RBH2", "BH3"]
LB_TYPE_ORDER = ["Amine/Aryl N", "Phosphine", "NHC"]
CL_TYPE_ORDER = ["CCl4", "CCl3", "CCl2", "CCl"]

B_COLORS = ["#64cd9f", "#94e9c6", "#dbf8ec", "#eefcf6"]
CL_COLORS = ["#fcc556", "#fdcf73", "#fadda2", "#f6e4c0"]
LB_COLORS = ["#82baef", "#f68e56", "#f990e7"]

LB_TYPE_RANGES = [
    (range(0, 69), "Amine/Aryl N"),
    (range(69, 144), "Phosphine"),
    (range(144, 234), "NHC"),
    (range(234, 388), "Amine/Aryl N"),
]

SOURCE_TABLES = [
    PROJECT_ROOT / "Figure" / "Specificity" / "observed_specificity_ranked.csv",
    PROJECT_ROOT / "Figure" / "Orthogonal_Reactivity" / "observed_DFT_orthogonal_pairs.csv",
    PROJECT_ROOT / "output" / "jacs_experiment" / "jacs_experiment_external_validation_predictions.csv",
]

CL_TYPE_OVERRIDES = {
    480: "CCl",
}


def record_mapping(mapping: dict[int, str], key: int, value: str, label: str) -> None:
    key = int(key)
    value = str(value)
    previous = mapping.get(key)
    if previous is not None and previous != value:
        raise ValueError(f"Conflicting {label} mapping for {key}: {previous!r} vs {value!r}")
    mapping[key] = value


def build_type_maps() -> tuple[dict[int, str], dict[int, str]]:
    b_type_map: dict[int, str] = {}
    cl_type_map: dict[int, str] = {}

    for path in SOURCE_TABLES:
        df = pd.read_csv(path)

        for idx_col, type_col in [
            ("B_Index", "B_type"),
            ("BN1_B_Index", "BN1_B_type"),
            ("BN2_B_Index", "BN2_B_type"),
        ]:
            if idx_col in df.columns and type_col in df.columns:
                subset = df[[idx_col, type_col]].dropna().drop_duplicates()
                for _, row in subset.iterrows():
                    record_mapping(b_type_map, row[idx_col], row[type_col], "B_type")

        for idx_col, type_col in [
            ("Cl_Index", "Cl_type"),
            ("Cl1_Cl_Index", "Cl1_Cl_type"),
            ("Cl2_Cl_Index", "Cl2_Cl_type"),
        ]:
            if idx_col in df.columns and type_col in df.columns:
                subset = df[[idx_col, type_col]].dropna().drop_duplicates()
                for _, row in subset.iterrows():
                    record_mapping(cl_type_map, row[idx_col], row[type_col], "Cl_type")

    for key, value in CL_TYPE_OVERRIDES.items():
        record_mapping(cl_type_map, key, value, "Cl_type")

    return b_type_map, cl_type_map


def get_lb_type(lb_index: int) -> str:
    lb_index = int(lb_index)
    for lb_range, lb_type in LB_TYPE_RANGES:
        if lb_index in lb_range:
            return lb_type
    raise ValueError(f"Unmapped Lewis base index: {lb_index}")


def plot_pie(ax: plt.Axes, counts: pd.Series, colors: list[str], title: str) -> None:
    counts = counts[counts > 0]
    total = int(counts.sum())

    def autopct(pct: float) -> str:
        count = round(total * pct / 100.0)
        return f"{pct:.1f}%\n(n={count})"

    ax.pie(
        counts.values,
        labels=counts.index,
        colors=colors[: len(counts)],
        autopct=autopct,
        startangle=0,
        wedgeprops={"edgecolor": "gray", "linewidth": 1},
        textprops={"fontsize": 10},
    )
    ax.set_title(title, fontsize=12)
    ax.axis("equal")


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    attempted = pd.read_csv(ATTEMPTED_CSV)
    final_ts = pd.read_csv(FINAL_TS_CSV)
    b_type_map, cl_type_map = build_type_maps()

    id_cols = ["B_Index", "N_Index", "Cl_Index", "B_N_Cl_conf", "Cl_r_conf"]
    missing = attempted.merge(final_ts[id_cols].drop_duplicates(), on=id_cols, how="left", indicator=True)
    missing = missing.loc[missing["_merge"] == "left_only"].drop(columns="_merge").copy()

    missing["B_type"] = missing["B_Index"].map(b_type_map)
    missing["LB_type"] = missing["N_Index"].map(get_lb_type)
    missing["Cl_type"] = missing["Cl_Index"].map(cl_type_map)

    if missing["B_type"].isna().any():
        missing_ids = sorted(missing.loc[missing["B_type"].isna(), "B_Index"].unique())
        raise ValueError(f"Missing B_type mapping for B indices: {missing_ids}")
    if missing["Cl_type"].isna().any():
        missing_ids = sorted(missing.loc[missing["Cl_type"].isna(), "Cl_Index"].unique())
        raise ValueError(f"Missing Cl_type mapping for Cl indices: {missing_ids}")

    counts_records: list[dict[str, object]] = []
    count_specs = [
        ("B_type", B_TYPE_ORDER),
        ("LB_type", LB_TYPE_ORDER),
        ("Cl_type", CL_TYPE_ORDER),
    ]

    for column, order in count_specs:
        counts = missing[column].value_counts().reindex(order).fillna(0).astype(int)
        for type_name, count in counts.items():
            counts_records.append(
                {
                    "axis": column,
                    "type": type_name,
                    "count": int(count),
                    "fraction": float(count / len(missing)),
                }
            )

    counts_df = pd.DataFrame(counts_records)
    combo_df = (
        missing.groupby(["B_type", "LB_type", "Cl_type"], dropna=False)
        .size()
        .reset_index(name="count")
        .sort_values("count", ascending=False)
    )

    counts_df.to_csv(COUNTS_PATH, index=False)
    combo_df.to_csv(COMBO_PATH, index=False)
    missing.to_csv(ANNOTATED_PATH, index=False)

    fig, axes = plt.subplots(1, 3, figsize=(16, 4.2), dpi=300)

    plot_pie(
        axes[0],
        counts_df.loc[counts_df["axis"] == "B_type"].set_index("type").loc[B_TYPE_ORDER, "count"],
        B_COLORS,
        "Missing TS by borane type",
    )
    plot_pie(
        axes[1],
        counts_df.loc[counts_df["axis"] == "Cl_type"].set_index("type").loc[CL_TYPE_ORDER, "count"],
        CL_COLORS,
        "Missing TS by chloride type",
    )
    plot_pie(
        axes[2],
        counts_df.loc[counts_df["axis"] == "LB_type"].set_index("type").loc[LB_TYPE_ORDER, "count"],
        LB_COLORS,
        "Missing TS by Lewis base type",
    )

    fig.suptitle(f"Missing TS in First_12000 absent from Borane_all.csv (n={len(missing)})", fontsize=14)
    plt.tight_layout()
    plt.savefig(FIGURE_PATH, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"Missing TS rows: {len(missing)}")
    print()
    print(counts_df.to_string(index=False))
    print()
    print("Top 15 combinations:")
    print(combo_df.head(15).to_string(index=False))
    print()
    print(f"Saved figure: {FIGURE_PATH}")
    print(f"Saved counts: {COUNTS_PATH}")
    print(f"Saved combinations: {COMBO_PATH}")
    print(f"Saved annotated rows: {ANNOTATED_PATH}")


if __name__ == "__main__":
    main()
