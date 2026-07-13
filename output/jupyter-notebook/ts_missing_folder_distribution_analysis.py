from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re

import matplotlib.pyplot as plt
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
FIRST_12000_CSV = ROOT / "Data" / "csvs" / "First_12000.csv"
BORANE_ALL_CSV = ROOT / "Data" / "TS" / "Borane_all.csv"
OUTPUT_DIR = ROOT / "output" / "jupyter-notebook"
FIGURE_PATH = OUTPUT_DIR / "ts_missing_type_distribution.png"
SUMMARY_CSV_PATH = OUTPUT_DIR / "ts_missing_type_summary.csv"
COMBO_CSV_PATH = OUTPUT_DIR / "ts_missing_type_top_combinations.csv"
MISSING_CSV_PATH = OUTPUT_DIR / "ts_missing_rows_with_types.csv"

MATCH_KEYS = ["B_Index", "N_Index", "Cl_Index", "B_N_Cl_conf", "Cl_r_conf"]

LB_TYPE_ORDER = ["Amine/Aryl N", "Phosphine", "NHC"]
B_TYPE_ORDER = ["BR3", "R2BH", "RBH2", "BH3"]
CL_TYPE_ORDER = ["CCl4", "CCl3", "CCl2", "CCl"]


@dataclass
class SmilesGraph:
    atoms: list[str]
    neighbors: dict[int, list[int]]


def _parse_bracket_atom(token: str) -> str:
    match = re.match(r"^\d*([A-Z][a-z]?|[a-z])", token)
    if not match:
        raise ValueError(f"Unsupported bracket atom token: [{token}]")
    symbol = match.group(1)
    if len(symbol) == 1:
        return symbol.upper()
    return symbol[0].upper() + symbol[1:].lower()


def parse_smiles_graph(smiles: str) -> SmilesGraph:
    atoms: list[str] = []
    neighbors: dict[int, list[int]] = {}
    branch_stack: list[int | None] = []
    ring_starts: dict[str, int] = {}
    current: int | None = None
    i = 0

    def add_atom(symbol: str) -> int:
        idx = len(atoms)
        atoms.append(symbol)
        neighbors[idx] = []
        return idx

    while i < len(smiles):
        ch = smiles[i]

        if ch in "-=#:$/\\":
            i += 1
            continue
        if ch == ".":
            current = None
            i += 1
            continue
        if ch == "(":
            branch_stack.append(current)
            i += 1
            continue
        if ch == ")":
            current = branch_stack.pop()
            i += 1
            continue
        if ch == "[":
            j = smiles.index("]", i)
            symbol = _parse_bracket_atom(smiles[i + 1 : j])
            new_atom = add_atom(symbol)
            if current is not None:
                neighbors[current].append(new_atom)
                neighbors[new_atom].append(current)
            current = new_atom
            i = j + 1
            continue
        if ch == "%":
            ring_id = smiles[i + 1 : i + 3]
            if current is None:
                raise ValueError(f"Ring closure without current atom in {smiles}")
            if ring_id in ring_starts:
                start = ring_starts.pop(ring_id)
                neighbors[start].append(current)
                neighbors[current].append(start)
            else:
                ring_starts[ring_id] = current
            i += 3
            continue
        if ch.isdigit():
            ring_id = ch
            if current is None:
                raise ValueError(f"Ring closure without current atom in {smiles}")
            if ring_id in ring_starts:
                start = ring_starts.pop(ring_id)
                neighbors[start].append(current)
                neighbors[current].append(start)
            else:
                ring_starts[ring_id] = current
            i += 1
            continue

        if smiles.startswith("Cl", i):
            symbol = "Cl"
            i += 2
        elif smiles.startswith("Br", i):
            symbol = "Br"
            i += 2
        else:
            symbol = ch
            i += 1

        if symbol in {"c", "n", "o", "p", "s", "b"}:
            symbol = symbol.upper()

        if not re.fullmatch(r"[A-Z][a-z]?", symbol):
            raise ValueError(f"Unsupported atom token '{symbol}' in {smiles}")

        new_atom = add_atom(symbol)
        if current is not None:
            neighbors[current].append(new_atom)
            neighbors[new_atom].append(current)
        current = new_atom

    if ring_starts:
        raise ValueError(f"Unclosed ring in SMILES: {smiles}")

    return SmilesGraph(atoms=atoms, neighbors=neighbors)


def classify_lb_type(lb_index: int) -> str:
    lb_index = int(lb_index)
    if 0 <= lb_index < 69 or 234 <= lb_index < 388:
        return "Amine/Aryl N"
    if 69 <= lb_index < 144:
        return "Phosphine"
    if 144 <= lb_index < 234:
        return "NHC"
    return "Other"


def classify_b_type(smiles: str) -> str:
    graph = parse_smiles_graph(smiles)
    for atom_idx, symbol in enumerate(graph.atoms):
        if symbol == "B":
            heavy_neighbors = len(graph.neighbors[atom_idx])
            hydrogen_count = max(0, 3 - heavy_neighbors)
            return {0: "BR3", 1: "R2BH", 2: "RBH2", 3: "BH3"}.get(
                hydrogen_count,
                f"BH{hydrogen_count}",
            )
    raise ValueError(f"No boron atom found in {smiles}")


def classify_cl_type(smiles: str) -> str:
    graph = parse_smiles_graph(smiles)
    max_cl_neighbors = 0
    for atom_idx, symbol in enumerate(graph.atoms):
        cl_neighbors = sum(
            1 for neighbor_idx in graph.neighbors[atom_idx]
            if graph.atoms[neighbor_idx] == "Cl"
        )
        max_cl_neighbors = max(max_cl_neighbors, cl_neighbors)
    return {
        4: "CCl4",
        3: "CCl3",
        2: "CCl2",
    }.get(max_cl_neighbors, "CCl")


def add_type_columns(df: pd.DataFrame) -> pd.DataFrame:
    typed = df.copy()
    typed["B_type"] = typed["B_smiles"].map(classify_b_type)
    typed["LB_type"] = typed["N_Index"].map(classify_lb_type)
    typed["Cl_type"] = typed["Cl_smiles"].map(classify_cl_type)
    return typed


def build_missing_dataframe() -> pd.DataFrame:
    first_12000 = pd.read_csv(FIRST_12000_CSV)
    borane_all = pd.read_csv(BORANE_ALL_CSV)

    merged = first_12000.merge(
        borane_all[MATCH_KEYS],
        on=MATCH_KEYS,
        how="left",
        indicator=True,
    )
    missing = merged.loc[merged["_merge"] == "left_only"].drop(columns=["_merge"])
    return add_type_columns(missing)


def build_summary_tables(missing: pd.DataFrame, full: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    summary_frames: list[pd.DataFrame] = []
    for dimension, order in [
        ("B_type", B_TYPE_ORDER),
        ("Cl_type", CL_TYPE_ORDER),
        ("LB_type", LB_TYPE_ORDER),
    ]:
        missing_counts = missing[dimension].value_counts().reindex(order).fillna(0).astype(int)
        full_counts = full[dimension].value_counts().reindex(order).fillna(0).astype(int)
        summary = pd.DataFrame(
            {
                "dimension": dimension,
                "category": missing_counts.index,
                "missing_count": missing_counts.values,
                "missing_fraction": (missing_counts / missing_counts.sum()).values,
                "first_12000_count": full_counts.values,
                "missing_rate_within_first_12000": (
                    missing_counts / full_counts.replace(0, pd.NA)
                ).values,
            }
        )
        summary_frames.append(summary)

    summary_df = pd.concat(summary_frames, ignore_index=True)
    combo_df = (
        missing.groupby(["B_type", "Cl_type", "LB_type"])
        .size()
        .reset_index(name="missing_count")
        .sort_values("missing_count", ascending=False)
        .head(20)
        .reset_index(drop=True)
    )
    combo_df["missing_fraction"] = combo_df["missing_count"] / len(missing)
    return summary_df, combo_df


def plot_missing_type_distribution(missing: pd.DataFrame) -> Path:
    fig, axes = plt.subplots(1, 3, figsize=(15, 3.6), dpi=300)
    plot_specs = [
        ("B_type", B_TYPE_ORDER, ["#64cd9f", "#94e9c6", "#dbf8ec", "#eefcf7"], "Missing TS by borane type"),
        ("Cl_type", CL_TYPE_ORDER, ["#fcc556", "#fdcf73", "#fadda2", "#f6e4c0"], "Missing TS by chloride type"),
        ("LB_type", LB_TYPE_ORDER, ["#82baef", "#f68e56", "#f990e7"], "Missing TS by Lewis base type"),
    ]

    for ax, (column, order, colors, title) in zip(axes, plot_specs):
        counts = missing[column].value_counts().reindex(order).fillna(0).astype(int)
        labels = [
            f"{category}\n{count} ({count / len(missing) * 100:.1f}%)"
            for category, count in counts.items()
        ]
        ax.pie(
            counts.values,
            labels=labels,
            colors=colors[: len(counts)],
            startangle=0,
            wedgeprops={"edgecolor": "gray", "linewidth": 1},
            textprops={"fontsize": 9},
        )
        ax.set_title(title, fontsize=11)
        ax.axis("equal")

    plt.tight_layout()
    fig.savefig(FIGURE_PATH, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return FIGURE_PATH


def run_analysis() -> dict[str, object]:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    first_12000 = add_type_columns(pd.read_csv(FIRST_12000_CSV))
    missing = build_missing_dataframe()
    summary_df, combo_df = build_summary_tables(missing=missing, full=first_12000)
    figure_path = plot_missing_type_distribution(missing)

    missing.to_csv(MISSING_CSV_PATH, index=False)
    summary_df.to_csv(SUMMARY_CSV_PATH, index=False)
    combo_df.to_csv(COMBO_CSV_PATH, index=False)

    return {
        "missing": missing,
        "summary_df": summary_df,
        "combo_df": combo_df,
        "figure_path": figure_path,
        "missing_count": len(missing),
        "total_count": len(first_12000),
    }


def main() -> None:
    result = run_analysis()
    summary_df = result["summary_df"]
    combo_df = result["combo_df"]

    print(f"Missing TS: {result['missing_count']} / {result['total_count']}")
    print()
    for dimension in ["B_type", "Cl_type", "LB_type"]:
        print(f"[{dimension}]")
        print(
            summary_df.loc[summary_df["dimension"] == dimension, [
                "category",
                "missing_count",
                "missing_fraction",
                "first_12000_count",
                "missing_rate_within_first_12000",
            ]].to_string(index=False)
        )
        print()

    print("[Top missing combinations]")
    print(combo_df.to_string(index=False))
    print()
    print(f"Figure saved to: {result['figure_path']}")
    print(f"Summary CSV saved to: {SUMMARY_CSV_PATH}")
    print(f"Top combinations CSV saved to: {COMBO_CSV_PATH}")
    print(f"Detailed missing rows saved to: {MISSING_CSV_PATH}")


if __name__ == "__main__":
    main()
