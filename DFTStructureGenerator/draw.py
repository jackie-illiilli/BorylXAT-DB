import os

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.interpolate import make_interp_spline
from scipy.stats import pearsonr
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score


# ===== Smooth histogram =====
def smooth_hist(data, bins=30, frequence=True):
    hist, edges = np.histogram(data, bins=bins)
    if frequence:
        hist = hist / len(data)

    centers = 0.5 * (edges[:-1] + edges[1:])
    spline = make_interp_spline(centers, hist, k=3)

    x_smooth = np.linspace(centers.min(), centers.max(), 300)
    y_smooth = spline(x_smooth)
    # Cubic spline smoothing can overshoot below zero for sparse histograms.
    y_smooth = np.clip(y_smooth, 0, None)

    return x_smooth, y_smooth


# ===== Plot distribution (independent function) =====
def plot_distribution(ax, data, orientation='top',
                      bins=30,
                      color_line='black',
                      color_fill='lightgray',
                      write_axis=False, 
                      frequence=True):

    xs, ys = smooth_hist(data, bins, frequence=frequence)

    if orientation == 'top':
        ax.plot(xs, ys, color=color_line, lw=1.2)
        ax.fill_between(xs, ys, 0, color=color_fill)

    elif orientation == 'right':
        ax.plot(ys, xs, color=color_line, lw=1.2)
        ax.fill_betweenx(xs, ys, 0, color=color_fill)

    if write_axis:
        ax.axis('on')
    else:
        ax.axis('off')


# ===== Scatter + linear fit =====
def plot_scatter_fit(ax, x, y,
                     color='black',
                     s=6,
                     linestyle='--'):

    # scatter
    ax.scatter(x, y, s=s, color=color, edgecolor='none', alpha=0.9)

    # linear fit
    coef = np.polyfit(x, y, 1)
    x_fit = np.linspace(x.min(), x.max(), 200)
    y_fit = np.polyval(coef, x_fit)

    ax.plot(x_fit, y_fit, color='black', lw=1.2, linestyle=linestyle)

    # statistics
    r, _ = pearsonr(x, y)
    y_pred = np.polyval(coef, x)
    r2 = r2_score(y, y_pred)

    print("coef:", coef)
    print(f"Pearson r: {r:.3f}")
    print(f"R^2: {r2:.3f}")

    return coef, r, r2


def plot_panel(ax_main, ax_top,
               x, y,
               xlabel,
               x_ticks,
               ylabel='',
               color_main="#215f9a",
               color_fill="#4e95d9"):

    # scatter + fit
    plot_scatter_fit(ax_main, x, y, color=color_main)

    # labels
    ax_main.set_xlabel(xlabel, fontsize=15)
    ax_main.set_ylabel(ylabel, fontsize=15)
    ax_main.set_xticks(x_ticks)
    ax_main.set_yticks(np.arange(0, 60, 10))
    ax_main.tick_params(labelsize=15)

    # distribution
    plot_distribution(ax_top, x,
                      orientation='top',
                      color_line=color_main,
                      color_fill=color_fill)


def observed_type_order(series, preferred_order):
    observed = list(series.dropna().unique())
    ordered = [each for each in preferred_order if each in observed]
    ordered.extend([each for each in sorted(observed) if each not in ordered])
    return ordered


def fit_line(x_values, y_values):
    if len(x_values) < 2 or len(np.unique(x_values)) < 2:
        return None, None, None
    coef = np.polyfit(x_values, y_values, 1)
    y_fit = np.polyval(coef, x_values)
    ss_res = np.sum((y_values - y_fit) ** 2)
    ss_tot = np.sum((y_values - np.mean(y_values)) ** 2)
    r2 = np.nan if ss_tot == 0 else 1 - ss_res / ss_tot
    return coef, y_fit, r2


def pearson_corr(x_values, y_values):
    if len(x_values) < 2 or np.std(x_values) == 0 or np.std(y_values) == 0:
        return np.nan
    return np.corrcoef(x_values, y_values)[0, 1]


def format_bep_annotation(coef, pearson_r, r2, residual_std):
    slope, intercept = coef
    sign = "+" if intercept >= 0 else "-"
    return (
        f"y = {slope:.2f}x {sign} {abs(intercept):.2f}\n"
        f"Pearson r = {pearson_r:.2f}, R$^2$ = {r2:.2f}\n"
        f"res. std. = {residual_std:.2f} kcal mol$^{{-1}}$"
    )


def top_bep_deviations(group_df, type_column, type_name, coef, x_col, y_col, top_n=3):
    deviation_df = group_df.copy()
    deviation_df["BEP_fit"] = np.polyval(coef, deviation_df[x_col].to_numpy())
    deviation_df["BEP_residual"] = deviation_df[y_col] - deviation_df["BEP_fit"]
    deviation_df["BEP_abs_residual"] = deviation_df["BEP_residual"].abs()
    deviation_df["type_axis"] = type_column
    deviation_df["type"] = type_name

    output_columns = [
        "type_axis",
        "type",
        "B_Index",
        "N_Index",
        "Cl_Index",
        "B_smiles",
        "N_smiles",
        "Cl_smiles",
        x_col,
        y_col,
        "BEP_fit",
        "BEP_residual",
        "BEP_abs_residual",
    ]
    return (
        deviation_df.sort_values("BEP_abs_residual", ascending=False)
        .head(top_n)[output_columns]
    )


def print_bep_deviation_table(top_df):
    type_axis = top_df["type_axis"].iloc[0]
    type_name = top_df["type"].iloc[0]
    print(f"  {type_axis} = {type_name}: top {len(top_df)} residual reactions")
    with pd.option_context("display.max_colwidth", None, "display.width", 240):
        print(
            top_df[
                [
                    "B_Index",
                    "N_Index",
                    "Cl_Index",
                    "B_smiles",
                    "N_smiles",
                    "Cl_smiles",
                    "BEP_residual",
                    "BEP_abs_residual",
                ]
            ].to_string(index=False)
        )
    print()


def safe_filename_token(value):
    token = "".join(ch if ch.isalnum() else "_" for ch in str(value))
    return token.strip("_") or "unknown"


def split_type_save_path(save_path, type_name):
    stem, ext = os.path.splitext(save_path)
    return f"{stem}_{safe_filename_token(type_name)}{ext}"


def plot_bep_single_type(
    group_df,
    type_column,
    type_name,
    color,
    save_path,
    title=None,
    top_n=3,
):
    x_col = "deltaG(kcal/mol)"
    y_col = "deltaGa(kcal/mol)"
    x_group = group_df[x_col].to_numpy()
    y_group = group_df[y_col].to_numpy()

    plt.rcParams["font.sans-serif"] = "Arial"
    fig = plt.figure(figsize=(4.5, 3.5), dpi=300)
    gs = fig.add_gridspec(5, 5)
    ax_main = fig.add_subplot(gs[1:4, 0:3])
    ax_top = fig.add_subplot(gs[0, 0:3], sharex=ax_main)
    ax_right = fig.add_subplot(gs[1:4, 3], sharey=ax_main)

    ax_main.scatter(
        x_group,
        y_group,
        s=8,
        color=color,
        edgecolor="none",
        alpha=0.75,
    )

    fit_row = None
    deviation_df = pd.DataFrame()
    coef, _, r2 = fit_line(x_group, y_group)
    if coef is not None:
        x_fit = np.linspace(x_group.min(), x_group.max(), 100)
        ax_main.plot(x_fit, np.polyval(coef, x_fit), color="black", lw=1.4)
        y_group_fit = np.polyval(coef, x_group)
        residual_std = (
            np.std(y_group - y_group_fit, ddof=1)
            if len(y_group) > 1
            else np.nan
        )
        pearson_r = pearson_corr(x_group, y_group)
        fit_row = {
            "type_axis": type_column,
            "type": type_name,
            "count": len(group_df),
            "slope": coef[0],
            "intercept": coef[1],
            "Pearson_r": pearson_r,
            "R2": r2,
            "residual_std": residual_std,
            "save_path": save_path,
        }
        annotation_text = format_bep_annotation(coef, pearson_r, r2, residual_std)
        ax_main.text(
            0.04,
            0.96,
            annotation_text,
            transform=ax_main.transAxes,
            ha="left",
            va="top",
            fontsize=10,
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.9, "pad": 3},
        )
        deviation_df = top_bep_deviations(
            group_df,
            type_column,
            type_name,
            coef,
            x_col,
            y_col,
            top_n=top_n,
        )
        print(f"Top {top_n} BEP residual reactions before plotting {title or type_column}: {type_name}")
        print_bep_deviation_table(deviation_df)
    else:
        print(f"Skip BEP fit for {type_column} = {type_name}: not enough x variation")

    plot_distribution(ax_top, x_group, "top", color_line=color, color_fill=color, frequence=True)
    plot_distribution(ax_right, y_group, "right", color_line=color, color_fill=color, frequence=True)

    ax_main.set_xlabel(r"$\Delta G_{rxn}$ (kcal/mol)", fontsize=15)
    ax_main.set_ylabel(r"$\Delta G^{\ddagger}$ (kcal/mol)", fontsize=15)
    ax_main.set_xticks(np.arange(-60, 15, 15))
    ax_main.set_yticks(np.arange(0, 60, 10))
    ax_main.tick_params(labelsize=15)

    plt.tight_layout()
    plt.subplots_adjust(hspace=0.02, wspace=0.02)
    save_dir = os.path.dirname(save_path)
    if save_dir:
        os.makedirs(save_dir, exist_ok=True)
    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.show()
    return fit_row, deviation_df


def plot_bep_by_type(
    df,
    type_column,
    type_order,
    colors,
    save_path,
    title=None,
    top_n=3,
):
    type_order = observed_type_order(df[type_column], type_order)
    color_map = {type_name: colors[i % len(colors)] for i, type_name in enumerate(type_order)}

    fit_rows = []
    deviation_tables = []
    for type_name in type_order:
        group_df = df.loc[df[type_column] == type_name]
        if group_df.empty:
            continue
        type_save_path = split_type_save_path(save_path, type_name)
        fit_row, deviation_df = plot_bep_single_type(
            group_df=group_df,
            type_column=type_column,
            type_name=type_name,
            color=color_map[type_name],
            save_path=type_save_path,
            title=title,
            top_n=top_n,
        )
        if fit_row is not None:
            fit_rows.append(fit_row)
        if not deviation_df.empty:
            deviation_tables.append(deviation_df)

    fit_df = pd.DataFrame(fit_rows)
    if deviation_tables:
        deviation_df = pd.concat(deviation_tables, ignore_index=True)
    else:
        deviation_df = pd.DataFrame()
    return fit_df, deviation_df


def calc_distribution_counts(
    y,
    eachsize=0.01,
    title=None,
    xlab=None,
    ylab="Count",
    y_max=None,
    y_min=None,
    figure_size=(4, 3),
    color="green",
    save_path=None,
    show=True,
):
    if y_max is None:
        y_max = np.max(y)
    if y_min is None:
        y_min = np.min(y)
    x_values = np.arange(y_min, y_max + eachsize, eachsize)
    counts = [0 for _ in x_values]
    z = (y - y_min) / eachsize
    for each in z:
        try:
            counts[int(each)] += 1
        except Exception:
            continue
    counts = np.array(counts)

    fig = plt.figure(figsize=figure_size)
    ax = fig.add_subplot(111)
    ax.patch.set_alpha(0.0)
    plt.bar(x_values, counts, width=eachsize / 2, color=color)
    plt.xlim(y_min - eachsize, y_max + eachsize)
    plt.ylim(0, np.max(counts) * 1.2)
    plt.xlabel(xlab, fontsize=30)
    plt.ylabel(ylab, fontsize=30)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    if title is not None:
        plt.title(title)
    plt.tight_layout()
    if save_path is not None:
        plt.savefig(save_path)
    if show:
        plt.show()
    return counts


def plot_scatter_with_metrics(
    x,
    y,
    title=None,
    min_=None,
    max_=None,
    figure_size=(5, 5),
    dpi=300,
    xlabel=r'Predict $\Delta G^{\ddagger}$',
    ylabel=r'Computed $\Delta G^{\ddagger}$',
    point_size=5,
    point_color="g",
    line_color="r",
    alpha=1,
    save_path=None,
    show=True,
):
    plt.figure(figsize=figure_size, dpi=dpi)
    plt.rcParams['font.sans-serif'] = 'Arial'
    if min_ is not None and max_ is not None:
        plt.xlim(min_, max_)
        plt.ylim(min_, max_)
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
    plt.scatter(x, y, s=point_size, c=point_color, alpha=alpha, edgecolors='none')
    plt.xlabel(xlabel, fontsize=25)
    plt.ylabel(ylabel, fontsize=25)
    if min_ is not None and max_ is not None:
        plt.plot([min_, max_], [min_, max_], c=line_color)
    if title is not None:
        plt.title(title, fontsize=24)
    if save_path is not None:
        plt.savefig(save_path)
    if show:
        plt.show()


def plot_line_with_metrics(
    x,
    y,
    title=None,
    min_=None,
    max_=None,
    figure_size=(7, 5),
    dpi=300,
    marker="*",
    color="g",
    save_path=None,
    show=True,
):
    r2 = r2_score(x, y)
    mae = mean_absolute_error(x, y)
    mse = mean_squared_error(x, y)

    plt.figure(figsize=figure_size, dpi=dpi)
    if min_ is not None and max_ is not None:
        plt.xlim(min_, max_)
        plt.ylim(min_, max_)
    plt.xticks(fontsize=24)
    plt.yticks(fontsize=24)
    if title is not None:
        plt.title(f"{title}\nR2:{r2:.3f}, MAE:{mae:.3f}, MSE:{mse:.3f}", fontsize=24)
    plt.plot(x, y, marker=marker, c=color)
    if save_path is not None:
        plt.savefig(save_path)
    if show:
        plt.show()


def draw_correlation_map(
    x,
    target=None,
    figure_size=(5, 5),
    colors='coolwarm',
    annot=True,
    show_label=False,
    dpi=300,
):
    import seaborn as sns

    df = pd.DataFrame(x)
    correlation_matrix = np.abs(df.corr())
    if target is not None:
        target_series = pd.Series(target, index=df.index)
        target_correlations = np.abs(df.corrwith(target_series))
        np.fill_diagonal(correlation_matrix.values, target_correlations.to_numpy())
    mask = np.triu(np.ones_like(correlation_matrix, dtype=bool))
    if target is not None:
        np.fill_diagonal(mask, False)
    print("Max correlation: ", np.max(correlation_matrix.to_numpy()[~mask]))
    fig, ax = plt.subplots(figsize=figure_size, dpi=dpi)
    plt.rcParams['font.sans-serif'] = 'Arial'
    annot_kws = {"fontsize": 10}
    ax = sns.heatmap(
        correlation_matrix,
        mask=mask,
        cmap=colors,
        annot=annot,
        fmt='.1f',
        center=0,
        cbar=1,
        annot_kws=annot_kws,
    )
    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize=25)
    cbar.ax.set_yticks(np.arange(0, 1, 0.2))
    if not show_label:
        ax.set_xticklabels('', fontsize=25)
        ax.set_yticklabels('', fontsize=25)
    plt.tight_layout()
    return correlation_matrix
