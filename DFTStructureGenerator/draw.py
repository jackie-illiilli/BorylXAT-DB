import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.interpolate import make_interp_spline
from scipy.stats import pearsonr
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score


# ===== Smooth histogram =====
def smooth_hist(data, bins=30):
    hist, edges = np.histogram(data, bins=bins)
    hist = hist / len(data)

    centers = 0.5 * (edges[:-1] + edges[1:])
    spline = make_interp_spline(centers, hist, k=3)

    x_smooth = np.linspace(centers.min(), centers.max(), 300)
    y_smooth = spline(x_smooth)

    return x_smooth, y_smooth


# ===== Plot distribution (independent function) =====
def plot_distribution(ax, data, orientation='top',
                      bins=30,
                      color_line='black',
                      color_fill='lightgray',
                      write_axis=False):

    xs, ys = smooth_hist(data, bins)

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
    figure_size=(5, 5),
    colors='coolwarm',
    annot=True,
    show_label=False,
    dpi=300,
):
    import seaborn as sns

    df = pd.DataFrame(x)
    correlation_matrix = np.abs(df.corr())
    mask = np.triu(np.ones_like(correlation_matrix, dtype=bool))
    print(np.max(correlation_matrix.to_numpy()[~mask]))
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
