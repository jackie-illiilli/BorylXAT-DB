import numpy as np
from scipy.interpolate import make_interp_spline
from scipy.stats import pearsonr
from sklearn.metrics import r2_score


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
    ax_main.set_yticks(np.arange(0,60,10))
    ax_main.tick_params(labelsize=15)

    # distribution
    plot_distribution(ax_top, x,
                      orientation='top',
                      color_line=color_main,
                      color_fill=color_fill)