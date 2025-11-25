from __future__ import annotations

from pathlib import Path
from typing import Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from .common import ensure_dir


def scatter_shell_vs_solid(df: pd.DataFrame, components: Sequence[str], title: str, outfile: Path) -> None:
    n_components = len(components)
    fig, axes = plt.subplots(
        1,
        n_components,
        figsize=(5 * n_components, 5),
        sharex=False,
        sharey=False,
        constrained_layout=True,
    )
    axes = np.atleast_1d(axes)
    for ax, comp in zip(axes, components):
        x = df[f"{comp}_solid"]
        y = df[f"{comp}_shell"]
        sns.scatterplot(x=x, y=y, s=12, alpha=0.7, ax=ax, linewidth=0)
        combined = np.concatenate([x.to_numpy(), y.to_numpy()])
        min_val = np.nanmin(combined)
        max_val = np.nanmax(combined)
        span = max_val - min_val
        padding = 0.05 * span if span > 0 else 1.0
        ax.plot([min_val - padding, max_val + padding], [min_val - padding, max_val + padding], "k--", linewidth=1)
        r = x.corr(y)
        ax.text(
            0.02,
            0.98,
            f"r = {r:.3f}",
            transform=ax.transAxes,
            ha="left",
            va="top",
            bbox={"boxstyle": "round", "fc": "white", "ec": "gray", "alpha": 0.9},
        )
        ax.set_title(comp)
        ax.set_xlabel("SOLID")
        ax.set_ylabel("SHELL")
    fig.suptitle(title)
    ensure_dir(outfile.parent)
    fig.savefig(outfile, dpi=300)
    plt.close(fig)


def relative_error_histograms(df: pd.DataFrame, components: Sequence[str], title: str, outfile: Path) -> None:
    melted = pd.concat(
        [
            pd.DataFrame({"component": comp, "rel_error_pct": df[f"{comp}_rel_error_pct"]})
            for comp in components
        ],
        ignore_index=True,
    )
    g = sns.FacetGrid(
        melted,
        col="component",
        col_wrap=3,
        sharex=False,
        sharey=False,
        height=3.4,
        aspect=1.1,
    )
    g.map_dataframe(sns.histplot, x="rel_error_pct", kde=True, bins=40, color="C1")
    g.set_titles(col_template="{col_name}")
    g.set_axis_labels("Relative error (%)", "Count")
    plt.subplots_adjust(top=0.85)
    g.fig.suptitle(title)
    ensure_dir(outfile.parent)
    g.fig.savefig(outfile, dpi=300)
    plt.close(g.fig)


def absolute_error_boxplot(df: pd.DataFrame, components: Sequence[str], title: str, outfile: Path) -> None:
    melted = pd.concat(
        [
            pd.DataFrame({"component": comp, "abs_error": df[f"{comp}_diff"].abs()})
            for comp in components
        ],
        ignore_index=True,
    )
    fig, ax = plt.subplots(figsize=(max(6, 1.8 * len(components)), 4.5))
    sns.boxplot(data=melted, x="component", y="abs_error", ax=ax, color="C0", showfliers=False)
    sns.stripplot(
        data=melted,
        x="component",
        y="abs_error",
        ax=ax,
        color="gray",
        alpha=0.35,
        size=3,
        jitter=0.2,
    )
    ax.set_title(title)
    ax.set_xlabel("Component")
    ax.set_ylabel("Absolute error (SOLID - SHELL)")
    ensure_dir(outfile.parent)
    fig.savefig(outfile, dpi=300, bbox_inches="tight")
    plt.close(fig)


def summary_heatmap(
    summary_df: pd.DataFrame,
    value_col: str,
    title: str,
    outfile: Path,
    cmap: str = "vlag",
    center: float | None = 0.0,
) -> None:
    pivot = summary_df.pivot(index="case", columns="component", values=value_col)
    fig, ax = plt.subplots(figsize=(max(6, 1.2 * len(pivot.columns)), max(4, 0.5 * len(pivot.index))))
    sns.heatmap(pivot, annot=True, fmt=".3g", cmap=cmap, center=center, ax=ax, cbar_kws={"shrink": 0.8})
    ax.set_title(title)
    ax.set_xlabel("Component")
    ax.set_ylabel("Case")
    ensure_dir(outfile.parent)
    fig.savefig(outfile, dpi=300, bbox_inches="tight")
    plt.close(fig)


def mae_barplot(summary_df: pd.DataFrame, value_col: str, title: str, outfile: Path) -> None:
    fig, ax = plt.subplots(figsize=(9, 4.5))
    sns.barplot(data=summary_df, x="component", y=value_col, hue="case", ax=ax)
    ax.set_title(title)
    ax.set_xlabel("Component")
    ax.set_ylabel(value_col.replace("_", " ").title())
    ax.legend(title="Case", bbox_to_anchor=(1.02, 1), loc="upper left")
    ensure_dir(outfile.parent)
    fig.savefig(outfile, dpi=300, bbox_inches="tight")
    plt.close(fig)
