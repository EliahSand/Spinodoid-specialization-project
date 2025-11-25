from __future__ import annotations

import argparse
from pathlib import Path
from typing import Sequence

import pandas as pd

from .common import (
    DEFAULT_DATA_ROOT,
    DEFAULT_OUTPUT_ROOT,
    add_error_columns,
    compute_error_summary,
    discover_cases,
    ensure_dir,
    merge_case_data,
    set_plot_style,
)
from .plotting import absolute_error_boxplot, mae_barplot, relative_error_histograms, scatter_shell_vs_solid, summary_heatmap


def run_component_group(group_name: str, components: Sequence[str]) -> None:
    parser = argparse.ArgumentParser(
        description=f"Compare SOLID vs SHELL results for {group_name}.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--data-root",
        type=Path,
        default=DEFAULT_DATA_ROOT,
        help="Root folder containing lamellar results (with trXX/.../FEA_shell).",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=DEFAULT_OUTPUT_ROOT,
        help="Where plots and metrics will be written.",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Process only the first N cases (useful for quick sanity checks).",
    )
    args = parser.parse_args()

    cases = discover_cases(args.data_root)
    set_plot_style()

    summaries = []
    for idx, case in enumerate(cases):
        if args.limit is not None and idx >= args.limit:
            break
        if case.solid_path is None:
            print(f"[skip] {case.case_name}: missing SOLID results; expected under {case.shell_path.parent.parent/'FEA_solid'}")
            continue
        try:
            merged = merge_case_data(case, components)
        except (FileNotFoundError, ValueError) as exc:
            print(f"[skip] {case.case_name}: {exc}")
            continue
        merged = add_error_columns(merged, components)

        plot_root = args.output_root / "plots" / group_name / case.slug
        title = f"{group_name} | {case.human_label}"
        scatter_shell_vs_solid(merged, components, title=title, outfile=plot_root / "shell_vs_solid_scatter.png")
        relative_error_histograms(merged, components, title=f"{title} | Relative errors", outfile=plot_root / "relative_error_hist.png")
        absolute_error_boxplot(merged, components, title=f"{title} | Absolute errors", outfile=plot_root / "absolute_error_boxplot.png")

        summary = compute_error_summary(merged, components)
        summary["case"] = case.slug
        summary["ratio"] = case.ratio
        summary["angle"] = case.angle
        summaries.append(summary)

    if not summaries:
        print("No cases with SOLID data available for comparison.")
        return

    summary_df = pd.concat(summaries, ignore_index=True)
    metrics_path = args.output_root / "metrics" / f"{group_name}_metrics.csv"
    ensure_dir(metrics_path.parent)
    summary_df.to_csv(metrics_path, index=False)
    print(f"[done] Wrote metrics to {metrics_path}")

    aggregate_root = args.output_root / "plots" / group_name / "aggregate"
    summary_heatmap(
        summary_df,
        value_col="pearson_r",
        title=f"{group_name} | Correlation (SOLID vs SHELL)",
        outfile=aggregate_root / "correlation_heatmap.png",
        cmap="coolwarm",
        center=0,
    )
    summary_heatmap(
        summary_df,
        value_col="mean_abs_error",
        title=f"{group_name} | Mean absolute error",
        outfile=aggregate_root / "mae_heatmap.png",
        cmap="crest",
        center=None,
    )
    mae_barplot(
        summary_df,
        value_col="mean_abs_error",
        title=f"{group_name} | Mean absolute error by case",
        outfile=aggregate_root / "mae_barplot.png",
    )

