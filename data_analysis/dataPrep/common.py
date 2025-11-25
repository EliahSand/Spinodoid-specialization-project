from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence

import numpy as np
import pandas as pd

DEFAULT_DATA_ROOT = Path(__file__).resolve().parents[2] / "Matlab" / "results" / "sheets" / "lamellar"
DEFAULT_OUTPUT_ROOT = Path(__file__).resolve().parents[1] / "results"

DISPLACEMENTS = ["U1", "U2", "U3"]
STRESSES = ["S11", "S22", "S33", "S12", "S13", "S23"]
VON_MISES = ["SMises"]

_EPS = 1e-12


@dataclass
class CaseInfo:
    ratio: str
    angle: str
    case_name: str
    shell_path: Path
    solid_path: Optional[Path]

    @property
    def slug(self) -> str:
        if self.ratio and self.angle:
            return f"{self.ratio}_{self.angle}"
        if self.ratio:
            return self.ratio
        return self.case_name or "case"

    @property
    def human_label(self) -> str:
        if self.ratio and self.angle:
            return f"{self.ratio} / {self.angle}"
        return self.case_name


def extract_angle(case_name: str) -> str:
    match = re.search(r"ang(\d+)", case_name, flags=re.IGNORECASE)
    return f"ang{match.group(1)}" if match else ""


def find_solid_csv(case_dir: Path) -> Optional[Path]:
    solid_dir = case_dir / "FEA_solid"
    if not solid_dir.exists():
        return None
    candidates = list(solid_dir.glob("*.csv"))
    if not candidates:
        return None
    preferred = [p for p in candidates if "midplane" in p.name.lower()]
    preferred = preferred or [p for p in candidates if "solid" in p.name.lower()]
    preferred = preferred or candidates
    return preferred[0]


def discover_cases(data_root: Path = DEFAULT_DATA_ROOT) -> List[CaseInfo]:
    shell_files = data_root.glob("**/FEA_shell/midplane_results_shell.csv")
    cases: List[CaseInfo] = []
    for shell_path in shell_files:
        case_dir = shell_path.parent.parent
        ratio = case_dir.parent.name
        angle = extract_angle(case_dir.name)
        case_name = case_dir.name
        solid_path = find_solid_csv(case_dir)
        cases.append(
            CaseInfo(
                ratio=ratio,
                angle=angle,
                case_name=case_name,
                shell_path=shell_path,
                solid_path=solid_path,
            )
        )
    return sorted(cases, key=lambda c: (c.ratio, c.angle, c.case_name))


def _clean_columns(df: pd.DataFrame) -> pd.DataFrame:
    cleaned = df.copy()
    cleaned.columns = [col.strip() for col in cleaned.columns]
    return cleaned


def merge_case_data(case: CaseInfo, components: Sequence[str]) -> pd.DataFrame:
    if case.solid_path is None:
        raise FileNotFoundError(f"No SOLID CSV found for {case.case_name}")
    shell_df = _clean_columns(pd.read_csv(case.shell_path))
    solid_df = _clean_columns(pd.read_csv(case.solid_path))

    required_cols = ["Label", *components]
    missing_shell = set(required_cols) - set(shell_df.columns)
    missing_solid = set(required_cols) - set(solid_df.columns)
    if missing_shell:
        raise ValueError(f"SHELL data missing columns {sorted(missing_shell)} for {case.case_name}")
    if missing_solid:
        raise ValueError(f"SOLID data missing columns {sorted(missing_solid)} for {case.case_name}")

    shell_df = shell_df[required_cols]
    solid_df = solid_df[required_cols]

    merged = pd.merge(
        solid_df,
        shell_df,
        on="Label",
        suffixes=("_solid", "_shell"),
        how="inner",
    )
    return merged


def add_error_columns(df: pd.DataFrame, components: Sequence[str]) -> pd.DataFrame:
    enriched = df.copy()
    for comp in components:
        diff_col = f"{comp}_diff"
        rel_col = f"{comp}_rel_error_pct"
        solid_col = f"{comp}_solid"
        shell_col = f"{comp}_shell"

        enriched[diff_col] = enriched[solid_col] - enriched[shell_col]
        denom = np.where(np.abs(enriched[solid_col].to_numpy()) > _EPS, np.abs(enriched[solid_col].to_numpy()), _EPS)
        enriched[rel_col] = 100.0 * enriched[diff_col].to_numpy() / denom
    return enriched


def compute_error_summary(df: pd.DataFrame, components: Sequence[str]) -> pd.DataFrame:
    rows: List[Dict[str, float]] = []
    for comp in components:
        diff = df[f"{comp}_diff"]
        rel = df[f"{comp}_rel_error_pct"]
        solid = df[f"{comp}_solid"]
        shell = df[f"{comp}_shell"]
        corr = solid.corr(shell)
        rows.append(
            {
                "component": comp,
                "mean_abs_error": diff.abs().mean(),
                "rmse": float(np.sqrt(np.mean(np.square(diff)))),
                "max_abs_error": diff.abs().max(),
                "mean_abs_rel_error_pct": rel.abs().mean(),
                "pearson_r": corr,
            }
        )
    return pd.DataFrame(rows)


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def set_plot_style() -> None:
    import seaborn as sns

    sns.set_theme(context="talk", style="whitegrid", palette="deep")


def melt_relative_errors(df: pd.DataFrame, components: Sequence[str]) -> pd.DataFrame:
    records = []
    for comp in components:
        records.append(
            pd.DataFrame(
                {
                    "Label": df["Label"],
                    "component": comp,
                    "rel_error_pct": df[f"{comp}_rel_error_pct"],
                    "abs_error": df[f"{comp}_diff"].abs(),
                }
            )
        )
    return pd.concat(records, ignore_index=True)
