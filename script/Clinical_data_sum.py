"""
Create one PDF table per cohort (SAMPLE SOURCE) from a prostate cancer metadata Excel file.

Outputs:
    ./Table_Cohort_<COHORT>.pdf      # e.g., Table_Cohort_BIDMC.pdf, Table_Cohort_PCBN.pdf, Table_Cohort_UM.pdf

Notes:
- Cohorts are taken from the `SAMPLE SOURCE` column.
- Patient-level summaries use a single record per PATIENT per cohort (first valid entry per field).
- Age and PSA are median (IQR). Categorical entries are % (N).
- Percentages use available data and may not sum to 100% due to rounding/missing values.

The column names expected are those in your Excel (e.g., SAMPLE SOURCE, PATIENT, AGE AT RP (YRS), PRE-OP PSA, PATIENT STAGE,
PATIENT ISUP GRADE, RACE, RECUR). Adjust `KEY_COLS` below if your columns differ.
"""

import os
import numpy as np
import pandas as pd

from reportlab.lib.pagesizes import letter
from reportlab.platypus import SimpleDocTemplate, Table, TableStyle, Paragraph, Spacer
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib.units import inch


# -----------------------------
# Configuration
# -----------------------------
INPUT_XLSX = "../data/metadata_for_eytan-BIDMC_PCBN_UM_tumor_benign.xlsx"  # <-- change path if needed
OUTPUT_DIR = "../results/clincal_data_sum"  # directory to save PDFs

# Columns we will use (present in your file); change if your schema changes.
KEY_COLS = [
    "SAMPLE ID",
    "PATIENT",
    "SAMPLE SOURCE",
    "SAMPLE TYPE",
    "PATIENT PRETREATED (NHT)",
    "PRESUMED SAMPLE PHENOTYPE",
    "PATIENT GLEASON SCORE",
    "PATIENT ISUP GRADE",
    "PATIENT STAGE",
    "RACE",
    "AGE AT RP (YRS)",
    "PRE-OP PSA",
    "RECUR",
    "TIME TO EVENT (BCR OR CENSORED, MONTHS)",
]

# Stage order for display
STAGE_ORDER = ["T2", "T2a", "T2b", "T2c", "T3", "T3a", "T3b", "T4", "Unk", "N/A"]


# -----------------------------
# Helpers
# -----------------------------
def detect_header_row(df_no_header: pd.DataFrame, search="SAMPLE ID", max_rows=20) -> int:
    """Detect the row index that contains the header (by finding a cell equal to `search`)."""
    for i in range(min(max_rows, len(df_no_header))):
        if any(str(v).strip().upper() == search.upper() for v in df_no_header.iloc[i].tolist()):
            return i
    return 0  # fallback


def first_valid(series):
    """Return the first non-NA value in a series (used for patient-level aggregation)."""
    for x in series:
        if pd.notna(x):
            return x
    return np.nan


def median_iqr(series) -> str:
    """Median (IQR) string for a numeric series."""
    s = pd.to_numeric(series, errors="coerce").dropna()
    if len(s) == 0:
        return "NA"
    med = np.median(s)
    q1 = np.percentile(s, 25)
    q3 = np.percentile(s, 75)
    return f"{med:.1f} ({q1:.1f}–{q3:.1f})"


def build_cohort_rows(sub: pd.DataFrame) -> list:
    """Build the 2-column rows for a single-cohort table."""
    rows = []

    # 1) High-level counts/medians
    n_pat = sub["PATIENT"].nunique()
    age_str = median_iqr(sub["AGE AT RP (YRS)"]) if "AGE AT RP (YRS)" in sub.columns else "NA"
    psa_str = median_iqr(sub["PRE-OP PSA"]) if "PRE-OP PSA" in sub.columns else "NA"

    rows.append(["Patients (N)", f"{n_pat}"])
    rows.append(["Age at RP (years), median (IQR)", age_str])
    rows.append(["Pre-op PSA (ng/mL), median (IQR)", psa_str])
    rows.append(["", ""])

    # 2) ISUP grade distribution
    rows.append(["ISUP grade, % (N)", ""])
    if "PATIENT ISUP GRADE" in sub.columns:
        isup_counts = (
            sub["PATIENT ISUP GRADE"]
            .round(0)
            .dropna()
            .astype(int)
            .value_counts()
            .sort_index()
        )
        for g in [1, 2, 3, 4, 5]:
            cnt = int(isup_counts.get(g, 0))
            pct = (cnt / n_pat * 100) if n_pat else 0
            rows.append([f"  Grade {g}", f"{pct:.1f}% ({cnt})"])
    else:
        rows.append(["  (no ISUP data)", "NA"])
    rows.append(["", ""])

    # 3) Clinical stage distribution
    rows.append(["Clinical stage, % (N)", ""])
    if "PATIENT STAGE" in sub.columns:
        stage_counts = (
            sub["PATIENT STAGE"].astype(str).str.strip().replace({"nan": np.nan}).value_counts(dropna=False)
        )
        ordered = [st for st in STAGE_ORDER if st in stage_counts.index]
        # add any extra stages we didn't anticipate
        for st in stage_counts.index:
            if st not in ordered and st not in ["nan", "NaN", "None"]:
                ordered.append(st)
        for st in ordered:
            cnt = int(stage_counts.get(st, 0))
            pct = (cnt / n_pat * 100) if n_pat else 0
            rows.append([f"  {st}", f"{pct:.1f}% ({cnt})"])
    else:
        rows.append(["  (no stage data)", "NA"])
    rows.append(["", ""])

    # 4) Race distribution (top categories + Other)
    rows.append(["Race, % (N)", ""])
    if "RACE" in sub.columns:
        race_counts = sub["RACE"].astype(str).str.upper().replace({"NAN": "UNKNOWN"}).value_counts()
        total = race_counts.sum()
        top5 = race_counts.head(5)
        for race, cnt in top5.items():
            pct = (cnt / total * 100) if total else 0
            rows.append([f"  {race.title()}", f"{pct:.1f}% ({int(cnt)})"])
        other = race_counts.iloc[5:].sum()
        if other > 0:
            pct = (other / total * 100) if total else 0
            rows.append(["  Other", f"{pct:.1f}% ({int(other)})"])
    else:
        rows.append(["  (no race data)", "NA"])
    rows.append(["", ""])

    # 5) Biochemical recurrence
    rows.append(["Biochemical recurrence, % (N)", ""])
    if "RECUR" in sub.columns:
        recur_counts = sub["RECUR"].round(0).dropna().astype(int).value_counts()
        yes_cnt = int(recur_counts.get(1, 0))
        no_cnt = int(recur_counts.get(0, 0))
        total = yes_cnt + no_cnt
        yes_pct = (yes_cnt / total * 100) if total else 0
        no_pct = (no_cnt / total * 100) if total else 0
        rows.append(["  Yes", f"{yes_pct:.1f}% ({yes_cnt})"])
        rows.append(["  No", f"{no_pct:.1f}% ({no_cnt})"])
    else:
        rows.append(["  (no recurrence data)", "NA"])

    return rows


def save_pdf_for_cohort(cohort_name: str, rows: list, out_dir: str = "."):
    """Render a single-cohort two-column PDF table."""
    os.makedirs(out_dir, exist_ok=True)
    pdf_path = os.path.join(out_dir, f"Table_Cohort_{cohort_name}.pdf")

    doc = SimpleDocTemplate(
        pdf_path, pagesize=letter, rightMargin=36, leftMargin=36, topMargin=36, bottomMargin=36
    )
    styles = getSampleStyleSheet()
    story = []

    title = Paragraph(f"Supplementary Table — Patient Characteristics ({cohort_name})", styles["Title"])
#    subtitle = Paragraph(
#        "Values computed at the <i>patient</i> level within this cohort only. "
#        "Age and PSA are reported as median (IQR).",
#        styles["Normal"],
#    )

    story.append(title)
    story.append(Spacer(1, 0.2 * inch))
#    story.append(subtitle)
#    story.append(Spacer(1, 0.2 * inch))

    tbl = Table(rows, colWidths=[3.4 * inch, 3.0 * inch], repeatRows=0)
    ts = TableStyle(
        [
            ("FONTNAME", (0, 0), (-1, -1), "Helvetica"),
            ("FONTSIZE", (0, 0), (-1, -1), 9),
            ("ALIGN", (0, 0), (0, -1), "LEFT"),
            ("ALIGN", (1, 0), (1, -1), "CENTER"),
            ("GRID", (0, 0), (-1, -1), 0.25, colors.HexColor("#CCCCCC")),
            ("LEFTPADDING", (0, 0), (-1, -1), 4),
            ("RIGHTPADDING", (0, 0), (-1, -1), 4),
            ("TOPPADDING", (0, 0), (-1, -1), 2),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 2),
        ]
    )

    # Emphasize section headers
    section_headers = {"ISUP grade, % (N)", "Clinical stage, % (N)", "Race, % (N)", "Biochemical recurrence, % (N)"}
    for i, row in enumerate(rows):
        if row[0] in section_headers:
            ts.add("FONTNAME", (0, i), (-1, i), "Helvetica-Bold")
            ts.add("BACKGROUND", (0, i), (-1, i), colors.HexColor("#F5F7FA"))

    tbl.setStyle(ts)
    story.append(tbl)

#    foot = Paragraph(
#        'Notes: Cohort derived from "SAMPLE SOURCE". One record per PATIENT in this cohort. '
#        "Percentages use available data and may not sum to 100% due to rounding or missing values.",
#        styles["Italic"],
#    )
#    story.append(Spacer(1, 0.15 * inch))
#    story.append(foot)

    doc.build(story)
    print(f"[Saved] {pdf_path}")


def main():
    # 1) Load Excel; detect header row (some files place headers below a title row)
    raw = pd.read_excel(INPUT_XLSX, header=None, engine="openpyxl")
    header_row_idx = detect_header_row(raw, search="SAMPLE ID", max_rows=20)

    header = raw.iloc[header_row_idx].astype(str).str.strip().tolist()
    df = raw.iloc[header_row_idx + 1 :].copy()
    df.columns = header

    # Keep only the columns we need (if present)
    keep = [c for c in KEY_COLS if c in df.columns]
    df = df[keep].copy()

    if "SAMPLE SOURCE" not in df.columns:
        raise ValueError('Column "SAMPLE SOURCE" is required but not found.')

    # 2) Collapse to patient-level per cohort
    agg_map = {}
    for col in [
        "AGE AT RP (YRS)",
        "PRE-OP PSA",
        "PATIENT STAGE",
        "PATIENT ISUP GRADE",
        "PATIENT GLEASON SCORE",
        "RACE",
        "RECUR",
        "TIME TO EVENT (BCR OR CENSORED, MONTHS)",
    ]:
        if col in df.columns:
            agg_map[col] = first_valid

    patient_df = (
        df.groupby(["SAMPLE SOURCE", "PATIENT"], dropna=False)
        .agg(agg_map)
        .reset_index()
    )

    # Normalize numeric columns
    for col in ["AGE AT RP (YRS)", "PRE-OP PSA", "PATIENT ISUP GRADE", "RECUR", "TIME TO EVENT (BCR OR CENSORED, MONTHS)"]:
        if col in patient_df.columns:
            patient_df[col] = pd.to_numeric(patient_df[col], errors="coerce")

    # Clean up stage string spacing
    if "PATIENT STAGE" in patient_df.columns:
        patient_df["PATIENT STAGE"] = patient_df["PATIENT STAGE"].astype(str).str.strip()

    # 3) Generate one PDF per cohort present in the data
    cohorts = (
        patient_df["SAMPLE SOURCE"]
        .dropna()
        .astype(str)
        .sort_values()
        .unique()
        .tolist()
    )

    for cohort in cohorts:
        sub = patient_df[patient_df["SAMPLE SOURCE"] == cohort].copy()
        rows = build_cohort_rows(sub)
        save_pdf_for_cohort(cohort, rows, out_dir=OUTPUT_DIR)


if __name__ == "__main__":
    main()