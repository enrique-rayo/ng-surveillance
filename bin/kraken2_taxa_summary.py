#!/usr/bin/env python3
"""Create batch-level Kraken2 taxa and descriptive QC summaries."""

import argparse
import base64
import csv
import html
import io
from pathlib import Path

from qc_flagging import MAD_MULTIPLIER, MAD_SCALE, robust_limits

RANKS = {"F": "family", "G": "genus", "S": "species"}


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--results-dir",
        required=True,
        type=Path,
        help="Directory recursively containing *_kraken_report.txt and *_kraken_summary.csv",
    )
    parser.add_argument("--batch-name", required=True)
    parser.add_argument("--csv-output", required=True, type=Path)
    parser.add_argument("--html-output", required=True, type=Path)
    return parser.parse_args()


def sample_from_path(path, suffix):
    return path.name[: -len(suffix)]


def parse_report(path):
    taxa = {rank: [] for rank in RANKS}
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 6:
                raise ValueError(
                    f"{path}:{line_number}: expected 6 tab-separated fields, got {len(fields)}"
                )
            rank = fields[3].strip()
            if rank in RANKS:
                taxa[rank].append(
                    {
                        "taxon_name": fields[5].strip(),
                        "taxid": fields[4].strip(),
                        "pct_of_reads": float(fields[0].strip().replace(",", ".")),
                        "reads_assigned_direct": int(fields[2].strip()),
                    }
                )
    for rank in taxa:
        taxa[rank].sort(
            key=lambda row: (-row["pct_of_reads"], -row["reads_assigned_direct"], row["taxid"])
        )
        taxa[rank] = taxa[rank][:5]
    return taxa


def read_qc_summaries(paths):
    qc = {}
    required = {
        "sample",
        "unclassified_pct",
        "species_pct",
        "genus_only_unresolved_pct",
    }
    for path in paths:
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            if not required.issubset(reader.fieldnames or []):
                missing = ", ".join(sorted(required - set(reader.fieldnames or [])))
                raise ValueError(f"{path}: missing required column(s): {missing}")
            rows = list(reader)
        if len(rows) != 1:
            raise ValueError(f"{path}: expected exactly one data row, got {len(rows)}")
        row = rows[0]
        sample = row["sample"].strip()
        qc[sample] = {
            "ng": float(row["species_pct"].replace(",", ".")),
            "unclassified": float(row["unclassified_pct"].replace(",", ".")),
            "genus_only": float(row["genus_only_unresolved_pct"].replace(",", ".")),
        }
    return qc


def classify(qc, baseline_samples=None):
    baseline_samples = set(qc) if baseline_samples is None else set(baseline_samples)
    if not baseline_samples:
        raise ValueError("Kraken2 baseline must contain at least one sample")
    unknown = baseline_samples - set(qc)
    if unknown:
        raise ValueError(f"Kraken2 baseline contains unknown samples: {sorted(unknown)}")
    baseline = [qc[sample] for sample in baseline_samples]
    ng_limits = robust_limits([row["ng"] for row in baseline])
    uc_limits = robust_limits([row["unclassified"] for row in baseline])
    genus_limits = robust_limits([row["genus_only"] for row in baseline])
    for row in qc.values():
        row["scatter_category"] = (
            "FLAG"
            if (
                abs(row["ng"] - ng_limits["median"])
                > MAD_MULTIPLIER * ng_limits["scaled_mad"]
                or abs(row["unclassified"] - uc_limits["median"])
                > MAD_MULTIPLIER * uc_limits["scaled_mad"]
            )
            else "TYPICAL"
        )
        row["genus_category"] = (
            "FLAG"
            if abs(row["genus_only"] - genus_limits["median"])
            > MAD_MULTIPLIER * genus_limits["scaled_mad"]
            else "TYPICAL"
        )
    return ng_limits, uc_limits, genus_limits


def make_plot(qc, batch_name, ng_limits, uc_limits):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(12, 8))
    colors = {"TYPICAL": "#2e8b57", "FLAG": "#d62728"}
    for category in ("TYPICAL", "FLAG"):
        rows = [
            (sample, row)
            for sample, row in qc.items()
            if row["scatter_category"] == category
        ]
        if rows:
            ax.scatter(
                [row["ng"] for _, row in rows],
                [row["unclassified"] for _, row in rows],
                color=colors[category],
                label=category,
                s=55,
                zorder=3,
            )
            for sample, row in rows:
                ax.annotate(
                    sample,
                    (row["ng"], row["unclassified"]),
                    xytext=(4, 4),
                    textcoords="offset points",
                    fontsize=7,
                )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("% N. gonorrhoeae (taxid 485)")
    ax.set_ylabel("% Unclassified (taxid 0)")
    ax.set_title(f"N. gonorrhoeae vs. Unclassified (batch: {batch_name})", pad=16)
    ax.axvline(ng_limits["median"], color="#444444", linestyle=":", linewidth=1.5)
    ax.axhline(uc_limits["median"], color="#444444", linestyle=":", linewidth=1.5)
    for value in (ng_limits["lower"], ng_limits["upper"]):
        if value > 0:
            ax.axvline(value, color="#777777", linestyle=":", linewidth=1, alpha=0.45)
    for value in (uc_limits["lower"], uc_limits["upper"]):
        if value > 0:
            ax.axhline(value, color="#777777", linestyle=":", linewidth=1, alpha=0.45)
    ax.grid(True, which="both", alpha=0.15)
    ax.legend()
    caption = (
        "Flagged if >3 scaled MAD from batch median "
        f"(median Ng%: {ng_limits['median']:.3f}, scaled MAD: {ng_limits['scaled_mad']:.3f}; "
        f"median Unclassified%: {uc_limits['median']:.3f}, "
        f"scaled MAD: {uc_limits['scaled_mad']:.3f})"
    )
    fig.text(0.5, 0.015, caption, ha="center", fontsize=9)
    fig.tight_layout(rect=(0, 0.045, 1, 1))
    image = io.BytesIO()
    fig.savefig(image, format="png", dpi=160)
    plt.close(fig)
    return base64.b64encode(image.getvalue()).decode("ascii")


def make_genus_only_plot(qc, batch_name, limits):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rows = sorted(qc.items(), key=lambda item: item[1]["genus_only"])
    fig_height = max(8, len(rows) * 0.28)
    fig, ax = plt.subplots(figsize=(12, fig_height))
    colors = {
        "TYPICAL": "#2e8b57",
        "FLAG": "#d62728",
    }
    positions = list(range(len(rows)))
    ax.barh(
        positions,
        [row["genus_only"] for _, row in rows],
        color=[colors[row["genus_category"]] for _, row in rows],
    )
    ax.set_yticks(positions)
    ax.set_yticklabels([sample for sample, _ in rows], fontsize=8)
    ax.set_xlabel("% reads assigned directly to genus Neisseria (taxid 482)")
    ax.set_title(f"Genus-only unresolved Neisseria reads (batch: {batch_name})", pad=14)
    ax.axvline(limits["median"], color="#444444", linestyle=":", linewidth=1.5)
    for value in (limits["lower"], limits["upper"]):
        if value >= 0:
            ax.axvline(value, color="#777777", linestyle=":", linewidth=1, alpha=0.45)
    for position, (_, row) in zip(positions, rows):
        ax.text(
            row["genus_only"],
            position,
            f" {row['genus_only']:.3f}%",
            va="center",
            fontsize=7,
        )
    ax.grid(True, axis="x", alpha=0.15)
    caption = (
        "Flagged if >3 scaled MAD from batch median "
        f"(median: {limits['median']:.4f}%, raw MAD: {limits['raw_mad']:.4f}, "
        f"scaled MAD: {limits['scaled_mad']:.4f}; "
        f"3-MAD interval: [{limits['lower']:.4f}, {limits['upper']:.4f}]%)"
    )
    fig.text(0.5, 0.01, caption, ha="center", fontsize=9)
    fig.tight_layout(rect=(0, 0.035, 1, 1))
    image = io.BytesIO()
    fig.savefig(image, format="png", dpi=160)
    plt.close(fig)
    return base64.b64encode(image.getvalue()).decode("ascii")


def write_csv(path, taxa_by_sample):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(
            [
                "sample",
                "rank_level",
                "taxon_name",
                "taxid",
                "pct_of_reads",
                "reads_assigned_direct",
            ]
        )
        for sample in sorted(taxa_by_sample):
            for rank in ("F", "G", "S"):
                for row in taxa_by_sample[sample][rank]:
                    writer.writerow(
                        [
                            sample,
                            RANKS[rank],
                            row["taxon_name"],
                            row["taxid"],
                            f"{row['pct_of_reads']:.2f}",
                            row["reads_assigned_direct"],
                        ]
                    )


def fmt_limits(label, limits):
    return (
        f"<strong>{html.escape(label)}</strong>: median {limits['median']:.3f}%; "
        f"raw MAD {limits['raw_mad']:.3f}; scaled MAD {limits['scaled_mad']:.3f}; "
        f"3-MAD interval [{limits['lower']:.3f}, {limits['upper']:.3f}]%"
    )


def write_html(
    path,
    batch_name,
    scatter_plot,
    genus_plot,
    taxa_by_sample,
    qc,
    ng_limits,
    uc_limits,
    genus_limits,
):
    out = [
        "<!doctype html><html><head><meta charset='utf-8'>",
        f"<title>Kraken2 batch summary — {html.escape(batch_name)}</title>",
        "<style>body{font-family:Arial,sans-serif;max-width:1400px;margin:2rem auto;"
        "padding:0 1rem;color:#222}img{max-width:100%;height:auto}table{border-collapse:"
        "collapse;margin-bottom:1.5rem}th,td{border:1px solid #bbb;padding:.35rem .55rem;"
        "text-align:right}th:first-child,td:first-child{text-align:left}.flag{color:#b00020;"
        "font-weight:bold}.typical{color:#16733a;font-weight:bold}.note{background:#fff4d6;"
        "padding:1rem;border-left:5px solid #d99b00}.plots{display:grid;grid-template-"
        "columns:repeat(auto-fit,minmax(min(100%,600px),1fr));gap:1rem;align-items:start}"
        ".plot-note{background:#eef5ff;padding:.8rem}</style></head><body>",
        f"<h1>Kraken2 batch summary: {html.escape(batch_name)}</h1>",
        "<p class='note'><strong>Interpretation:</strong> These are descriptive, batch-relative "
        "flags for review, not validated diagnostic thresholds. A flagged sample needs human "
        "review, not automatic exclusion.</p>",
        "<div class='plots'><div>",
        f"<img alt='Kraken2 QC scatter plot' src='data:image/png;base64,{scatter_plot}'>",
        "</div><div>",
        "<p class='plot-note'><strong>Genus-only unresolved:</strong> reads confidently placed "
        "in genus <em>Neisseria</em> but not resolved to any species. This differs from "
        "unclassified reads, which have no confident placement, and from "
        "<em>N. gonorrhoeae</em> reads resolved at species level.</p>",
        f"<img alt='Genus-only unresolved Neisseria bar chart' src='data:image/png;base64,{genus_plot}'>",
        "</div></div>",
        "<h2>QC statistics and categories</h2>",
        "<p>Flag rule: absolute distance from the batch median greater than three scaled MAD, "
        f"where scaled MAD = {MAD_SCALE} × raw MAD.</p>",
        f"<p>{fmt_limits('% N. gonorrhoeae', ng_limits)}<br>"
        f"{fmt_limits('% Unclassified', uc_limits)}<br>"
        f"{fmt_limits('% Genus-only unresolved Neisseria', genus_limits)}</p>",
        "<table><thead><tr><th>Sample</th><th>Ng %</th><th>Unclassified %</th>"
        "<th>Scatter category</th><th>Genus-only %</th>"
        "<th>Genus-only category</th></tr></thead><tbody>",
    ]
    for sample in sorted(qc):
        row = qc[sample]
        scatter_css = row["scatter_category"].lower()
        genus_css = row["genus_category"].lower()
        out.append(
            f"<tr><td>{html.escape(sample)}</td><td>{row['ng']:.2f}</td>"
            f"<td>{row['unclassified']:.2f}</td>"
            f"<td class='{scatter_css}'>{row['scatter_category']}</td>"
            f"<td>{row['genus_only']:.6f}</td>"
            f"<td class='{genus_css}'>{row['genus_category']}</td></tr>"
        )
    out.append("</tbody></table><h2>Top taxa by sample</h2>")
    for sample in sorted(taxa_by_sample):
        out.append(f"<section><h3>{html.escape(sample)}</h3>")
        for rank in ("F", "G", "S"):
            out.append(
                f"<h4>{RANKS[rank].title()}</h4><table><thead><tr><th>Taxon</th>"
                "<th>Taxid</th><th>% reads</th><th>Reads assigned directly</th>"
                "</tr></thead><tbody>"
            )
            for row in taxa_by_sample[sample][rank]:
                out.append(
                    f"<tr><td>{html.escape(row['taxon_name'])}</td>"
                    f"<td>{html.escape(row['taxid'])}</td><td>{row['pct_of_reads']:.2f}</td>"
                    f"<td>{row['reads_assigned_direct']}</td></tr>"
                )
            out.append("</tbody></table>")
        out.append("</section>")
    out.append("</body></html>")
    path.write_text("\n".join(out), encoding="utf-8")


def main():
    args = parse_args()
    reports = sorted(args.results_dir.rglob("*_kraken_report.txt"))
    summaries = sorted(args.results_dir.rglob("*_kraken_summary.csv"))
    if not reports:
        raise ValueError(f"No *_kraken_report.txt files found under {args.results_dir}")
    if not summaries:
        raise ValueError(f"No *_kraken_summary.csv files found under {args.results_dir}")

    taxa = {
        sample_from_path(path, "_kraken_report.txt"): parse_report(path) for path in reports
    }
    qc = read_qc_summaries(summaries)
    if set(taxa) != set(qc):
        missing_qc = sorted(set(taxa) - set(qc))
        missing_reports = sorted(set(qc) - set(taxa))
        raise ValueError(
            f"Report/summary sample mismatch; missing summaries={missing_qc}, "
            f"missing reports={missing_reports}"
        )
    if any(row["ng"] <= 0 or row["unclassified"] <= 0 for row in qc.values()):
        raise ValueError("Log-scale QC plot requires positive percentages for every sample")

    ng_limits, uc_limits, genus_limits = classify(qc)
    scatter_plot = make_plot(qc, args.batch_name, ng_limits, uc_limits)
    genus_plot = make_genus_only_plot(qc, args.batch_name, genus_limits)
    write_csv(args.csv_output, taxa)
    write_html(
        args.html_output,
        args.batch_name,
        scatter_plot,
        genus_plot,
        taxa,
        qc,
        ng_limits,
        uc_limits,
        genus_limits,
    )


if __name__ == "__main__":
    main()
