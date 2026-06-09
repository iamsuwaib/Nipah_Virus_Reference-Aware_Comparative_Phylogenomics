from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
VARIANT_DIR = ROOT / "09_country_wise_genome_variants"
OUTDIR = VARIANT_DIR / "simple_significant_summary"
STATS_DIR = VARIANT_DIR / "statistical_significance"

COUNTRY_COLORS = {
    "Bangladesh": "#1f77b4",
    "India": "#ff7f0e",
    "Malaysia": "#2ca02c",
}

REFERENCE_ARTIFACTS = pd.DataFrame(
    [
        {
            "country": "Bangladesh",
            "reference_id": "AY988601.1",
            "gene": "F",
            "aa_change": "S207L",
            "interpretation": "Reference-sequence correction, not a biological variant for main-text interpretation.",
            "literature_note": "Lo et al. Scientific Reports 2019 reported S207L as one of two changes correcting amino-acid sequence errors in deposited AY988601.",
        },
        {
            "country": "Bangladesh",
            "reference_id": "AY988601.1",
            "gene": "F",
            "aa_change": "G252D",
            "interpretation": "Reference-sequence correction, not a biological variant for main-text interpretation.",
            "literature_note": "Lo et al. Scientific Reports 2019 reported G252D as one of two changes correcting amino-acid sequence errors in deposited AY988601.",
        },
    ]
)


def save_table(df: pd.DataFrame, name: str) -> Path:
    path = OUTDIR / name
    df.to_csv(path, index=False)
    return path


def p_label(p_value: float) -> str:
    if pd.isna(p_value):
        return "P = NA"
    if p_value < 0.001:
        return "P < 0.001"
    if p_value < 0.01:
        return "P < 0.01"
    if p_value < 0.05:
        return "P < 0.05"
    return "P > 0.05"


def add_sig_bracket(ax, x1, x2, y, h, label):
    ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], color="#222222", lw=1.2)
    ax.text((x1 + x2) / 2, y + h + 3, label, ha="center", va="bottom", fontsize=10)


def load_pairwise_snp_pvalues() -> dict[str, float]:
    path = STATS_DIR / "02_snp_burden_pairwise_mannwhitney_fdr.csv"
    if not path.exists():
        return {}
    stats = pd.read_csv(path)
    return dict(zip(stats["comparison"], stats["p_fdr"]))


def load_variant_enrichment_pvalues() -> dict[tuple[str, str], float]:
    path = STATS_DIR / "07_significant_recurrent_variant_country_enrichment.csv"
    if not path.exists():
        return {}
    stats = pd.read_csv(path)
    return {
        (row["focal_country"], row["variant"]): row["p_fdr"]
        for _, row in stats.iterrows()
    }


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)

    counts = pd.read_csv(VARIANT_DIR / "all_country_SNP_counts_per_sequence.csv")
    burden = pd.read_csv(VARIANT_DIR / "gene_level_variant_burden_by_country.csv")
    recurrent_ns = pd.read_csv(VARIANT_DIR / "recurrent_nonsynonymous_SNPs_by_country_reference.csv")
    high_freq = pd.read_csv(VARIANT_DIR / "high_frequency_recurrent_SNPs_by_country_reference.csv")
    qc = pd.read_csv(VARIANT_DIR / "QC_flagged_sequences_for_sensitivity_review.csv")
    pairwise_snp_pvalues = load_pairwise_snp_pvalues()
    variant_pvalues = load_variant_enrichment_pvalues()

    country_summary = (
        counts.groupby(["country", "reference_id"])
        .agg(
            sequences=("sequence_id", "nunique"),
            median_snps_per_sequence=("snp_count", "median"),
            mean_snps_per_sequence=("snp_count", "mean"),
            min_snps=("snp_count", "min"),
            max_snps=("snp_count", "max"),
            median_snps_per_10kb=("snps_per_10kb", "median"),
        )
        .reset_index()
    )
    country_summary["mean_snps_per_sequence"] = country_summary["mean_snps_per_sequence"].round(2)
    country_summary["median_snps_per_10kb"] = country_summary["median_snps_per_10kb"].round(2)
    country_summary = country_summary.sort_values("mean_snps_per_sequence", ascending=False)
    save_table(country_summary, "01_country_SNP_burden_simple_summary.csv")

    coding_burden = burden[burden["gene"].ne("noncoding")].copy()
    gene_signal = coding_burden[
        [
            "country",
            "reference_id",
            "gene",
            "unique_snp_sites",
            "mean_snp_events_per_sequence",
            "synonymous_sites",
            "nonsynonymous_sites",
            "nonsense_sites",
        ]
    ].sort_values(["country", "mean_snp_events_per_sequence"], ascending=[True, False])
    save_table(gene_signal, "02_gene_level_signal_simple_summary.csv")

    top_recurrent_ns = recurrent_ns[
        [
            "country",
            "reference_id",
            "nt_change",
            "gene",
            "aa_change",
            "sequence_count",
            "percent",
            "frequency_class",
        ]
    ].sort_values(["country", "sequence_count"], ascending=[True, False])
    top_recurrent_ns = top_recurrent_ns.groupby("country", group_keys=False).head(12)
    save_table(top_recurrent_ns, "03_top_recurrent_nonsynonymous_variants.csv")

    high_freq_simple = high_freq[
        [
            "country",
            "reference_id",
            "nt_change",
            "gene",
            "variant_effect",
            "aa_change",
            "sequence_count",
            "percent",
        ]
    ].sort_values(["country", "sequence_count"], ascending=[True, False])
    save_table(high_freq_simple, "04_high_frequency_recurrent_variants.csv")
    save_table(REFERENCE_ARTIFACTS, "00_reference_artifacts_do_not_interpret_as_biology.csv")

    recurrent_ns_biological = recurrent_ns.merge(
        REFERENCE_ARTIFACTS[["country", "reference_id", "gene", "aa_change"]],
        on=["country", "reference_id", "gene", "aa_change"],
        how="left",
        indicator=True,
    )
    recurrent_ns_biological = recurrent_ns_biological[recurrent_ns_biological["_merge"].eq("left_only")].drop(columns="_merge")
    top_recurrent_ns_biological = recurrent_ns_biological[
        [
            "country",
            "reference_id",
            "nt_change",
            "gene",
            "aa_change",
            "sequence_count",
            "percent",
            "frequency_class",
        ]
    ].sort_values(["country", "sequence_count"], ascending=[True, False])
    top_recurrent_ns_biological = top_recurrent_ns_biological.groupby("country", group_keys=False).head(12)
    save_table(top_recurrent_ns_biological, "03b_top_recurrent_nonsynonymous_variants_reference_artifacts_removed.csv")

    qc_simple = qc[
        [
            "country",
            "reference_id",
            "sequence_id",
            "compared_sites",
            "snp_count",
            "snps_per_10kb",
            "country_95th_percentile_snp_count",
            "country_5th_percentile_compared_sites",
        ]
    ].sort_values(["country", "snp_count"], ascending=[True, False])
    save_table(qc_simple, "05_QC_sequences_to_review_before_final_claims.csv")

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "axes.titlesize": 15,
            "axes.labelsize": 12,
            "xtick.labelsize": 11,
            "ytick.labelsize": 10,
        }
    )

    order = ["Bangladesh", "India", "Malaysia"]
    fig, ax = plt.subplots(figsize=(8.6, 6.0))
    values = [counts[counts["country"].eq(country)]["snp_count"].values for country in order]
    box = ax.boxplot(values, patch_artist=True, labels=order, showfliers=True)
    for patch, country in zip(box["boxes"], order):
        patch.set_facecolor(COUNTRY_COLORS[country])
        patch.set_alpha(0.75)
    ax.set_title("Genome-wide SNP burden by country-specific reference", fontsize=15, pad=26)
    ax.set_ylabel("SNPs per sequence")
    ax.set_ylim(-20, 535)
    ax.grid(axis="y", alpha=0.25)
    ax.spines[["top", "right"]].set_visible(False)
    add_sig_bracket(
        ax,
        1,
        2,
        425,
        8,
        p_label(pairwise_snp_pvalues.get("Bangladesh vs India", float("nan"))),
    )
    add_sig_bracket(
        ax,
        2,
        3,
        460,
        8,
        p_label(pairwise_snp_pvalues.get("India vs Malaysia", float("nan"))),
    )
    add_sig_bracket(
        ax,
        1,
        3,
        495,
        8,
        p_label(pairwise_snp_pvalues.get("Bangladesh vs Malaysia", float("nan"))),
    )
    fig.subplots_adjust(top=0.82, left=0.12, right=0.98, bottom=0.12)
    fig.savefig(OUTDIR / "Figure_1_country_SNP_burden_boxplot.png", dpi=300)
    fig.savefig(OUTDIR / "Figure_1_country_SNP_burden_boxplot.pdf")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9.4, 5.5))
    ns = coding_burden.pivot(index="gene", columns="country", values="nonsynonymous_sites").fillna(0)
    ns = ns.reindex(["N", "P", "M", "F", "G", "L"])
    ns[order].plot(kind="bar", ax=ax, color=[COUNTRY_COLORS[country] for country in order], width=0.78)
    ax.set_title("Nonsynonymous SNP sites by gene", fontsize=15, pad=18)
    ax.set_xlabel("Gene")
    ax.set_ylabel("Unique nonsynonymous SNP sites")
    ax.set_ylim(0, 86)
    ax.grid(axis="y", alpha=0.25)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(title="Country", frameon=False)
    fig.tight_layout()
    fig.savefig(OUTDIR / "Figure_2_gene_nonsynonymous_variant_burden.png", dpi=300)
    fig.savefig(OUTDIR / "Figure_2_gene_nonsynonymous_variant_burden.pdf")
    plt.close(fig)

    top_plot = top_recurrent_ns.copy()
    top_plot["label"] = top_plot["country"] + " " + top_plot["gene"] + ":" + top_plot["aa_change"].fillna("")
    top_plot = top_plot.sort_values("percent").tail(18)
    fig, ax = plt.subplots(figsize=(9.6, 6.8))
    colors = [COUNTRY_COLORS.get(country, "#777777") for country in top_plot["country"]]
    ax.barh(top_plot["label"], top_plot["percent"], color=colors, alpha=0.85)
    ax.set_title("Most recurrent nonsynonymous variants")
    ax.set_xlabel("Sequences carrying variant (%)")
    ax.set_xlim(0, 105)
    ax.grid(axis="x", alpha=0.25)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUTDIR / "Figure_3_top_recurrent_nonsynonymous_variants.png", dpi=300)
    fig.savefig(OUTDIR / "Figure_3_top_recurrent_nonsynonymous_variants.pdf")
    plt.close(fig)

    top_plot = top_recurrent_ns_biological.copy()
    top_plot["label"] = top_plot["country"] + " " + top_plot["gene"] + ":" + top_plot["aa_change"].fillna("")
    top_plot["variant"] = top_plot["gene"] + ":" + top_plot["aa_change"].fillna("")
    top_plot["p_fdr_for_plot"] = top_plot.apply(
        lambda row: variant_pvalues.get((row["country"], row["variant"]), float("nan")),
        axis=1,
    )
    top_plot = top_plot.sort_values("percent").tail(18)
    fig, ax = plt.subplots(figsize=(9.6, 6.8))
    colors = [COUNTRY_COLORS.get(country, "#777777") for country in top_plot["country"]]
    bars = ax.barh(top_plot["label"], top_plot["percent"], color=colors, alpha=0.85)
    ax.set_title("Most recurrent nonsynonymous variants, reference artifacts removed")
    ax.set_xlabel("Sequences carrying variant (%)")
    ax.set_xlim(0, 122)
    for bar, p_value in zip(bars, top_plot["p_fdr_for_plot"]):
        ax.text(
            bar.get_width() + 1.2,
            bar.get_y() + bar.get_height() / 2,
            p_label(p_value),
            ha="left",
            va="center",
            fontsize=8.5,
        )
    ax.grid(axis="x", alpha=0.25)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUTDIR / "Figure_3b_top_recurrent_nonsynonymous_variants_reference_artifacts_removed.png", dpi=300)
    fig.savefig(OUTDIR / "Figure_3b_top_recurrent_nonsynonymous_variants_reference_artifacts_removed.pdf")
    plt.close(fig)

    narrative = [
        "Simple significant summary of country-wise genome-wide SNP analysis",
        "",
        "Main comparisons:",
        "- Bangladesh sequences vs AY988601.1",
        "- India sequences vs MH396625.1",
        "- Malaysia sequences vs AF212302.2",
        "",
        "Most important takeaways:",
        "1. Bangladesh and India show substantially higher within-country SNP burden than Malaysia when each is compared to its own reference.",
        "2. Malaysia sequences are mostly very close to AF212302.2; AJ627196.1 is the main Malaysia sequence flagged for review.",
        "3. Bangladesh F G252D and S207L are reference-sequence correction artifacts reported for AY988601.1, not biological variants to emphasize.",
        "4. After removing those reference artifacts, recurrent nonsynonymous signals remain in Bangladesh N/P/L/G/F and India P/L/G.",
        "5. India has a high-SNP tail driven especially by FJ513078.1 and PP554504.1, so those should be included in sensitivity analysis before final biological wording.",
        "6. These outputs are country-reference-aware. They are not prototype-relative lineage marker tables.",
        "",
        "Literature interpretation:",
        "- Bangladesh/India and Malaysia represent distinct major NiV genotypes/lineages in prior phylogenetic work.",
        "- Malaysia outbreak genomes are reported as highly similar to one another, matching the low SNP burden observed here.",
        "- Bangladesh has documented co-circulating lineages and frequent spillover, matching higher within-country diversity.",
        "- Kerala/India sequences have been reported as B genotype/Bangladesh-lineage-related but locally distinct, matching our India distribution and outliers.",
        "- A 2024 human-isolate comparative genomics paper reported P as highly variable and M as relatively stable, matching our gene-level nonsynonymous burden.",
        "",
        "Recommended paper-use outputs:",
        "- Figure_1_country_SNP_burden_boxplot.png",
        "- Figure_2_gene_nonsynonymous_variant_burden.png",
        "- Figure_3b_top_recurrent_nonsynonymous_variants_reference_artifacts_removed.png",
        "- 03b_top_recurrent_nonsynonymous_variants_reference_artifacts_removed.csv",
        "- 00_reference_artifacts_do_not_interpret_as_biology.csv",
        "- 05_QC_sequences_to_review_before_final_claims.csv",
        "",
        "P-value notes:",
        "- Figure 1 uses FDR-adjusted pairwise Mann-Whitney U tests.",
        "- Figure 2 is descriptive; gene-level Fisher tests are kept in the statistical_significance folder to avoid implying pairwise country comparisons.",
        "- Figure 3b uses FDR-adjusted Fisher exact tests for country enrichment of each recurrent nonsynonymous variant.",
    ]
    (OUTDIR / "README_simple_significant_summary.txt").write_text("\n".join(narrative) + "\n", encoding="utf-8")

    print(f"Wrote simplified summary to {OUTDIR}")
    print(country_summary.to_string(index=False))


if __name__ == "__main__":
    main()
