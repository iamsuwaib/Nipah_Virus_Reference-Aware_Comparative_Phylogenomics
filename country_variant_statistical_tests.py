from itertools import combinations
from pathlib import Path

import pandas as pd
from scipy.stats import fisher_exact, kruskal, mannwhitneyu


ROOT = Path(__file__).resolve().parents[1]
VARIANT_DIR = ROOT / "09_country_wise_genome_variants"
SIMPLE_DIR = VARIANT_DIR / "simple_significant_summary"
OUTDIR = VARIANT_DIR / "statistical_significance"

COUNTRIES = ["Bangladesh", "India", "Malaysia"]

REFERENCE_ARTIFACTS = {
    ("Bangladesh", "AY988601.1", "F", "S207L"),
    ("Bangladesh", "AY988601.1", "F", "G252D"),
}


def bh_fdr(pvalues):
    values = pd.Series(pvalues, dtype="float64")
    ranked = values.rank(method="first")
    n = values.notna().sum()
    adjusted = values * n / ranked
    adjusted = adjusted.sort_values(ascending=False).cummin().sort_index()
    return adjusted.clip(upper=1.0)


def significance_label(p):
    if pd.isna(p):
        return "NA"
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def snp_burden_tests(counts: pd.DataFrame, label: str):
    rows = []
    groups = [counts[counts["country"].eq(country)]["snp_count"] for country in COUNTRIES]
    h_stat, p_value = kruskal(*groups)
    omnibus = pd.DataFrame(
        [
            {
                "analysis": label,
                "test": "Kruskal-Wallis",
                "comparison": "Bangladesh vs India vs Malaysia",
                "statistic": h_stat,
                "p_value": p_value,
                "p_fdr": p_value,
                "significance": significance_label(p_value),
            }
        ]
    )

    for country_a, country_b in combinations(COUNTRIES, 2):
        a = counts[counts["country"].eq(country_a)]["snp_count"]
        b = counts[counts["country"].eq(country_b)]["snp_count"]
        stat, p = mannwhitneyu(a, b, alternative="two-sided")
        rows.append(
            {
                "analysis": label,
                "test": "Mann-Whitney U",
                "comparison": f"{country_a} vs {country_b}",
                "statistic": stat,
                "p_value": p,
            }
        )
    pairwise = pd.DataFrame(rows)
    pairwise["p_fdr"] = bh_fdr(pairwise["p_value"])
    pairwise["significance"] = pairwise["p_fdr"].map(significance_label)
    return omnibus, pairwise


def variant_enrichment_tests(events: pd.DataFrame, counts: pd.DataFrame):
    recurrent = pd.read_csv(SIMPLE_DIR / "03b_top_recurrent_nonsynonymous_variants_reference_artifacts_removed.csv")
    recurrent = recurrent.drop_duplicates(["gene", "aa_change"])
    all_sequences = counts[["country", "sequence_id"]].drop_duplicates()
    country_sizes = all_sequences.groupby("country")["sequence_id"].nunique().to_dict()
    rows = []

    for _, variant in recurrent.iterrows():
        key = (variant["country"], variant["reference_id"], variant["gene"], variant["aa_change"])
        if key in REFERENCE_ARTIFACTS:
            continue
        variant_events = events[
            events["gene"].eq(variant["gene"])
            & events["aa_change"].eq(variant["aa_change"])
            & events["variant_effect"].eq("nonsynonymous")
        ][["country", "sequence_id"]].drop_duplicates()

        for focal_country in COUNTRIES:
            focal_total = country_sizes.get(focal_country, 0)
            other_total = sum(size for country, size in country_sizes.items() if country != focal_country)
            focal_present = variant_events[variant_events["country"].eq(focal_country)]["sequence_id"].nunique()
            other_present = variant_events[~variant_events["country"].eq(focal_country)]["sequence_id"].nunique()
            focal_absent = focal_total - focal_present
            other_absent = other_total - other_present
            if min(focal_total, other_total) == 0:
                continue
            odds_ratio, p = fisher_exact(
                [[focal_present, focal_absent], [other_present, other_absent]],
                alternative="two-sided",
            )
            rows.append(
                {
                    "variant": f"{variant['gene']}:{variant['aa_change']}",
                    "gene": variant["gene"],
                    "aa_change": variant["aa_change"],
                    "focal_country": focal_country,
                    "focal_present": focal_present,
                    "focal_total": focal_total,
                    "other_present": other_present,
                    "other_total": other_total,
                    "odds_ratio": odds_ratio,
                    "p_value": p,
                }
            )

    df = pd.DataFrame(rows)
    df["p_fdr"] = bh_fdr(df["p_value"])
    df["significance"] = df["p_fdr"].map(significance_label)
    df["focal_percent"] = (df["focal_present"] / df["focal_total"] * 100).round(2)
    df["other_percent"] = (df["other_present"] / df["other_total"] * 100).round(2)
    df["enrichment_direction"] = df.apply(
        lambda row: "enriched"
        if row["focal_percent"] > row["other_percent"]
        else ("depleted" if row["focal_percent"] < row["other_percent"] else "equal"),
        axis=1,
    )
    return df.sort_values(["p_fdr", "p_value", "variant"])


def gene_effect_tests(events: pd.DataFrame):
    coding = events[events["variant_effect"].isin(["synonymous", "nonsynonymous"])].copy()
    rows = []
    for country in COUNTRIES:
        country_events = coding[coding["country"].eq(country)]
        for gene in ["N", "P", "M", "F", "G", "L"]:
            gene_events = country_events[country_events["gene"].eq(gene)]
            other_events = country_events[~country_events["gene"].eq(gene)]
            gene_ns = gene_events[gene_events["variant_effect"].eq("nonsynonymous")]["reference_position"].nunique()
            gene_syn = gene_events[gene_events["variant_effect"].eq("synonymous")]["reference_position"].nunique()
            other_ns = other_events[other_events["variant_effect"].eq("nonsynonymous")]["reference_position"].nunique()
            other_syn = other_events[other_events["variant_effect"].eq("synonymous")]["reference_position"].nunique()
            odds_ratio, p = fisher_exact([[gene_ns, gene_syn], [other_ns, other_syn]], alternative="two-sided")
            rows.append(
                {
                    "country": country,
                    "gene": gene,
                    "test": "Fisher exact",
                    "comparison": f"{gene} nonsynonymous:synonymous vs other genes",
                    "gene_nonsynonymous_sites": gene_ns,
                    "gene_synonymous_sites": gene_syn,
                    "other_nonsynonymous_sites": other_ns,
                    "other_synonymous_sites": other_syn,
                    "odds_ratio": odds_ratio,
                    "p_value": p,
                }
            )
    df = pd.DataFrame(rows)
    df["p_fdr"] = bh_fdr(df["p_value"])
    df["significance"] = df["p_fdr"].map(significance_label)
    return df.sort_values(["country", "p_fdr", "gene"])


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    counts = pd.read_csv(VARIANT_DIR / "all_country_SNP_counts_per_sequence.csv")
    events = pd.read_csv(VARIANT_DIR / "all_country_reference_aware_SNP_events.csv")
    qc = pd.read_csv(VARIANT_DIR / "QC_flagged_sequences_for_sensitivity_review.csv")

    omnibus, pairwise = snp_burden_tests(counts, "all_sequences")
    omnibus.to_csv(OUTDIR / "01_snp_burden_omnibus_kruskal_wallis.csv", index=False)
    pairwise.to_csv(OUTDIR / "02_snp_burden_pairwise_mannwhitney_fdr.csv", index=False)

    flagged_ids = set(qc["sequence_id"])
    sensitivity_counts = counts[~counts["sequence_id"].isin(flagged_ids)].copy()
    sens_omnibus, sens_pairwise = snp_burden_tests(sensitivity_counts, "QC_flagged_sequences_removed")
    sens_omnibus.to_csv(OUTDIR / "03_snp_burden_sensitivity_omnibus_kruskal_wallis.csv", index=False)
    sens_pairwise.to_csv(OUTDIR / "04_snp_burden_sensitivity_pairwise_mannwhitney_fdr.csv", index=False)

    gene_tests = gene_effect_tests(events)
    gene_tests.to_csv(OUTDIR / "05_gene_synonymous_nonsynonymous_fisher_tests_fdr.csv", index=False)

    variant_tests = variant_enrichment_tests(events, counts)
    variant_tests.to_csv(OUTDIR / "06_recurrent_variant_country_enrichment_fisher_fdr.csv", index=False)
    variant_tests[
        variant_tests["p_fdr"].lt(0.05) & variant_tests["enrichment_direction"].eq("enriched")
    ].to_csv(
        OUTDIR / "07_significant_recurrent_variant_country_enrichment.csv", index=False
    )

    readme = [
        "Statistical significance tests for country-wise NiV variant analysis",
        "",
        "Tests performed:",
        "1. Kruskal-Wallis test for SNP burden across Bangladesh, India, and Malaysia.",
        "2. Pairwise Mann-Whitney U tests with Benjamini-Hochberg FDR correction.",
        "3. Sensitivity repeat after removing QC-flagged high-SNP/low-compared-site sequences.",
        "4. Fisher exact tests for synonymous:nonsynonymous structure by gene within country.",
        "5. Fisher exact tests for country enrichment of recurrent nonsynonymous variants.",
        "",
        "Interpretation cautions:",
        "- P values test distributional/enrichment patterns only.",
        "- They do not prove functional effect, virulence, transmission advantage, or host adaptation.",
        "- Bangladesh AY988601.1 F S207L and G252D are excluded as known reference-sequence correction artifacts.",
    ]
    (OUTDIR / "README_statistical_tests.txt").write_text("\n".join(readme) + "\n", encoding="utf-8")

    print(f"Wrote statistical tests to {OUTDIR}")
    print("\nSNP burden omnibus:")
    print(omnibus.to_string(index=False))
    print("\nSNP burden pairwise:")
    print(pairwise.to_string(index=False))
    print("\nSensitivity pairwise:")
    print(sens_pairwise.to_string(index=False))
    print("\nTop significant recurrent variant enrichments:")
    print(
        variant_tests[
            variant_tests["p_fdr"].lt(0.05) & variant_tests["enrichment_direction"].eq("enriched")
        ]
        .head(12)
        .to_string(index=False)
    )


if __name__ == "__main__":
    main()
