from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq


ROOT = Path(__file__).resolve().parents[1]
ALIGNED = (
    ROOT
    / "06_phylogeny_all_sequences"
    / "rooted_reference_outgroup_117"
    / "Nipah_rooted_reference_outgroup_117_aligned.fasta"
)
META = (
    ROOT
    / "07_itol_annotations"
    / "FINAL_iTOL_rooted_reference_outgroup_117"
    / "00_tip_metadata_used_for_itol.csv"
)
OUTDIR = ROOT / "09_country_wise_genome_variants"
GENBANK_ANNOTATIONS = ROOT / "04_reference_sequences" / "country_reference_annotations.gb"

COUNTRY_REFERENCES = {
    "Bangladesh": "AY988601.1",
    "India": "MH396625.1",
    "Malaysia": "AF212302.2",
}

VALID_BASES = set("ACGT")

CANONICAL_CDS_PRODUCTS = {
    "nucleocapsid protein": "N",
    "phosphoprotein": "P",
    "matrix protein": "M",
    "fusion protein": "F",
    "glycoprotein": "G",
    "attachment glycoprotein": "G",
    "l protein": "L",
    "polymerase": "L",
}


def load_alignment(path: Path) -> dict[str, str]:
    return {record.id: str(record.seq).upper() for record in SeqIO.parse(str(path), "fasta")}


def coordinate_map(ref_aligned: str) -> list[int | None]:
    coord = 0
    mapping = []
    for char in ref_aligned:
        if char == "-":
            mapping.append(None)
        else:
            coord += 1
            mapping.append(coord)
    return mapping


def index_by_ref_coord(ref_map: list[int | None]) -> dict[int, int]:
    return {coord: idx for idx, coord in enumerate(ref_map) if coord is not None}


def load_reference_cds_features(path: Path) -> dict[str, list[dict]]:
    if not path.exists():
        raise FileNotFoundError(
            f"Missing {path}. Download GenBank annotations for country references before running."
        )

    features_by_reference = {}
    for record in SeqIO.parse(str(path), "genbank"):
        features = []
        for feature in record.features:
            if feature.type != "CDS":
                continue
            product = feature.qualifiers.get("product", [""])[0].strip()
            gene = CANONICAL_CDS_PRODUCTS.get(product.lower())
            if gene is None:
                continue
            features.append(
                {
                    "gene": gene,
                    "product": product,
                    "start": int(feature.location.start) + 1,
                    "end": int(feature.location.end),
                }
            )
        features_by_reference[record.id] = sorted(features, key=lambda item: item["start"])
    return features_by_reference


def feature_for_position(position: int, cds_features: list[dict]):
    for feature in cds_features:
        if feature["start"] <= position <= feature["end"]:
            return feature
    return None


def translate_codon(codon: str) -> str:
    if len(codon) != 3 or any(base not in VALID_BASES for base in codon):
        return "NA"
    return str(Seq(codon).translate(to_stop=False))


def coding_effect(
    position: int,
    seq_aligned: str,
    ref_aligned: str,
    ref_index: dict[int, int],
    cds_features: list[dict],
):
    feature = feature_for_position(position, cds_features)
    if feature is None:
        return {
            "genomic_region": "noncoding/intergenic",
            "gene": "noncoding",
            "product": "noncoding/intergenic",
            "codon_position_in_gene": "",
            "codon_number": "",
            "ref_codon": "",
            "alt_codon": "",
            "ref_aa": "",
            "alt_aa": "",
            "aa_change": "",
            "variant_effect": "noncoding",
        }

    offset = position - feature["start"]
    codon_start = feature["start"] + (offset // 3) * 3
    coords = [codon_start, codon_start + 1, codon_start + 2]
    if any(coord not in ref_index for coord in coords):
        return {
            "genomic_region": "CDS",
            "gene": feature["gene"],
            "product": feature["product"],
            "codon_position_in_gene": offset % 3 + 1,
            "codon_number": offset // 3 + 1,
            "ref_codon": "",
            "alt_codon": "",
            "ref_aa": "",
            "alt_aa": "",
            "aa_change": "",
            "variant_effect": "coding_unresolved",
        }

    idxs = [ref_index[coord] for coord in coords]
    ref_codon = "".join(ref_aligned[idx] for idx in idxs)
    alt_codon = "".join(seq_aligned[idx] for idx in idxs)
    ref_aa = translate_codon(ref_codon)
    alt_aa = translate_codon(alt_codon)
    if ref_aa == "NA" or alt_aa == "NA":
        effect = "coding_unresolved"
        aa_change = ""
    elif ref_aa == alt_aa:
        effect = "synonymous"
        aa_change = f"{ref_aa}{offset // 3 + 1}{alt_aa}"
    elif alt_aa == "*":
        effect = "nonsense"
        aa_change = f"{ref_aa}{offset // 3 + 1}*"
    else:
        effect = "nonsynonymous"
        aa_change = f"{ref_aa}{offset // 3 + 1}{alt_aa}"

    return {
        "genomic_region": "CDS",
        "gene": feature["gene"],
        "product": feature["product"],
        "codon_position_in_gene": offset % 3 + 1,
        "codon_number": offset // 3 + 1,
        "ref_codon": ref_codon,
        "alt_codon": alt_codon,
        "ref_aa": ref_aa,
        "alt_aa": alt_aa,
        "aa_change": aa_change,
        "variant_effect": effect,
    }


def call_country_snps(
    country: str,
    reference_id: str,
    alignment: dict[str, str],
    metadata: pd.DataFrame,
    cds_features: list[dict],
):
    ref_aligned = alignment[reference_id]
    ref_map = coordinate_map(ref_aligned)
    ref_index = index_by_ref_coord(ref_map)
    country_meta = metadata[
        metadata["country"].eq(country)
        & metadata["tip_id"].isin(alignment.keys())
        & ~metadata["tip_id"].eq(reference_id)
    ].copy()
    meta_by_id = country_meta.set_index("tip_id").to_dict("index")

    rows = []
    qc_rows = []
    for seq_id in country_meta["tip_id"].tolist():
        seq_aligned = alignment[seq_id]
        compared_sites = 0
        ambiguous_or_gap_sites = 0
        snp_count = 0
        transition_count = 0
        transversion_count = 0

        for aln_idx, ref_pos in enumerate(ref_map):
            if ref_pos is None:
                continue
            ref_base = ref_aligned[aln_idx]
            seq_base = seq_aligned[aln_idx]
            if ref_base not in VALID_BASES:
                continue
            if seq_base not in VALID_BASES:
                ambiguous_or_gap_sites += 1
                continue
            compared_sites += 1
            if seq_base == ref_base:
                continue
            snp_count += 1
            change = f"{ref_base}{ref_pos}{seq_base}"
            mut_type = "transition" if {ref_base, seq_base} in [set("AG"), set("CT")] else "transversion"
            if mut_type == "transition":
                transition_count += 1
            else:
                transversion_count += 1
            meta = meta_by_id.get(seq_id, {})
            rows.append(
                {
                    "country": country,
                    "reference_id": reference_id,
                    "sequence_id": seq_id,
                    "host": meta.get("host", "Unknown"),
                    "sample_type": meta.get("sample_type", "Unknown"),
                    "year": meta.get("year", ""),
                    "major_lineage": meta.get("major_lineage", ""),
                    "sequence_category": meta.get("sequence_category", ""),
                    "reference_position": ref_pos,
                    "reference_base": ref_base,
                    "alternate_base": seq_base,
                    "nt_change": change,
                    "mutation_type": mut_type,
                    **coding_effect(ref_pos, seq_aligned, ref_aligned, ref_index, cds_features),
                }
            )

        qc_rows.append(
            {
                "country": country,
                "reference_id": reference_id,
                "sequence_id": seq_id,
                "host": meta_by_id.get(seq_id, {}).get("host", "Unknown"),
                "sample_type": meta_by_id.get(seq_id, {}).get("sample_type", "Unknown"),
                "year": meta_by_id.get(seq_id, {}).get("year", ""),
                "compared_sites": compared_sites,
                "ambiguous_or_gap_reference_sites": ambiguous_or_gap_sites,
                "snp_count": snp_count,
                "transition_count": transition_count,
                "transversion_count": transversion_count,
                "snps_per_10kb": round(snp_count / compared_sites * 10000, 3) if compared_sites else 0,
            }
        )

    return pd.DataFrame(rows), pd.DataFrame(qc_rows)


def recurrent_summary(snp_df: pd.DataFrame, group_sizes: pd.DataFrame):
    if snp_df.empty:
        return pd.DataFrame()
    summary = (
        snp_df.groupby(
            [
                "country",
                "reference_id",
                "reference_position",
                "reference_base",
                "alternate_base",
                "nt_change",
                "genomic_region",
                "gene",
                "product",
                "variant_effect",
                "codon_number",
                "aa_change",
            ],
            dropna=False,
        )["sequence_id"]
        .nunique()
        .reset_index(name="sequence_count")
    )
    summary = summary.merge(group_sizes, on=["country", "reference_id"], how="left")
    summary["percent"] = (summary["sequence_count"] / summary["group_size"] * 100).round(2)
    summary["frequency_class"] = summary.apply(
        lambda row: "high_frequency_recurrent"
        if row["percent"] >= 80
        else ("recurrent" if row["sequence_count"] >= 5 else ("low_frequency" if row["sequence_count"] >= 2 else "singleton")),
        axis=1,
    )
    return summary.sort_values(["country", "sequence_count", "reference_position"], ascending=[True, False, True])


def gene_burden(snp_df: pd.DataFrame, group_sizes: pd.DataFrame):
    if snp_df.empty:
        return pd.DataFrame()
    rows = []
    for (country, reference_id, gene), group in snp_df.groupby(["country", "reference_id", "gene"]):
        group_size = int(group_sizes[(group_sizes["country"].eq(country)) & (group_sizes["reference_id"].eq(reference_id))]["group_size"].iloc[0])
        rows.append(
            {
                "country": country,
                "reference_id": reference_id,
                "gene": gene,
                "unique_snp_sites": group["reference_position"].nunique(),
                "total_sequence_level_snp_events": len(group),
                "mean_snp_events_per_sequence": round(len(group) / group_size, 3) if group_size else 0,
                "synonymous_sites": group[group["variant_effect"].eq("synonymous")]["reference_position"].nunique(),
                "nonsynonymous_sites": group[group["variant_effect"].eq("nonsynonymous")]["reference_position"].nunique(),
                "nonsense_sites": group[group["variant_effect"].eq("nonsense")]["reference_position"].nunique(),
                "noncoding_sites": group[group["variant_effect"].eq("noncoding")]["reference_position"].nunique(),
            }
        )
    return pd.DataFrame(rows).sort_values(["country", "gene"])


def qc_flags(qc_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for country, group in qc_df.groupby("country"):
        high_snp_cutoff = group["snp_count"].quantile(0.95)
        low_coverage_cutoff = group["compared_sites"].quantile(0.05)
        flagged = group[
            group["snp_count"].ge(high_snp_cutoff)
            | group["compared_sites"].le(low_coverage_cutoff)
        ].copy()
        flagged["country_95th_percentile_snp_count"] = round(high_snp_cutoff, 3)
        flagged["country_5th_percentile_compared_sites"] = round(low_coverage_cutoff, 3)
        rows.append(flagged)
    if not rows:
        return pd.DataFrame()
    return pd.concat(rows, ignore_index=True).sort_values(["country", "snp_count"], ascending=[True, False])


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    alignment = load_alignment(ALIGNED)
    metadata = pd.read_csv(META)
    reference_features = load_reference_cds_features(GENBANK_ANNOTATIONS)

    all_snps = []
    all_qc = []
    for country, ref_id in COUNTRY_REFERENCES.items():
        if ref_id not in alignment:
            raise ValueError(f"Reference {ref_id} is missing from alignment")
        if ref_id not in reference_features:
            raise ValueError(f"Reference {ref_id} is missing from {GENBANK_ANNOTATIONS}")
        snps, qc = call_country_snps(country, ref_id, alignment, metadata, reference_features[ref_id])
        country_dir = OUTDIR / country
        country_dir.mkdir(exist_ok=True)
        snps.to_csv(country_dir / f"{country}_SNPs_vs_{ref_id.replace('.', '_')}.csv", index=False)
        qc.to_csv(country_dir / f"{country}_SNP_counts_per_sequence.csv", index=False)
        all_snps.append(snps)
        all_qc.append(qc)

    snp_df = pd.concat(all_snps, ignore_index=True)
    qc_df = pd.concat(all_qc, ignore_index=True)
    group_sizes = (
        qc_df.groupby(["country", "reference_id"])["sequence_id"]
        .nunique()
        .reset_index(name="group_size")
    )
    recurrent = recurrent_summary(snp_df, group_sizes)
    burden = gene_burden(snp_df, group_sizes)
    recurrent_nonsynonymous = recurrent[
        recurrent["variant_effect"].eq("nonsynonymous")
        & recurrent["sequence_count"].ge(2)
    ].copy()
    high_frequency_recurrent = recurrent[
        recurrent["frequency_class"].eq("high_frequency_recurrent")
    ].copy()
    qc_flag_table = qc_flags(qc_df)

    snp_df.to_csv(OUTDIR / "all_country_reference_aware_SNP_events.csv", index=False)
    qc_df.to_csv(OUTDIR / "all_country_SNP_counts_per_sequence.csv", index=False)
    recurrent.to_csv(OUTDIR / "recurrent_SNPs_by_country_reference.csv", index=False)
    recurrent_nonsynonymous.to_csv(OUTDIR / "recurrent_nonsynonymous_SNPs_by_country_reference.csv", index=False)
    high_frequency_recurrent.to_csv(OUTDIR / "high_frequency_recurrent_SNPs_by_country_reference.csv", index=False)
    burden.to_csv(OUTDIR / "gene_level_variant_burden_by_country.csv", index=False)
    group_sizes.to_csv(OUTDIR / "country_variant_group_sizes.csv", index=False)
    qc_flag_table.to_csv(OUTDIR / "QC_flagged_sequences_for_sensitivity_review.csv", index=False)

    readme = [
        "Country-wise genome-wide SNP/variant profile",
        "",
        "Reference-aware comparisons:",
        "- Bangladesh vs AY988601.1",
        "- India vs MH396625.1",
        "- Malaysia vs AF212302.2",
        "",
        "Only clean A/C/G/T SNPs are called. Sites with sequence gaps or ambiguous bases are not counted as SNPs.",
        "Hendra NC_001906.3 is excluded.",
        "",
        "Main outputs:",
        "- all_country_reference_aware_SNP_events.csv",
        "- all_country_SNP_counts_per_sequence.csv",
        "- recurrent_SNPs_by_country_reference.csv",
        "- recurrent_nonsynonymous_SNPs_by_country_reference.csv",
        "- high_frequency_recurrent_SNPs_by_country_reference.csv",
        "- gene_level_variant_burden_by_country.csv",
        "- QC_flagged_sequences_for_sensitivity_review.csv",
        "",
        "Coding annotation uses the matching GenBank CDS coordinates for each country reference.",
        "Variant effects are synonymous, nonsynonymous, nonsense, noncoding, or coding_unresolved.",
        "",
        "Interpretation notes:",
        "- SNP counts are country-reference-aware and should not be compared as prototype-relative lineage markers.",
        "- QC-flagged sequences are not automatically removed; they are candidates for sensitivity analysis.",
        "- Insertions/deletions are not called in this SNP-only output.",
    ]
    (OUTDIR / "README_country_wise_genome_variants.txt").write_text("\n".join(readme) + "\n", encoding="utf-8")

    print("Wrote country-wise genome variant outputs to", OUTDIR)
    print(qc_df.groupby("country")["sequence_id"].nunique().reset_index(name="sequences_compared").to_string(index=False))
    print("\nSNP count summary:")
    print(qc_df.groupby("country")["snp_count"].describe().round(2).to_string())


if __name__ == "__main__":
    main()
