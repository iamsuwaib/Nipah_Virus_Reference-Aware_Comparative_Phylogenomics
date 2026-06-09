from pathlib import Path
from collections import Counter, defaultdict

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
ITOL_METADATA = (
    ROOT
    / "07_itol_annotations"
    / "FINAL_iTOL_rooted_reference_outgroup_117"
    / "00_tip_metadata_used_for_itol.csv"
)
OUTDIR = ROOT / "08_reports_and_manuscript" / "analysis_outputs" / "FG_protein_analysis"


# NC_002728/AF212302 reference coordinates, 1-based inclusive.
# The rooted dataset excludes NC_002728.1 and keeps AF212302.2 as the
# identical Malaysia/prototype reference representative.
GENE_COORDS = {
    # GenBank CDS coordinates from AF212302.2 / NC_002728.1.
    # Biopython reports F as [6653:8294](+) and G as [8942:10751](+),
    # which converts to 1-based inclusive coordinates below.
    "F": (6654, 8294),
    "G": (8943, 10751),
}
REFERENCE_ID = "AF212302.2"

# Country/genotype-aware reference definitions. Hendra (NC_001906.3) is only
# a phylogenetic outgroup and is intentionally not used as a Nipah protein
# substitution reference.
REFERENCE_COMPARISONS = {
    "Bangladesh_vs_AY988601": {
        "reference_id": "AY988601.1",
        "country": "Bangladesh",
        "description": "Bangladesh sequences compared with Bangladesh reference AY988601.1",
    },
    "India_vs_MH396625": {
        "reference_id": "MH396625.1",
        "country": "India",
        "description": "India sequences compared with India/Kerala reference MH396625.1",
    },
    "Malaysia_vs_AF212302": {
        "reference_id": "AF212302.2",
        "country": "Malaysia",
        "description": "Malaysia sequences compared with Malaysia/prototype reference AF212302.2",
    },
}

REFERENCE_PROTEIN_IDS = {
    "AF212302.2": "Malaysia/prototype",
    "AY988601.1": "Bangladesh",
    "MH396625.1": "India/Kerala",
    "AJ627196.1": "Malaysia pig",
    "MK801755.1": "Cambodia bat",
}


def load_alignment(path: Path):
    return {record.id: str(record.seq).upper() for record in SeqIO.parse(str(path), "fasta")}


def reference_coordinate_map(ref_aligned: str):
    coord = 0
    mapping = []
    for char in ref_aligned:
        if char == "-":
            mapping.append(None)
        else:
            coord += 1
            mapping.append(coord)
    return mapping


def extract_gene_from_alignment(aligned_seq: str, ref_map: list[int | None], start: int, end: int) -> str:
    chars = [
        aligned_seq[i]
        for i, coord in enumerate(ref_map)
        if coord is not None and start <= coord <= end
    ]
    return "".join(chars).replace("-", "")


def translate_cds(nucleotide: str) -> str:
    trimmed_len = len(nucleotide) - (len(nucleotide) % 3)
    trimmed = nucleotide[:trimmed_len]
    if not trimmed:
        return ""
    return str(Seq(trimmed).translate(to_stop=False))


def usable_protein(protein: str) -> bool:
    if not protein:
        return False
    if "X" in protein:
        return False
    if "*" in protein[:-1]:
        return False
    return True


def write_protein_fasta(records: dict[str, str], path: Path):
    lines = []
    for seq_id, protein in records.items():
        lines.append(f">{seq_id}")
        for i in range(0, len(protein), 80):
            lines.append(protein[i : i + 80])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def substitution_table(
    gene: str,
    proteins: dict[str, str],
    metadata: pd.DataFrame,
    reference_id: str,
    sequence_ids: set[str] | None = None,
):
    ref = proteins[reference_id]
    meta_by_tip = metadata.set_index("tip_id").to_dict("index")
    rows = []
    for seq_id, protein in proteins.items():
        if seq_id == reference_id:
            continue
        if sequence_ids is not None and seq_id not in sequence_ids:
            continue
        max_len = min(len(ref), len(protein))
        for idx in range(max_len):
            ref_aa = ref[idx]
            aa = protein[idx]
            if aa == ref_aa:
                continue
            if aa in {"X", "-"} or ref_aa in {"X", "-"}:
                continue
            meta = meta_by_tip.get(seq_id, {})
            rows.append(
                {
                    "gene": gene,
                    "reference_id": reference_id,
                    "sequence_id": seq_id,
                    "aa_position": idx + 1,
                    "reference_aa": ref_aa,
                    "alternate_aa": aa,
                    "substitution": f"{ref_aa}{idx + 1}{aa}",
                    "country": meta.get("country", "Unknown"),
                    "host": meta.get("host", "Unknown"),
                    "major_lineage": meta.get("major_lineage", "Unknown"),
                    "sequence_category": meta.get("sequence_category", "Study sequence"),
                }
            )
    return pd.DataFrame(rows)


def pairwise_reference_differences(gene: str, proteins: dict[str, str]):
    rows = []
    refs = [ref_id for ref_id in REFERENCE_PROTEIN_IDS if ref_id in proteins]
    for i, ref_a in enumerate(refs):
        for ref_b in refs[i + 1 :]:
            protein_a = proteins[ref_a]
            protein_b = proteins[ref_b]
            for idx in range(min(len(protein_a), len(protein_b))):
                aa_a = protein_a[idx]
                aa_b = protein_b[idx]
                if aa_a == aa_b or aa_a in {"X", "-"} or aa_b in {"X", "-"}:
                    continue
                rows.append(
                    {
                        "gene": gene,
                        "reference_a": ref_a,
                        "reference_a_label": REFERENCE_PROTEIN_IDS[ref_a],
                        "reference_b": ref_b,
                        "reference_b_label": REFERENCE_PROTEIN_IDS[ref_b],
                        "aa_position": idx + 1,
                        "reference_a_aa": aa_a,
                        "reference_b_aa": aa_b,
                        "difference": f"{aa_a}{idx + 1}{aa_b}",
                    }
                )
    return pd.DataFrame(rows)


def country_specific_summary(subs: pd.DataFrame, group_size: int):
    if subs.empty:
        return pd.DataFrame()
    rows = []
    for (gene, reference_id, substitution), group in subs.groupby(["gene", "reference_id", "substitution"]):
        n = group["sequence_id"].nunique()
        rows.append(
            {
                "gene": gene,
                "reference_id": reference_id,
                "substitution": substitution,
                "aa_position": int(group["aa_position"].iloc[0]),
                "reference_aa": group["reference_aa"].iloc[0],
                "alternate_aa": group["alternate_aa"].iloc[0],
                "sequence_count": n,
                "group_size": group_size,
                "percent": round(n / group_size * 100, 2) if group_size else 0,
            }
        )
    return pd.DataFrame(rows).sort_values(["gene", "aa_position", "substitution"])


def summarize_substitutions(subs: pd.DataFrame, metadata: pd.DataFrame, usable_by_gene: dict[str, set[str]]):
    meta_by_gene = {
        gene: metadata[metadata["tip_id"].isin(usable_ids)].copy()
        for gene, usable_ids in usable_by_gene.items()
    }

    lineage_rows = []
    for (gene, substitution), group in subs.groupby(["gene", "substitution"]):
        lineage_counts = meta_by_gene[gene]["major_lineage"].value_counts().to_dict()
        counts = group.groupby("major_lineage")["sequence_id"].nunique().to_dict()
        row = {
            "gene": gene,
            "substitution": substitution,
            "aa_position": int(group["aa_position"].iloc[0]),
            "reference_aa": group["reference_aa"].iloc[0],
            "alternate_aa": group["alternate_aa"].iloc[0],
        }
        for lineage, total in lineage_counts.items():
            n = counts.get(lineage, 0)
            row[f"{lineage}_count"] = n
            row[f"{lineage}_percent"] = round(n / total * 100, 2) if total else 0
        lineage_rows.append(row)

    country_rows = []
    for (gene, substitution), group in subs.groupby(["gene", "substitution"]):
        country_counts = meta_by_gene[gene]["country"].value_counts().to_dict()
        counts = group.groupby("country")["sequence_id"].nunique().to_dict()
        row = {
            "gene": gene,
            "substitution": substitution,
            "aa_position": int(group["aa_position"].iloc[0]),
            "reference_aa": group["reference_aa"].iloc[0],
            "alternate_aa": group["alternate_aa"].iloc[0],
        }
        for country, total in country_counts.items():
            n = counts.get(country, 0)
            row[f"{country}_count"] = n
            row[f"{country}_percent"] = round(n / total * 100, 2) if total else 0
        country_rows.append(row)

    return pd.DataFrame(lineage_rows), pd.DataFrame(country_rows)


def lineage_enriched(lineage_summary: pd.DataFrame):
    rows = []
    bi = "Bangladesh/India-enriched NiV"
    my = "Malaysia/prototype-related NiV"
    for _, row in lineage_summary.iterrows():
        bi_pct = float(row.get(f"{bi}_percent", 0))
        my_pct = float(row.get(f"{my}_percent", 0))
        bi_count = int(row.get(f"{bi}_count", 0))
        my_count = int(row.get(f"{my}_count", 0))
        label = ""
        if bi_pct >= 70 and my_pct <= 10 and bi_count >= 5:
            label = "Bangladesh/India-enriched"
        elif my_pct >= 70 and bi_pct <= 10 and my_count >= 5:
            label = "Malaysia/prototype-enriched"
        elif abs(bi_pct - my_pct) >= 50 and (bi_count + my_count) >= 5:
            label = "lineage-skewed"
        if label:
            out = row.to_dict()
            out["enrichment_label"] = label
            out["Bangladesh_India_minus_Malaysia_percent"] = round(bi_pct - my_pct, 2)
            rows.append(out)
    return pd.DataFrame(rows)


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    alignment = load_alignment(ALIGNED)
    metadata = pd.read_csv(ITOL_METADATA)
    metadata = metadata[metadata["tip_id"].isin(alignment.keys())].copy()

    if REFERENCE_ID not in alignment:
        raise RuntimeError(f"{REFERENCE_ID} not found in alignment.")

    ref_map = reference_coordinate_map(alignment[REFERENCE_ID])
    gene_qc_rows = []
    all_subs = []
    usable_by_gene = {}

    for gene, (start, end) in GENE_COORDS.items():
        nuc_records = {}
        protein_records = {}
        qc_rows = []

        for seq_id, aligned_seq in alignment.items():
            if seq_id == "NC_001906.3":
                # Hendra is useful for rooting, not for Nipah F/G substitution calls.
                continue
            cds = extract_gene_from_alignment(aligned_seq, ref_map, start, end)
            protein = translate_cds(cds)
            nuc_records[seq_id] = cds
            protein_records[seq_id] = protein
            qc_rows.append(
                {
                    "gene": gene,
                    "sequence_id": seq_id,
                    "cds_length_nt": len(cds),
                    "protein_length_aa": len(protein.rstrip("*")),
                    "has_terminal_stop": protein.endswith("*"),
                    "internal_stop_count": protein[:-1].count("*"),
                    "ambiguous_nt_count": sum(c not in "ACGT" for c in cds),
                    "usable_for_substitution_analysis": usable_protein(protein),
                }
            )

        usable = {
            seq_id: protein.rstrip("*")
            for seq_id, protein in protein_records.items()
            if usable_protein(protein)
        }
        usable_by_gene[gene] = set(usable.keys())
        write_protein_fasta(usable, OUTDIR / f"{gene}_proteins_usable.fasta")
        pd.DataFrame(qc_rows).to_csv(OUTDIR / f"{gene}_extraction_translation_qc.csv", index=False)
        subs = substitution_table(gene, usable, metadata, REFERENCE_ID)
        subs.to_csv(OUTDIR / f"{gene}_aa_substitutions_vs_AF212302.csv", index=False)
        all_subs.append(subs)

        ref_diffs = pairwise_reference_differences(gene, usable)
        ref_diffs.to_csv(OUTDIR / f"{gene}_reference_protein_pairwise_differences.csv", index=False)

        for comparison_name, config in REFERENCE_COMPARISONS.items():
            ref_id = config["reference_id"]
            country = config["country"]
            if ref_id not in usable:
                continue
            country_ids = set(
                metadata[
                    metadata["country"].eq(country)
                    & metadata["tip_id"].isin(usable.keys())
                ]["tip_id"]
            )
            country_subs = substitution_table(gene, usable, metadata, ref_id, country_ids)
            country_subs.to_csv(
                OUTDIR / f"{gene}_{comparison_name}_aa_substitutions.csv",
                index=False,
            )
            country_summary = country_specific_summary(country_subs, len(country_ids))
            country_summary.to_csv(
                OUTDIR / f"{gene}_{comparison_name}_substitution_summary.csv",
                index=False,
            )

        qc_df = pd.DataFrame(qc_rows)
        gene_qc_rows.append(
            {
                "gene": gene,
                "reference_coordinates": f"{start}-{end}",
                "expected_cds_length_nt": end - start + 1,
                "records_extracted": len(qc_df),
                "usable_proteins": int(qc_df["usable_for_substitution_analysis"].sum()),
                "records_with_internal_stop": int((qc_df["internal_stop_count"] > 0).sum()),
                "records_with_ambiguous_nt": int((qc_df["ambiguous_nt_count"] > 0).sum()),
            }
        )

    all_subs_df = pd.concat(all_subs, ignore_index=True)
    all_subs_df.to_csv(OUTDIR / "FG_all_aa_substitutions_vs_AF212302_MalaysiaPrototypeOnly.csv", index=False)

    pairwise_refs = []
    for gene in GENE_COORDS:
        p = OUTDIR / f"{gene}_reference_protein_pairwise_differences.csv"
        if p.exists():
            pairwise_refs.append(pd.read_csv(p))
    if pairwise_refs:
        pd.concat(pairwise_refs, ignore_index=True).to_csv(
            OUTDIR / "FG_reference_protein_pairwise_differences.csv",
            index=False,
        )

    lineage_summary, country_summary = summarize_substitutions(all_subs_df, metadata, usable_by_gene)
    lineage_summary.to_csv(OUTDIR / "FG_substitution_summary_by_lineage.csv", index=False)
    country_summary.to_csv(OUTDIR / "FG_substitution_summary_by_country.csv", index=False)
    enriched = lineage_enriched(lineage_summary)
    enriched.to_csv(OUTDIR / "FG_lineage_enriched_substitutions.csv", index=False)
    pd.DataFrame(gene_qc_rows).to_csv(OUTDIR / "FG_extraction_translation_qc_summary.csv", index=False)

    notes = [
        "F/G protein analysis",
        "",
        f"Prototype-relative reference sequence: {REFERENCE_ID}",
        "This prototype-relative comparison is not the only reference framework.",
        "Country-specific outputs are also generated using AY988601.1 for Bangladesh, MH396625.1 for India, and AF212302.2 for Malaysia.",
        "Hendra outgroup NC_001906.3 was excluded from Nipah F/G substitution calls.",
        "",
        "Outputs:",
        "- F_proteins_usable.fasta",
        "- G_proteins_usable.fasta",
        "- F_aa_substitutions_vs_AF212302.csv",
        "- G_aa_substitutions_vs_AF212302.csv",
        "- FG_all_aa_substitutions_vs_AF212302_MalaysiaPrototypeOnly.csv",
        "- F/G country-specific *_aa_substitutions.csv",
        "- F/G country-specific *_substitution_summary.csv",
        "- FG_reference_protein_pairwise_differences.csv",
        "- FG_substitution_summary_by_lineage.csv",
        "- FG_substitution_summary_by_country.csv",
        "- FG_lineage_enriched_substitutions.csv",
        "- FG_extraction_translation_qc_summary.csv",
    ]
    (OUTDIR / "README_FG_PROTEIN_ANALYSIS.txt").write_text("\n".join(notes) + "\n", encoding="utf-8")

    print("Wrote F/G protein analysis outputs to", OUTDIR)
    print(pd.DataFrame(gene_qc_rows).to_string(index=False))
    if not enriched.empty:
        print("\nTop enriched substitutions:")
        print(enriched.sort_values(["gene", "aa_position"]).head(30).to_string(index=False))


if __name__ == "__main__":
    main()
