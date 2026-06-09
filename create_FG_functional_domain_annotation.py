from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


ROOT = Path(__file__).resolve().parents[1]
FG_DIR = ROOT / "08_reports_and_manuscript" / "analysis_outputs" / "FG_protein_analysis"
OUT_DIR = ROOT / "11_FG_functional_domain_annotation"
FIG_DIR = OUT_DIR / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)


INPUTS = [
    ("Bangladesh", "F", FG_DIR / "F_Bangladesh_vs_AY988601_substitution_summary.csv"),
    ("Bangladesh", "G", FG_DIR / "G_Bangladesh_vs_AY988601_substitution_summary.csv"),
    ("India", "F", FG_DIR / "F_India_vs_MH396625_substitution_summary.csv"),
    ("India", "G", FG_DIR / "G_India_vs_MH396625_substitution_summary.csv"),
    ("Malaysia", "F", FG_DIR / "F_Malaysia_vs_AF212302_substitution_summary.csv"),
    ("Malaysia", "G", FG_DIR / "G_Malaysia_vs_AF212302_substitution_summary.csv"),
]


def f_region(pos: int):
    # Approximate screening regions for NiV F0. Coordinates are used for
    # interpretation triage only, not as experimentally defined boundaries.
    if 1 <= pos <= 26:
        return "Signal peptide / N-terminal region", "Mostly trafficking/processing context; not directly interpretable without experiments."
    if 27 <= pos <= 108:
        return "F2 ectodomain, upstream of cathepsin-L cleavage site", "F2 region; variants here may affect precursor stability or processing, but evidence is required."
    if pos == 109 or pos == 108:
        return "Cleavage-site-adjacent region", "Close to the F0 activation boundary; potentially interesting but requires experimental validation."
    if 110 <= pos <= 136:
        return "Fusion peptide / immediate post-cleavage region", "Functionally important fusion-initiation region."
    if 137 <= pos <= 184:
        return "HR1 / fusion-core region", "Heptad-repeat region involved in six-helix bundle formation during fusion."
    if 185 <= pos <= 452:
        return "F1 ectodomain core", "Structured ectodomain; interpret as structural-context/hypothesis only unless supported by specific literature."
    if 453 <= pos <= 488:
        return "HR2 / membrane-proximal ectodomain", "Heptad-repeat region near membrane; relevant to fusion-core formation."
    if 489 <= pos <= 511:
        return "Transmembrane region", "Membrane anchor; may affect trafficking or fusion regulation, but functional effect is not inferred here."
    if pos >= 512:
        return "Cytoplasmic tail", "C-terminal tail; may affect trafficking/endocytosis, but functional effect is not inferred here."
    return "Unassigned", "Position outside expected annotated screening range."


G_INTERFACE_RESIDUES = {505, 533, 559, 579, 581, 588}


def g_region(pos: int):
    # NiV G is a type-II membrane glycoprotein. Stalk/head boundaries follow
    # published structural summaries of the ectodomain; N-terminal membrane
    # topology is approximate and used only for triage.
    if 1 <= pos <= 44:
        region = "Cytoplasmic tail / N-terminal region"
        note = "Type-II G N terminus; not part of the resolved receptor-binding head."
    elif 45 <= pos <= 70:
        region = "Transmembrane region"
        note = "Membrane anchor; not directly interpretable for receptor binding."
    elif 71 <= pos <= 95:
        region = "Membrane-proximal ectodomain"
        note = "Proximal ectodomain before the published stalk boundary."
    elif 96 <= pos <= 147:
        region = "Stalk domain"
        note = "Stalk region can influence attachment protein architecture and F triggering."
    elif 148 <= pos <= 165:
        region = "Neck domain"
        note = "Neck region connects stalk to the globular head."
    elif 166 <= pos <= 177:
        region = "Linker to globular head"
        note = "Connector region between neck and receptor-binding head."
    elif 178 <= pos <= 602:
        region = "Globular head / receptor-binding domain"
        note = "Within the beta-propeller head that binds ephrin-B2/B3."
    else:
        region = "Unassigned"
        note = "Position outside expected annotated screening range."

    if pos in G_INTERFACE_RESIDUES:
        note += " This residue is reported in or directly adjacent to the ephrin-binding interface/pocket."
    elif 570 <= pos <= 590:
        note += " This position is close to reported ephrin-binding pocket residues; interpret as structural proximity only."
    elif 500 <= pos <= 590:
        note += " This position lies in the broader receptor-binding head region; functional effect is hypothesis-generating only."
    return region, note


def interpretation_level(row):
    sub = row["substitution"]
    gene = row["gene"]
    pos = int(row["aa_position"])
    country = row["country"]
    percent = float(row["percent"])

    if country == "Bangladesh" and gene == "F" and sub in {"S207L", "G252D"}:
        return (
            "known_reference_artifact",
            "Exclude from biological interpretation",
            "Lo et al. 2019 reported AY988601.1 F sequence correction issues at these positions; detected signal is useful QC validation, not a biological mutation.",
        )
    if gene == "G" and (pos in G_INTERFACE_RESIDUES or 570 <= pos <= 590):
        return (
            "structural_proximity",
            "Supplement; mention only cautiously if needed",
            "Located at or near published ephrin-binding pocket/interface positions; receptor-binding effect is not inferred without modeling/experiments.",
        )
    if gene == "F" and (108 <= pos <= 136 or 137 <= pos <= 184 or 453 <= pos <= 488 or 489 <= pos <= 511):
        return (
            "functional_region_proximity",
            "Supplement; discuss only if recurrent",
            "Falls in or near a fusion-relevant region; effect remains hypothesis-generating.",
        )
    if percent >= 20:
        return (
            "recurrent_country_signal",
            "Supplement; main text only if biologically useful",
            "Recurrent within the country-specific reference comparison, but no direct functional evidence was assigned.",
        )
    return (
        "low_frequency_context",
        "Supplement only",
        "Low-frequency or singleton variant; retain for transparency and future surveillance comparison.",
    )


def load_inputs():
    frames = []
    for country, gene, path in INPUTS:
        df = pd.read_csv(path)
        if df.empty:
            continue
        df.insert(0, "country", country)
        frames.append(df)
    return pd.concat(frames, ignore_index=True)


def annotate(df):
    rows = []
    for _, row in df.iterrows():
        pos = int(row["aa_position"])
        if row["gene"] == "F":
            region, region_note = f_region(pos)
        else:
            region, region_note = g_region(pos)
        level, placement, interp_note = interpretation_level(row)
        out = row.to_dict()
        out["domain_or_region"] = region
        out["domain_note"] = region_note
        out["interpretation_level"] = level
        out["recommended_manuscript_placement"] = placement
        out["interpretation_note"] = interp_note
        rows.append(out)
    annotated = pd.DataFrame(rows)
    annotated["is_recurrent_focus"] = (
        (annotated["sequence_count"] >= 2)
        | (annotated["percent"] >= 10)
        | (annotated["interpretation_level"].isin(["known_reference_artifact", "structural_proximity"]))
    )
    return annotated


def make_region_summary(annotated):
    return (
        annotated.groupby(["country", "gene", "domain_or_region"], dropna=False)
        .agg(
            substitutions=("substitution", "count"),
            recurrent_focus_substitutions=("is_recurrent_focus", "sum"),
            max_percent=("percent", "max"),
        )
        .reset_index()
        .sort_values(["gene", "country", "domain_or_region"])
    )


def write_readme(annotated, recurrent):
    text = f"""# F/G Functional Domain Annotation

Generated: 2026-05-08

## Purpose

This analysis maps country-reference-aware F and G amino-acid substitutions onto literature-guided functional regions. It is intended as a biological interpretation layer, not as proof of altered receptor binding, fusion activity, virulence, antigenicity, or host range.

## Reference-Aware Comparisons Used

- Bangladesh: F/G substitutions relative to `AY988601.1`
- India: F/G substitutions relative to `MH396625.1`
- Malaysia: F/G substitutions relative to `AF212302.2`

## Important Interpretation Rule

The Bangladesh F substitutions `S207L` and `G252D` are retained in the annotation table as known reference-artifact/QC-validation signals, but they must not be interpreted as biological mutations. Lo et al. (2019) reported correction issues in the corresponding `AY988601.1` F sequence positions.

## Output Files

- `FG_functional_domain_annotation_all_substitutions.csv`: all country-reference-aware F/G substitutions with domain and interpretation notes.
- `FG_functional_domain_annotation_recurrent_focus.csv`: recurrent, high-frequency, artifact, or structural-proximity variants prioritized for review.
- `FG_functional_domain_region_summary.csv`: compact count summary by country, protein, and region.
- `FG_functional_domain_annotation.xlsx`: same tables in one workbook for manuscript/supplement use.
- `figures/F_G_domain_variant_map.png`: domain schematic showing recurrent-focus substitutions along F and G.

## High-Level Summary

- Total annotated substitutions: {len(annotated)}
- Recurrent-focus substitutions: {len(recurrent)}
- Interpretation levels represented: {", ".join(sorted(annotated["interpretation_level"].unique()))}

## Literature Basis

- NiV F is a class I fusion protein cleaved by cathepsin L into F1/F2 and contains a fusion peptide plus heptad-repeat regions required for membrane fusion.
- NiV G is a type-II attachment glycoprotein with stalk/neck/linker/globular-head organization; the globular head binds ephrin-B2/B3.
- Bowden et al. (2008), Negrete et al. (2005), Wong et al. (2017), and later structural summaries support the G receptor-binding interpretation.
- Lo et al. (2019) supports removal of `AY988601.1` F `S207L` and `G252D` from biological interpretation.

## Safe Manuscript Wording

Protein-level interpretation was limited to consensus amino-acid substitutions mapped onto published functional domains and structural regions. No structural modeling, receptor-binding prediction, fusion assay, or antigenicity prediction was performed; therefore, substitutions are interpreted as surveillance and hypothesis-generating signals unless directly supported by published experimental evidence.
"""
    (OUT_DIR / "README_FG_FUNCTIONAL_DOMAIN_ANNOTATION.md").write_text(text, encoding="utf-8")


LABEL_SUBSTITUTIONS = {
    ("Bangladesh", "F", "M19I"),
    ("Bangladesh", "F", "S207L"),
    ("Bangladesh", "F", "G252D"),
    ("Bangladesh", "F", "A503V"),
    ("India", "F", "Y10N"),
    ("India", "F", "L13F"),
    ("India", "G", "T135A"),
    ("India", "G", "N423D"),
    ("India", "G", "L577P"),
    ("India", "G", "E579Q"),
    ("Bangladesh", "G", "K145R"),
    ("Bangladesh", "G", "K376T"),
    ("Bangladesh", "G", "K571N"),
}

LABEL_ORDER = [
    ("India", "F", "Y10N"),
    ("India", "F", "L13F"),
    ("Bangladesh", "F", "M19I"),
    ("Bangladesh", "F", "S207L"),
    ("Bangladesh", "F", "G252D"),
    ("Bangladesh", "F", "A503V"),
    ("India", "G", "T135A"),
    ("Bangladesh", "G", "K145R"),
    ("Bangladesh", "G", "K376T"),
    ("India", "G", "N423D"),
    ("Bangladesh", "G", "K571N"),
    ("India", "G", "L577P"),
    ("India", "G", "E579Q"),
]

LABEL_NUMBERS = {key: idx + 1 for idx, key in enumerate(LABEL_ORDER)}

LABEL_X_OFFSETS = {
    ("India", "F", "Y10N"): -8,
    ("India", "F", "L13F"): 0,
    ("Bangladesh", "F", "M19I"): 9,
    ("Bangladesh", "G", "K571N"): -13,
    ("India", "G", "L577P"): 0,
    ("India", "G", "E579Q"): 13,
}


def draw_domain(ax, y, length, domains, title, variants, color):
    ax.add_patch(Rectangle((1, y - 0.08), length, 0.16, facecolor="#f3f4f6", edgecolor="#666666", lw=0.8))
    for start, end, label, dcolor in domains:
        ax.add_patch(Rectangle((start, y - 0.12), end - start + 1, 0.24, facecolor=dcolor, edgecolor="white", lw=0.8))
        if (end - start) > 18:
            ax.text((start + end) / 2, y + 0.18, label, ha="center", va="bottom", fontsize=8)
    offsets = {}
    for _, row in variants.iterrows():
        pos = int(row["aa_position"])
        offsets[pos] = offsets.get(pos, 0) + 1
        dy = 0.34 + 0.12 * (offsets[pos] - 1)
        marker = "x" if row["interpretation_level"] == "known_reference_artifact" else "o"
        ax.scatter([pos], [y - dy], s=35 + float(row["percent"]) * 0.8, c=color[row["country"]], marker=marker, edgecolors="black", linewidths=0.5, zorder=5)
        label_key = (row["country"], row["gene"], row["substitution"])
        if label_key in LABEL_SUBSTITUTIONS:
            label = f'{row["country"][0]}:{row["substitution"]}'
            if row["interpretation_level"] == "known_reference_artifact":
                label += "*"
            ax.annotate(
                label,
                xy=(pos, y - dy),
                xytext=(pos, y - dy - 0.22 - 0.05 * (offsets[pos] - 1)),
                ha="center",
                va="top",
                fontsize=8,
                rotation=35,
                arrowprops={"arrowstyle": "-", "lw": 0.4, "color": "#555555"},
            )
    ax.text(1, y + 0.45, title, ha="left", va="center", fontsize=12, fontweight="bold")


def make_variant_map(recurrent):
    colors = {"Bangladesh": "#1f77b4", "India": "#ff7f0e", "Malaysia": "#2ca02c"}
    f_domains = [
        (1, 26, "SP/N", "#9ecae1"),
        (27, 108, "F2", "#c7e9c0"),
        (109, 136, "Cleavage/FP", "#fdd0a2"),
        (137, 184, "HR1", "#fdae6b"),
        (185, 452, "F1 core", "#dadaeb"),
        (453, 488, "HR2", "#bcbddc"),
        (489, 511, "TM", "#969696"),
        (512, 546, "CT", "#636363"),
    ]
    g_domains = [
        (1, 44, "CT/N", "#fbb4ae"),
        (45, 70, "TM", "#b3cde3"),
        (71, 95, "prox.", "#ccebc5"),
        (96, 147, "stalk", "#decbe4"),
        (148, 165, "neck", "#fed9a6"),
        (166, 177, "linker", "#ffffcc"),
        (178, 602, "head/RBD", "#e5d8bd"),
    ]
    fig = plt.figure(figsize=(16, 10))
    gs = fig.add_gridspec(3, 1, height_ratios=[1, 1, 0.72], hspace=0.42)
    ax_f = fig.add_subplot(gs[0, 0])
    ax_g = fig.add_subplot(gs[1, 0], sharex=ax_f)
    ax_key = fig.add_subplot(gs[2, 0])
    ax_key.set_axis_off()

    def setup_axis(ax, title, length, domains, variants):
        ax.set_xlim(1, 610)
        ax.set_ylim(0, 1)
        ax.set_yticks([])
        ax.spines[["left", "right", "top"]].set_visible(False)
        ax.grid(axis="x", color="#dddddd", alpha=0.45)
        ax.text(1, 0.93, title, ha="left", va="top", fontsize=13, fontweight="bold")
        y = 0.58
        ax.add_patch(Rectangle((1, y - 0.055), length, 0.11, facecolor="#f3f4f6", edgecolor="#666666", lw=0.7))
        for start, end, label, dcolor in domains:
            ax.add_patch(Rectangle((start, y - 0.07), end - start + 1, 0.14, facecolor=dcolor, edgecolor="white", lw=0.8))
            if end - start > 15:
                ax.text((start + end) / 2, y + 0.115, label, ha="center", va="bottom", fontsize=9)

        lanes = {}
        for _, row in variants.iterrows():
            pos = int(row["aa_position"])
            lane = lanes.get(pos, 0)
            lanes[pos] = lane + 1
            marker_y = 0.30 - 0.060 * min(lane, 3)
            key = (row["country"], row["gene"], row["substitution"])
            is_artifact = row["interpretation_level"] == "known_reference_artifact"
            marker = "X" if is_artifact else "o"
            size = 42 + min(float(row["percent"]), 100) * 0.75
            ax.scatter(
                [pos],
                [marker_y],
                s=size,
                c=[colors[row["country"]]],
                marker=marker,
                edgecolors="black",
                linewidths=0.8,
                zorder=4,
            )
            if key in LABEL_NUMBERS:
                n = LABEL_NUMBERS[key]
                label_y = marker_y - 0.12 - 0.035 * min(lane, 3)
                label_x = pos + LABEL_X_OFFSETS.get(key, 0)
                ax.plot([pos, label_x], [marker_y - 0.025, label_y + 0.025], color="#777777", lw=0.6, zorder=3)
                ax.text(
                    label_x,
                    label_y,
                    str(n),
                    ha="center",
                    va="center",
                    fontsize=8,
                    fontweight="bold",
                    bbox={"boxstyle": "circle,pad=0.22", "fc": "white", "ec": "#555555", "lw": 0.6},
                    zorder=5,
                )

    setup_axis(ax_f, "F protein", 546, f_domains, recurrent[recurrent["gene"] == "F"])
    setup_axis(ax_g, "G protein", 602, g_domains, recurrent[recurrent["gene"] == "G"])
    ax_g.set_xlabel("Amino-acid position", fontsize=12)
    ax_f.tick_params(labelbottom=False)

    fig.suptitle("F/G substitutions mapped to literature-guided functional regions", fontsize=17, fontweight="bold", y=0.97)

    # Bottom key in two columns.
    ax_key.text(0.02, 0.98, "Numbered substitution key", fontsize=13, fontweight="bold", va="top", transform=ax_key.transAxes)
    y_positions = [0.82, 0.68, 0.54, 0.40, 0.26, 0.12, -0.02]
    x_positions = [0.02, 0.52]
    lookup = recurrent.set_index(["country", "gene", "substitution"], drop=False)
    entries = []
    for key in LABEL_ORDER:
        if key not in lookup.index:
            continue
        row = lookup.loc[key]
        if isinstance(row, pd.DataFrame):
            row = row.iloc[0]
        n = LABEL_NUMBERS[key]
        country, gene, sub = key
        artifact = " *" if row["interpretation_level"] == "known_reference_artifact" else ""
        entries.append((n, country, f"{n}. {country} {gene}:{sub}{artifact}"))
    for i, (_, country, text) in enumerate(entries):
        col = 0 if i < 7 else 1
        row_i = i if i < 7 else i - 7
        x = x_positions[col]
        y = y_positions[row_i]
        ax_key.scatter([x], [y - 0.015], s=65, c=[colors[country]], edgecolors="black", transform=ax_key.transAxes, clip_on=False)
        ax_key.text(x + 0.025, y, text, fontsize=10, va="top", transform=ax_key.transAxes)

    ax_key.text(0.78, 0.82, "Color key", fontsize=11, fontweight="bold", va="top", transform=ax_key.transAxes)
    y = 0.68
    for country, col in colors.items():
        ax_key.scatter([0.78], [y - 0.015], s=70, c=[col], edgecolors="black", transform=ax_key.transAxes, clip_on=False)
        ax_key.text(0.805, y, country, fontsize=10, va="top", transform=ax_key.transAxes)
        y -= 0.14
    ax_key.text(0.78, 0.20, "* Known AY988601.1 F artifact", fontsize=9, va="top", color="#555555", transform=ax_key.transAxes)

    fig.subplots_adjust(left=0.07, right=0.985, top=0.92, bottom=0.08)
    fig.savefig(FIG_DIR / "F_G_domain_variant_map.png", dpi=300)
    fig.savefig(FIG_DIR / "F_G_domain_variant_map.pdf")
    plt.close(fig)


def make_region_focus_barplot(recurrent):
    plot_df = (
        recurrent.groupby(["gene", "domain_or_region"])
        .size()
        .reset_index(name="recurrent_focus_substitutions")
        .sort_values(["gene", "recurrent_focus_substitutions"], ascending=[True, False])
    )
    order = plot_df["domain_or_region"].tolist()
    colors = {"F": "#4c78a8", "G": "#59a14f"}

    fig, ax = plt.subplots(figsize=(11, 5.5))
    labels = [f"{g}: {r}" for g, r in zip(plot_df["gene"], plot_df["domain_or_region"])]
    bars = ax.barh(labels, plot_df["recurrent_focus_substitutions"], color=[colors[g] for g in plot_df["gene"]])
    ax.invert_yaxis()
    ax.set_xlabel("Number of recurrent-focus substitutions")
    ax.set_title("Functional-region distribution of prioritized F/G substitutions", fontsize=14, fontweight="bold")
    ax.grid(axis="x", alpha=0.25)
    for bar, value in zip(bars, plot_df["recurrent_focus_substitutions"]):
        ax.text(value + 0.1, bar.get_y() + bar.get_height() / 2, str(value), va="center", fontsize=10)
    fig.subplots_adjust(left=0.38, right=0.98, top=0.86, bottom=0.12)
    fig.savefig(FIG_DIR / "F_G_recurrent_focus_by_functional_region.png", dpi=300)
    fig.savefig(FIG_DIR / "F_G_recurrent_focus_by_functional_region.pdf")
    plt.close(fig)


def main():
    raw = load_inputs()
    annotated = annotate(raw)
    recurrent = annotated[annotated["is_recurrent_focus"]].copy()
    region_summary = make_region_summary(annotated)

    annotated.to_csv(OUT_DIR / "FG_functional_domain_annotation_all_substitutions.csv", index=False)
    recurrent.to_csv(OUT_DIR / "FG_functional_domain_annotation_recurrent_focus.csv", index=False)
    region_summary.to_csv(OUT_DIR / "FG_functional_domain_region_summary.csv", index=False)

    try:
        with pd.ExcelWriter(OUT_DIR / "FG_functional_domain_annotation.xlsx", engine="openpyxl") as writer:
            recurrent.to_excel(writer, sheet_name="recurrent_focus", index=False)
            annotated.to_excel(writer, sheet_name="all_substitutions", index=False)
            region_summary.to_excel(writer, sheet_name="region_summary", index=False)
    except PermissionError:
        with pd.ExcelWriter(OUT_DIR / "FG_functional_domain_annotation_updated.xlsx", engine="openpyxl") as writer:
            recurrent.to_excel(writer, sheet_name="recurrent_focus", index=False)
            annotated.to_excel(writer, sheet_name="all_substitutions", index=False)
            region_summary.to_excel(writer, sheet_name="region_summary", index=False)

    write_readme(annotated, recurrent)
    make_variant_map(recurrent)
    make_region_focus_barplot(recurrent)
    print(OUT_DIR)


if __name__ == "__main__":
    main()
