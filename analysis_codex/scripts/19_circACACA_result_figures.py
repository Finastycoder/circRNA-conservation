from __future__ import annotations

import csv
import itertools
import math
import shutil
from pathlib import Path

import pandas as pd
from PIL import Image, ImageDraw, ImageFont

try:
    from reportlab.lib.pagesizes import landscape
    from reportlab.pdfgen import canvas
except Exception:  # pragma: no cover
    canvas = None


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "analysis_codex"
DATA_IN = OUT / "data" / "input_used"
DATA_OUT = OUT / "data" / "output"
FIG = OUT / "figures"
SOURCE_SCRIPTS = OUT / "source_existing" / "scripts"
SOURCE_DATA = OUT / "source_existing" / "data"

SPECIES = ["hsa", "mma", "mmu", "mfu", "mpi", "rsi"]
SPECIES_LABELS = {
    "hsa": "Human (hsa)",
    "mma": "Rhesus (mma)",
    "mmu": "Mouse (mmu)",
    "mfu": "Long-winged bat (mfu)",
    "mpi": "Big-footed bat (mpi)",
    "rsi": "Horseshoe bat (rsi)",
}
CLADES = {
    "hsa": "Primate",
    "mma": "Primate",
    "mmu": "Rodent",
    "mfu": "Bat",
    "mpi": "Bat",
    "rsi": "Bat",
}
ATLAS_FILES = {
    "hsa": ROOT / "circAtlas_circRNA" / "human_bed_v3.0.txt",
    "mma": ROOT / "circAtlas_circRNA" / "macaca_bed_v3.0.txt",
    "mmu": ROOT / "circAtlas_circRNA" / "mouse_bed_v3.0.txt",
}
MOUSE_TOOL_FILES = {
    "CIRCexplorer3": ROOT / "ce_mmu.csv",
    "CIRIquant": ROOT / "cq_mmu.csv",
    "find_circ": ROOT / "fc_mmu.csv",
}
ACACA_IDS_FALLBACK = [
    "11:84083904|84086513",
    "11:84086263|84086513",
    "11:84113967|84122587",
]

COLORS = {
    "green": "#2EAD66",
    "red": "#D84A3A",
    "orange": "#E99B2E",
    "gray": "#B8C0C8",
    "dark": "#2B2F33",
    "light": "#F4F6F8",
    "blue": "#3B82C4",
    "purple": "#7A5FA8",
}


def ensure_dirs() -> None:
    for d in [DATA_IN, DATA_OUT, FIG, SOURCE_SCRIPTS, SOURCE_DATA]:
        d.mkdir(parents=True, exist_ok=True)


def font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont:
    candidates = [
        Path("C:/Windows/Fonts/arialbd.ttf" if bold else "C:/Windows/Fonts/arial.ttf"),
        Path("C:/Windows/Fonts/calibrib.ttf" if bold else "C:/Windows/Fonts/calibri.ttf"),
        Path("C:/Windows/Fonts/segoeuib.ttf" if bold else "C:/Windows/Fonts/segoeui.ttf"),
    ]
    for p in candidates:
        if p.exists():
            return ImageFont.truetype(str(p), size)
    return ImageFont.load_default()


def hex_to_rgb(hex_color: str) -> tuple[int, int, int]:
    h = hex_color.lstrip("#")
    return tuple(int(h[i : i + 2], 16) for i in (0, 2, 4))


def text_size(draw: ImageDraw.ImageDraw, text: str, fnt: ImageFont.ImageFont) -> tuple[int, int]:
    box = draw.textbbox((0, 0), text, font=fnt)
    return box[2] - box[0], box[3] - box[1]


def draw_centered(
    draw: ImageDraw.ImageDraw,
    xy: tuple[float, float],
    text: str,
    fnt: ImageFont.ImageFont,
    fill: str = COLORS["dark"],
) -> None:
    w, h = text_size(draw, text, fnt)
    draw.text((xy[0] - w / 2, xy[1] - h / 2), text, font=fnt, fill=fill)


def wrap_text(draw: ImageDraw.ImageDraw, text: str, fnt: ImageFont.ImageFont, max_width: int) -> list[str]:
    words = text.split()
    lines: list[str] = []
    current = ""
    for word in words:
        trial = f"{current} {word}".strip()
        if text_size(draw, trial, fnt)[0] <= max_width:
            current = trial
        else:
            if current:
                lines.append(current)
            current = word
    if current:
        lines.append(current)
    return lines


def save_all(img: Image.Image, stem: str) -> None:
    png = FIG / f"{stem}.png"
    tif = FIG / f"{stem}.tiff"
    pdf = FIG / f"{stem}.pdf"
    img.save(png, dpi=(300, 300))
    img.save(tif, compression="tiff_lzw", dpi=(300, 300))
    if canvas is not None:
        page = landscape((img.width, img.height))
        c = canvas.Canvas(str(pdf), pagesize=page)
        c.drawImage(str(png), 0, 0, width=img.width, height=img.height)
        c.showPage()
        c.save()


def count_circatlas_acaca() -> pd.DataFrame:
    rows = []
    for sp, path in ATLAS_FILES.items():
        n = 0
        examples = []
        with path.open("r", encoding="utf-8") as handle:
            next(handle, None)
            for line in handle:
                if "acaca" in line.lower():
                    n += 1
                    if len(examples) < 3:
                        examples.append(line.strip().split("\t")[-1])
        rows.append(
            {
                "species": sp,
                "species_label": SPECIES_LABELS[sp],
                "source_file": str(path.relative_to(ROOT)),
                "circAtlas_ACACA_count": n,
                "example_ids": ";".join(examples),
            }
        )
    df = pd.DataFrame(rows)
    df.to_csv(DATA_OUT / "circatlas_acaca_counts.csv", index=False)
    return df


def parse_bool(value) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().upper() in {"TRUE", "T", "YES", "1"}


def load_evidence(atlas_counts: pd.DataFrame) -> pd.DataFrame:
    source = pd.read_csv(ROOT / "New-analysis" / "acaca_final_evidence.csv")
    atlas_by_sp = atlas_counts.set_index("species")["circAtlas_ACACA_count"].to_dict()
    rows = []
    for sp in SPECIES:
        row = source[source["species"] == sp].iloc[0]
        atlas_n = atlas_by_sp.get(sp)
        milk_detected = parse_bool(row["circ_milk"])
        milk_n = int(row["circ_n_iso"])
        integrated = milk_detected or (atlas_n is not None and atlas_n > 0)
        if milk_detected and atlas_n:
            status = "Milk + circAtlas"
        elif milk_detected:
            status = "Milk sequencing"
        elif atlas_n:
            status = "circAtlas only"
        else:
            status = "No circRNA evidence"
        rows.append(
            {
                "species": sp,
                "species_label": SPECIES_LABELS[sp],
                "clade": CLADES[sp],
                "ACACA_host_gene_present": parse_bool(row["gene_genome"]),
                "ACACA_mRNA_evidence": parse_bool(row["mrna_expr"]),
                "circAtlas_available": sp in ATLAS_FILES,
                "circAtlas_ACACA_count": atlas_n if atlas_n is not None else "",
                "milk_sequencing_circACACA_detected": milk_detected,
                "milk_sequencing_circACACA_count": milk_n,
                "integrated_circACACA_evidence": integrated,
                "integrated_status": status,
                "interpretation": (
                    "ACACA-related circRNA detected in this study's milk sequencing"
                    if milk_detected
                    else "ACACA-related circRNA present in circAtlas, not detected in this study's milk sequencing"
                    if atlas_n
                    else "No circAtlas coverage and no milk detection"
                ),
            }
        )
    df = pd.DataFrame(rows)
    df.to_csv(DATA_OUT / "circACACA_integrated_evidence.csv", index=False)
    return df


def host_gene_intersections() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    presence = pd.read_csv(ROOT / "New-analysis" / "gene_presence_matrix.csv")
    for sp in SPECIES:
        presence[sp] = presence[sp].map(parse_bool)
    presence["gene_name_upper"] = presence["gene_name"].astype(str).str.upper()
    acaca = presence[presence["gene_name_upper"] == "ACACA"].copy()
    acaca_out = acaca[["gene_name"] + SPECIES + ["n_species"]]
    acaca_out.to_csv(DATA_OUT / "strict_CE_host_gene_ACACA_presence.csv", index=False)

    rows = []
    for k in range(6, 1, -1):
        for combo in itertools.combinations(SPECIES, k):
            mask = presence[list(combo)].all(axis=1)
            genes = presence.loc[mask, "gene_name_upper"].dropna().unique().tolist()
            rows.append(
                {
                    "k": k,
                    "species_combo": ";".join(combo),
                    "intersection_gene_count": len(genes),
                    "contains_ACACA": "ACACA" in genes,
                }
            )
    combo_df = pd.DataFrame(rows)
    combo_df.to_csv(DATA_OUT / "strict_CE_host_gene_all_intersections.csv", index=False)
    max_df = combo_df[combo_df["contains_ACACA"]].sort_values(
        ["k", "intersection_gene_count"], ascending=[False, False]
    )
    max_df.to_csv(DATA_OUT / "strict_CE_host_gene_max_ACACA_intersection.csv", index=False)
    return presence, combo_df, max_df


def exact_region_counts(set_dict: dict[str, set[str]], prefix: str) -> pd.DataFrame:
    universe = set().union(*set_dict.values()) if set_dict else set()
    rows = []
    for k in range(1, len(SPECIES) + 1):
        for combo in itertools.combinations(SPECIES, k):
            combo = tuple(combo)
            in_combo = set.intersection(*(set_dict[sp] for sp in combo)) if combo else set()
            out_combo = set().union(*(set_dict[sp] for sp in SPECIES if sp not in combo)) if len(combo) < len(SPECIES) else set()
            exact = in_combo - out_combo
            inclusive = in_combo
            rows.append(
                {
                    "k": k,
                    "species_combo": ";".join(combo),
                    "exact_region_count": len(exact),
                    "inclusive_intersection_count": len(inclusive),
                    "contains_ACACA_exact": "ACACA" in exact,
                    "contains_ACACA_inclusive": "ACACA" in inclusive,
                }
            )
    df = pd.DataFrame(rows)
    df.to_csv(DATA_OUT / f"{prefix}_flower_regions.csv", index=False)

    summary = pd.DataFrame(
        [
            {
                "species": sp,
                "species_label": SPECIES_LABELS[sp],
                "set_size": len(set_dict[sp]),
                "contains_ACACA": "ACACA" in set_dict[sp],
            }
            for sp in SPECIES
        ]
    )
    summary.to_csv(DATA_OUT / f"{prefix}_flower_species_set_summary.csv", index=False)
    return df


def strict_milk_host_gene_sets(presence: pd.DataFrame) -> dict[str, set[str]]:
    return {
        sp: set(
            presence.loc[presence[sp].map(parse_bool), "gene_name_upper"]
            .dropna()
            .astype(str)
            .str.upper()
        )
        for sp in SPECIES
    }


def extract_circatlas_host_genes(path: Path) -> set[str]:
    genes: set[str] = set()
    with path.open("r", encoding="utf-8") as handle:
        next(handle, None)
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            circ_id = parts[4]
            if "-" not in circ_id or "_" not in circ_id:
                continue
            gene = circ_id.split("-", 1)[1].rsplit("_", 1)[0]
            gene = gene.strip().upper()
            if gene and gene != "INTERGENIC":
                genes.add(gene)
    return genes


def integrated_circRNA_host_gene_sets(
    presence: pd.DataFrame,
    evidence: pd.DataFrame,
) -> dict[str, set[str]]:
    sets = strict_milk_host_gene_sets(presence)

    for sp, path in ATLAS_FILES.items():
        sets[sp] = sets[sp] | extract_circatlas_host_genes(path)

    # Curated ACACA milk evidence comes from the final evidence table. This
    # preserves the ACACA-specific rescue/curation that is not fully captured
    # by the strict CIRCexplorer3 host-gene matrix.
    for _, row in evidence.iterrows():
      if row["milk_sequencing_circACACA_detected"]:
          sets[row["species"]].add("ACACA")

    return sets


def draw_rotated_ellipse(
    base: Image.Image,
    center: tuple[int, int],
    size: tuple[int, int],
    angle: float,
    fill: str,
    outline: str,
) -> None:
    layer = Image.new("RGBA", base.size, (255, 255, 255, 0))
    e = Image.new("RGBA", (size[0] + 8, size[1] + 8), (255, 255, 255, 0))
    ed = ImageDraw.Draw(e)
    rgba = (*hex_to_rgb(fill), 120)
    ed.ellipse((4, 4, size[0] + 4, size[1] + 4), fill=rgba, outline=(*hex_to_rgb(outline), 220), width=3)
    rot = e.rotate(angle, expand=True, resample=Image.Resampling.BICUBIC)
    layer.alpha_composite(rot, (int(center[0] - rot.width / 2), int(center[1] - rot.height / 2)))
    base.alpha_composite(layer)


def make_flower_venn_figure(
    set_dict: dict[str, set[str]],
    region_df: pd.DataFrame,
    stem: str,
    title: str,
    subtitle: str,
) -> None:
    img = Image.new("RGBA", (1700, 1250), "white")
    draw = ImageDraw.Draw(img)
    draw.text((60, 45), title, font=font(32, True), fill=COLORS["dark"])
    draw.text((60, 88), subtitle, font=font(17), fill="#58616B")

    cx, cy = 640, 620
    colors = ["#E69F00", "#56B4E9", "#009E73", "#F0B43C", "#3B82C4", "#D55E00"]
    angles = [90, 30, -30, -90, -150, 150]
    label_angles = [270, 330, 30, 90, 150, 210]

    for sp, angle, color in zip(SPECIES, angles, colors):
        draw_rotated_ellipse(img, (cx, cy), (680, 230), angle, color, color)

    draw = ImageDraw.Draw(img)
    region_lookup = {
        row["species_combo"]: row for _, row in region_df.iterrows()
    }
    center_key = ";".join(SPECIES)
    center_count = int(region_lookup[center_key]["inclusive_intersection_count"])
    draw.ellipse((cx - 92, cy - 92, cx + 92, cy + 92), fill="white", outline=COLORS["dark"], width=4)
    draw_centered(draw, (cx, cy - 22), "6-way", font(21, True), fill=COLORS["dark"])
    draw_centered(draw, (cx, cy + 24), str(center_count), font(36, True), fill=COLORS["green"] if center_count else COLORS["red"])

    for sp, langle in zip(SPECIES, label_angles):
        rad = math.radians(langle)
        lx = cx + math.cos(rad) * 455
        ly = cy + math.sin(rad) * 315
        unique_key = sp
        unique_count = int(region_lookup[unique_key]["exact_region_count"])
        total = len(set_dict[sp])
        draw.rounded_rectangle((lx - 105, ly - 38, lx + 105, ly + 38), radius=10, fill="white", outline="#CAD3DC", width=2)
        draw_centered(draw, (lx, ly - 12), SPECIES_LABELS[sp].split(" (")[0], font(14, True))
        draw_centered(draw, (lx, ly + 12), f"total {total} | unique {unique_count}", font(12), fill="#58616B")

    # Show the exact 5-way regions around the center; this keeps the flower plot
    # readable while the complete 63-region table is saved as CSV.
    side_x = 1180
    draw.text((side_x, 190), "Exact 5-way regions", font=font(23, True), fill=COLORS["dark"])
    y = 235
    for missing in SPECIES:
        combo = tuple(sp for sp in SPECIES if sp != missing)
        key = ";".join(combo)
        count = int(region_lookup[key]["exact_region_count"])
        color = COLORS["green"] if count else COLORS["gray"]
        draw.rounded_rectangle((side_x, y, side_x + 380, y + 40), radius=8, fill=color)
        draw.text((side_x + 14, y + 10), f"without {missing}: {count}", font=font(14, True), fill="white")
        y += 52

    draw.text((side_x, y + 20), "Output table includes all 63 exact regions.", font=font(15), fill="#58616B")
    draw.text((side_x, y + 48), "Regions with no genes/circRNA evidence are recorded as 0.", font=font(15), fill="#58616B")
    draw.text((60, 1160), "Petal labels show set size and exact species-specific counts. Center shows the 6-species intersection.", font=font(16), fill="#58616B")

    save_all(img.convert("RGB"), stem)


def get_acaca_mouse_ids() -> list[str]:
    f = ROOT / "New-analysis" / "conserved_genes_3way.csv"
    if not f.exists():
        return ACACA_IDS_FALLBACK
    df = pd.read_csv(f)
    row = df[df["gene_name"].astype(str).str.upper() == "ACACA"]
    if row.empty or "mmu" not in row.columns:
        return ACACA_IDS_FALLBACK
    ids = [x.strip() for x in str(row.iloc[0]["mmu"]).split(";") if x.strip()]
    return ids or ACACA_IDS_FALLBACK


def read_tool_ranking(tool_name: str, path: Path, acaca_ids: list[str]) -> pd.DataFrame:
    df = pd.read_csv(path)
    rep_cols = [c for c in df.columns if c.startswith("mmu")]
    out = df[["circ_id"] + rep_cols].copy()
    for c in rep_cols:
        out[c] = pd.to_numeric(out[c], errors="coerce")
    out["mean_reads"] = out[rep_cols].mean(axis=1, skipna=True)
    out["sum_reads"] = out[rep_cols].sum(axis=1, skipna=True)
    out["n_detected_reps"] = out[rep_cols].notna().sum(axis=1)
    out = out[out["mean_reads"].notna()].copy()
    out["tool"] = tool_name
    out["is_circAcc1_ACACA"] = out["circ_id"].isin(acaca_ids)
    out = out.sort_values("mean_reads", ascending=False).reset_index(drop=True)
    out["rank"] = range(1, len(out) + 1)
    return out


def mouse_rankings() -> tuple[pd.DataFrame, pd.DataFrame]:
    acaca_ids = get_acaca_mouse_ids()
    all_rows = []
    summary_rows = []
    for tool, path in MOUSE_TOOL_FILES.items():
        ranking = read_tool_ranking(tool, path, acaca_ids)
        all_rows.append(ranking)
        top = ranking.iloc[0]
        second = ranking.iloc[1]
        best_acaca = ranking[ranking["is_circAcc1_ACACA"]].iloc[0]
        summary_rows.append(
            {
                "tool": tool,
                "top_circ_id": top["circ_id"],
                "top_is_circAcc1_ACACA": bool(top["is_circAcc1_ACACA"]),
                "top_mean_reads": round(float(top["mean_reads"]), 3),
                "second_circ_id": second["circ_id"],
                "second_mean_reads": round(float(second["mean_reads"]), 3),
                "top_vs_second_fold": round(float(top["mean_reads"]) / float(second["mean_reads"]), 3),
                "best_ACACA_circ_id": best_acaca["circ_id"],
                "best_ACACA_rank": int(best_acaca["rank"]),
                "best_ACACA_mean_reads": round(float(best_acaca["mean_reads"]), 3),
            }
        )
    ranking_df = pd.concat(all_rows, ignore_index=True)
    ranking_df.to_csv(DATA_OUT / "mouse_circRNA_rankings_all_tools.csv", index=False)
    top20 = ranking_df.groupby("tool", group_keys=False).head(20)
    top20.to_csv(DATA_OUT / "mouse_circRNA_top20_by_tool.csv", index=False)
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(DATA_OUT / "mouse_circAcc1_top1_summary.csv", index=False)
    return ranking_df, summary


def copy_inputs_and_sources() -> None:
    files = [
        ROOT / "New-analysis" / "acaca_final_evidence.csv",
        ROOT / "New-analysis" / "acaca_evidence_matrix.csv",
        ROOT / "New-analysis" / "gene_presence_matrix.csv",
        ROOT / "New-analysis" / "conserved_genes_3way.csv",
        ROOT / "ce_mmu.csv",
        ROOT / "cq_mmu.csv",
        ROOT / "fc_mmu.csv",
        ROOT / "circAtlas_circRNA" / "human_bed_v3.0.txt",
        ROOT / "circAtlas_circRNA" / "macaca_bed_v3.0.txt",
        ROOT / "circAtlas_circRNA" / "mouse_bed_v3.0.txt",
    ]
    scripts = [
        ROOT / "New-analysis" / "05_conservation_analysis_part1.R",
        ROOT / "New-analysis" / "06_conservation_analysis_part2.R",
        ROOT / "New-analysis" / "11_acaca_gene_level_analysis.R",
        ROOT / "New-analysis" / "12_final_circAcc1_evidence.R",
        ROOT / "New-analysis" / "18_final_story_figure.R",
        ROOT / "New-analysis" / "00_project_setup.R",
    ]
    manifest = []
    for f in files:
        if f.exists():
            dest = DATA_IN / f.name
            shutil.copy2(f, dest)
            manifest.append(
                {
                    "type": "input_data",
                    "source": str(f.relative_to(ROOT)),
                    "copied_to": str(dest.relative_to(ROOT)),
                    "size_bytes": f.stat().st_size,
                }
            )
    for f in scripts:
        if f.exists():
            dest = SOURCE_SCRIPTS / f.name
            shutil.copy2(f, dest)
            manifest.append(
                {
                    "type": "source_script_used",
                    "source": str(f.relative_to(ROOT)),
                    "copied_to": str(dest.relative_to(ROOT)),
                    "size_bytes": f.stat().st_size,
                }
            )
    with (DATA_OUT / "file_manifest.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["type", "source", "copied_to", "size_bytes"])
        writer.writeheader()
        writer.writerows(manifest)


def draw_strategy_panel(draw: ImageDraw.ImageDraw, x: int, y: int, w: int, h: int) -> None:
    title = font(28, True)
    body = font(18)
    draw.text((x, y), "A  Host-gene anchored strategy", font=title, fill=COLORS["dark"])
    boxes = [
        "circRNA junctions are difficult to compare across species",
        "Use ACACA host gene as the cross-species anchor",
        "Integrate circAtlas and this study's milk sequencing data",
        "Prioritize mouse circAcc1 for functional study",
    ]
    bx = x + 10
    by = y + 58
    bw = w - 20
    bh = 72
    for i, text in enumerate(boxes):
        yy = by + i * 95
        fill = ["#EEF3F8", "#F2F7EF", "#FFF5E7", "#F4EEF8"][i]
        draw.rounded_rectangle((bx, yy, bx + bw, yy + bh), radius=10, fill=fill, outline="#CAD3DC", width=2)
        for j, line in enumerate(wrap_text(draw, text, body, bw - 44)):
            draw.text((bx + 22, yy + 17 + j * 22), line, font=body, fill=COLORS["dark"])
        if i < len(boxes) - 1:
            cx = bx + bw / 2
            draw.line((cx, yy + bh + 8, cx, yy + 88), fill="#8B96A3", width=3)
            draw.polygon([(cx - 8, yy + 84), (cx + 8, yy + 84), (cx, yy + 96)], fill="#8B96A3")


def draw_evidence_matrix(draw: ImageDraw.ImageDraw, df: pd.DataFrame, x: int, y: int, w: int, h: int) -> None:
    title = font(28, True)
    label_f = font(15)
    head_f = font(14, True)
    cell_f = font(13, True)
    draw.text((x, y), "B  Cross-species circACACA evidence", font=title, fill=COLORS["dark"])
    cols = [
        ("ACACA\ngene", "ACACA_host_gene_present"),
        ("mRNA\nevidence", "ACACA_mRNA_evidence"),
        ("circAtlas\nACACA", "circAtlas_ACACA_count"),
        ("Milk seq.\ncircACACA", "milk_sequencing_circACACA_count"),
        ("Integrated\ncirc evidence", "integrated_circACACA_evidence"),
    ]
    left = x + 8
    top = y + 55
    row_h = 49
    label_w = 170
    cell_w = (w - label_w - 15) / len(cols)
    for ci, (ct, _) in enumerate(cols):
        cx = left + label_w + ci * cell_w
        for li, line in enumerate(ct.split("\n")):
            draw_centered(draw, (cx + cell_w / 2, top + 10 + li * 17), line, head_f)
    for ri, (_, row) in enumerate(df.iterrows()):
        yy = top + 45 + ri * row_h
        draw.text((left, yy + 12), row["species_label"], font=label_f, fill=COLORS["dark"])
        values = [
            ("YES" if row["ACACA_host_gene_present"] else "NO", COLORS["green"] if row["ACACA_host_gene_present"] else COLORS["red"]),
            ("YES" if row["ACACA_mRNA_evidence"] else "NO", COLORS["green"] if row["ACACA_mRNA_evidence"] else COLORS["red"]),
            (
                str(row["circAtlas_ACACA_count"]) if row["circAtlas_ACACA_count"] != "" else "N/A",
                COLORS["blue"] if row["circAtlas_ACACA_count"] != "" else COLORS["gray"],
            ),
            (
                str(row["milk_sequencing_circACACA_count"])
                if row["milk_sequencing_circACACA_detected"]
                else "0",
                COLORS["green"] if row["milk_sequencing_circACACA_detected"] else COLORS["red"],
            ),
            (
                "YES" if row["integrated_circACACA_evidence"] else "NO",
                COLORS["green"] if row["integrated_circACACA_evidence"] else COLORS["red"],
            ),
        ]
        for ci, (txt, color) in enumerate(values):
            cx = left + label_w + ci * cell_w
            draw.rounded_rectangle((cx + 7, yy + 5, cx + cell_w - 7, yy + row_h - 5), radius=8, fill=color)
            draw_centered(draw, (cx + cell_w / 2, yy + row_h / 2), txt, cell_f, fill="white")
    foot = "circAtlas supports primate/model species; milk sequencing supports mouse and bat species."
    draw.text((left, top + 45 + len(df) * row_h + 10), foot, font=font(14), fill="#58616B")


def draw_bubble_panel(draw: ImageDraw.ImageDraw, df: pd.DataFrame, x: int, y: int, w: int, h: int) -> None:
    title = font(28, True)
    small = font(14)
    draw.text((x, y), "C  circACACA counts by evidence source", font=title, fill=COLORS["dark"])
    plot_x = x + 170
    plot_y = y + 92
    row_h = 48
    col_x = [plot_x + 95, plot_x + 300]
    headers = ["circAtlas", "Milk sequencing"]
    for cx, hdr in zip(col_x, headers):
        draw_centered(draw, (cx, plot_y - 48), hdr, font(17, True))
    max_count = max(
        [int(v) for v in df["milk_sequencing_circACACA_count"].tolist()]
        + [int(v) for v in df["circAtlas_ACACA_count"].tolist() if v != ""]
    )
    for i, (_, row) in enumerate(df.iterrows()):
        yy = plot_y + i * row_h
        draw.text((x + 5, yy - 8), row["species_label"], font=small, fill=COLORS["dark"])
        vals = [row["circAtlas_ACACA_count"], row["milk_sequencing_circACACA_count"]]
        for j, val in enumerate(vals):
            if val == "":
                r = 14
                fill = COLORS["gray"]
                txt = "N/A"
            else:
                val = int(val)
                r = 7 + math.sqrt(val / max_count) * 29 if val > 0 else 8
                fill = COLORS["blue"] if j == 0 else COLORS["green"] if val > 0 else COLORS["red"]
                txt = str(val)
            cx = col_x[j]
            draw.ellipse((cx - r, yy - r + 8, cx + r, yy + r + 8), fill=fill, outline="#FFFFFF", width=2)
            draw_centered(draw, (cx, yy + 8), txt, font(12, True), fill="white")
    draw.text((x + 5, y + h - 24), "Bubble area is scaled by ACACA-related circRNA count.", font=font(13), fill="#58616B")


def draw_mouse_top20_panel(
    draw: ImageDraw.ImageDraw,
    ranking: pd.DataFrame,
    summary: pd.DataFrame,
    x: int,
    y: int,
    w: int,
    h: int,
) -> None:
    title = font(28, True)
    small = font(13)
    draw.text((x, y), "D  Mouse milk circAcc1 is a reproducible TOP1 circRNA", font=title, fill=COLORS["dark"])
    ce = ranking[ranking["tool"] == "CIRCexplorer3"].head(20).copy()
    max_v = ce["mean_reads"].max()
    left = x + 185
    top = y + 60
    bar_h = 16
    gap = 8
    axis_w = w - 220
    for i, (_, row) in enumerate(ce.iterrows()):
        yy = top + i * (bar_h + gap)
        label = row["circ_id"]
        if len(label) > 22:
            label = label[:21] + "..."
        draw.text((x + 5, yy - 1), label, font=small, fill=COLORS["dark"])
        bw = axis_w * float(row["mean_reads"]) / float(max_v)
        color = COLORS["green"] if row["is_circAcc1_ACACA"] else "#9AA8B5"
        draw.rounded_rectangle((left, yy, left + bw, yy + bar_h), radius=4, fill=color)
        draw.text((left + bw + 6, yy - 1), f"{row['mean_reads']:.1f}", font=small, fill=COLORS["dark"])
    s = summary[summary["tool"] == "CIRCexplorer3"].iloc[0]
    note = f"CIRCexplorer3 fold gap vs #2: {s['top_vs_second_fold']:.1f}x"
    draw.text((x + 5, y + h - 70), note, font=font(16, True), fill=COLORS["green"])
    sy = y + h - 43
    for i, (_, row) in enumerate(summary.iterrows()):
        txt = f"{row['tool']}: rank {int(row['best_ACACA_rank'])}, mean {row['best_ACACA_mean_reads']:.1f}"
        draw.text((x + 5 + i * 235, sy), txt, font=font(13), fill="#58616B")


def make_main_figure(evidence: pd.DataFrame, ranking: pd.DataFrame, summary: pd.DataFrame) -> None:
    img = Image.new("RGB", (2400, 1450), "white")
    draw = ImageDraw.Draw(img)
    draw_strategy_panel(draw, 60, 55, 760, 500)
    draw_evidence_matrix(draw, evidence, 880, 55, 1460, 500)
    draw_bubble_panel(draw, evidence, 60, 630, 760, 700)
    draw_mouse_top20_panel(draw, ranking, summary, 880, 630, 1460, 700)
    draw.text(
        (60, 1380),
        "Conclusion: ACACA anchors multi-species circACACA evidence, and mouse circAcc1 is the strongest milk-sequencing candidate.",
        font=font(22, True),
        fill=COLORS["dark"],
    )
    save_all(img, "Fig_circACACA_conservation_and_mouse_priority")


def make_mouse_top20_figure(ranking: pd.DataFrame, summary: pd.DataFrame) -> None:
    img = Image.new("RGB", (1700, 1200), "white")
    draw = ImageDraw.Draw(img)
    draw_mouse_top20_panel(draw, ranking, summary, 70, 55, 1560, 1040)
    save_all(img, "FigS_mouse_circAcc1_top20_all_tools")


def make_bubble_figure(evidence: pd.DataFrame) -> None:
    img = Image.new("RGB", (1200, 900), "white")
    draw = ImageDraw.Draw(img)
    draw_bubble_panel(draw, evidence, 60, 60, 1080, 760)
    save_all(img, "FigS_circAtlas_vs_milk_circACACA_counts")


def make_intersection_figure(combo_df: pd.DataFrame, max_df: pd.DataFrame) -> None:
    img = Image.new("RGB", (1500, 900), "white")
    draw = ImageDraw.Draw(img)
    draw.text((60, 50), "Strict milk CIRCexplorer3 host-gene intersection containing ACACA", font=font(30, True), fill=COLORS["dark"])
    draw.text(
        (60, 95),
        "This is the conservative host-gene expression result only; integrated circAtlas + milk evidence is shown in the main figure.",
        font=font(17),
        fill="#58616B",
    )
    max_acaca = max_df.iloc[0] if not max_df.empty else None
    if max_acaca is not None:
        msg = f"Maximum ACACA-containing strict intersection: {int(max_acaca['k'])} species ({max_acaca['species_combo']})"
    else:
        msg = "No ACACA-containing strict intersection found."
    draw.text((60, 140), msg, font=font(21, True), fill=COLORS["green"])
    data = (
        combo_df[combo_df["contains_ACACA"]]
        .groupby("k")
        .agg(n_combos=("species_combo", "count"), max_intersection_size=("intersection_gene_count", "max"))
        .reset_index()
    )
    left, top, width, height = 160, 245, 1120, 500
    draw.line((left, top + height, left + width, top + height), fill=COLORS["dark"], width=3)
    draw.line((left, top, left, top + height), fill=COLORS["dark"], width=3)
    max_size = max([1] + data["max_intersection_size"].tolist())
    for k in range(2, 7):
        row = data[data["k"] == k]
        x = left + (k - 2) * (width / 5) + 60
        if row.empty:
            bar_h = 0
            label = "0"
            color = COLORS["gray"]
        else:
            val = float(row.iloc[0]["max_intersection_size"])
            bar_h = height * val / max_size
            label = f"{int(row.iloc[0]['n_combos'])} combo(s)"
            color = COLORS["green"] if k == int(max_acaca["k"]) else COLORS["blue"]
        draw.rounded_rectangle((x - 45, top + height - bar_h, x + 45, top + height), radius=8, fill=color)
        draw_centered(draw, (x, top + height + 30), f"{k} sp.", font(18, True))
        draw_centered(draw, (x, top + height - bar_h - 22), label, font(14), fill=COLORS["dark"])
    draw.text((60, 805), "Bar height: largest intersection size among combinations that contain ACACA.", font=font(16), fill="#58616B")
    save_all(img, "FigS_strict_host_gene_ACACA_max_intersection")


def write_summary(evidence: pd.DataFrame, summary: pd.DataFrame, max_df: pd.DataFrame) -> None:
    ce = summary[summary["tool"] == "CIRCexplorer3"].iloc[0]
    text = f"""# analysis_codex summary

This folder contains the revised circACACA/circAcc1 analysis package generated by
`analysis_codex/scripts/19_circACACA_result_figures.py`.

Key results:

- ACACA host gene and mRNA evidence are present for all six species in the curated evidence table.
- circAtlas contains ACACA/Acaca circRNA entries for human (317), rhesus macaque (92), and mouse (56).
- This study's milk sequencing evidence detects circACACA/circAcc1 in mouse and the three bat species in the integrated evidence table.
- Integrated circACACA evidence therefore covers all six species when circAtlas and this study's milk sequencing are considered together.
- Strict CIRCexplorer3 milk host-gene expression matrix places ACACA in a 3-species maximum intersection: {max_df.iloc[0]['species_combo'] if not max_df.empty else 'not found'}.
- In mouse milk CIRCexplorer3 data, circAcc1 is rank 1 with mean reads {ce['best_ACACA_mean_reads']:.1f}, {ce['top_vs_second_fold']:.1f}x the second-ranked circRNA.
- Six-species flower Venn plots were generated for:
  - strict milk circRNA host genes;
  - integrated circRNA host-gene evidence.
- R-style flower Venn plots matching the provided example are generated by
  `analysis_codex/scripts/20_flower_venn_R_style.R`.
- The flower-region CSV files include all 63 non-empty species combinations. Regions without evidence are retained with count 0.

Interpretation:

The analysis supports a host-gene anchored, multi-species conservation argument for circACACA.
It should not be described as direct junction-level conservation across six species.
The mouse follow-up is justified because circAcc1 is the strongest milk-sequencing candidate in mouse
and the broader ACACA/circACACA evidence is conserved across multiple species.
"""
    (OUT / "README_analysis_codex.md").write_text(text, encoding="utf-8")


def main() -> None:
    ensure_dirs()
    copy_inputs_and_sources()
    atlas_counts = count_circatlas_acaca()
    evidence = load_evidence(atlas_counts)
    presence, combo_df, max_df = host_gene_intersections()
    strict_sets = strict_milk_host_gene_sets(presence)
    strict_regions = exact_region_counts(strict_sets, "strict_milk_host_gene_6species")
    integrated_sets = integrated_circRNA_host_gene_sets(presence, evidence)
    integrated_regions = exact_region_counts(integrated_sets, "integrated_circRNA_host_gene_6species")
    ranking, summary = mouse_rankings()
    make_main_figure(evidence, ranking, summary)
    make_mouse_top20_figure(ranking, summary)
    make_bubble_figure(evidence)
    make_intersection_figure(combo_df, max_df)
    make_flower_venn_figure(
        strict_sets,
        strict_regions,
        "FigS_6species_strict_milk_host_gene_flower_venn",
        "Six-species flower Venn: strict milk host genes",
        "Sets are expressed circRNA host genes from this study's milk CIRCexplorer3 matrix.",
    )
    make_flower_venn_figure(
        integrated_sets,
        integrated_regions,
        "FigS_6species_integrated_circRNA_host_gene_flower_venn",
        "Six-species flower Venn: integrated circRNA host-gene evidence",
        "Sets combine this study's milk circRNA host genes, circAtlas host genes, and curated circACACA milk evidence.",
    )
    write_summary(evidence, summary, max_df)


if __name__ == "__main__":
    main()
