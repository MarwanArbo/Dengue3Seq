#!/usr/bin/env python3
"""
Cartographie des amorces sur un genome de reference (+ test sur d'autres genomes).
Version 4 : regroupement F/R explicite (fiable, plus d'heuristique fragile),
diagramme allege et corrige (fini les artefacts de rendu).

Usage :
  python3 primer_mapping_v4.py primers.tsv reference.fasta [autre_genome1.fasta ...]
"""
import sys, re

IUPAC_SET = {
    'A':{'A'},'C':{'C'},'G':{'G'},'T':{'T'},
    'R':{'A','G'},'Y':{'C','T'},'S':{'G','C'},'W':{'A','T'},
    'K':{'G','T'},'M':{'A','C'},'B':{'C','G','T'},'D':{'A','G','T'},
    'H':{'A','C','T'},'V':{'A','C','G'},'N':{'A','C','G','T'}
}
IUPAC_REGEX = {k: ("[" + "".join(v) + "]" if len(v) > 1 else next(iter(v))) for k, v in IUPAC_SET.items()}
COMP = {'A':'T','T':'A','C':'G','G':'C','R':'Y','Y':'R','S':'S','W':'W',
        'K':'M','M':'K','B':'V','V':'B','D':'H','H':'D','N':'N'}

# Regroupement EXPLICITE des amorces par amplicon (pool, [noms des amorces qui le composent]).
# Remplace l'ancienne heuristique par numero, qui laissait DENV3-R-8 orpheline.

# =========================================================================
# CONFIGURATION DU DIAGRAMME — modifie librement les valeurs ci-dessous
# =========================================================================
FIGURE_WIDTH = 8       # largeur de l'image en pouces
FIGURE_HEIGHT = 3.2    # hauteur de l'image en pouces

# Couleur de fond de chaque pool (format hexadecimal #RRGGBB)
POOL_COLORS = {
    "1": "#083248",   # bleu navy
    "2": "#8C0E0F",   # rouge navy
    "3": "#E89C31",   # jaune navy
}

# Couleur des croix signalant les mismatches
MISMATCH_COLOR = "#C0392B"   # rouge

# Couleur de la ligne de fond (backbone) qui traverse chaque rangee
BACKBONE_COLOR = "#E5E4DF"   # gris tres clair

# Couleur du texte a l'INTERIEUR des amplicons (les numeros 1, 2, 3...)
AMPLICON_TEXT_COLOR = "#ffffff"   # blanc

# Couleur du titre du graphique et des labels de pool (Pool 1/2/3 sur l'axe Y)
TITLE_COLOR = "#3A3A38"           # gris fonce

# Couleur des graduations de l'axe des abscisses (0, 2000, 4000, ...)
AXIS_LABEL_COLOR = "#9A9A94"      # gris clair

# Couleur du texte indiquant la taille (pb) de chaque amplicon, sous le rectangle
LENGTH_LABEL_COLOR = "#9A9A94"    # gris clair, italique

# Taille des croix de mismatch (reduite)
MISMATCH_MARKER_SIZE = 6          # etait 13
# =========================================================================

AMPLICON_GROUPS = [
    ("1", "1", ["DENV3-F-1", "CDC-D3-R-1"]),
    ("3", "1", ["CDC-D3-F-3", "CDC-D3-R-3"]),
    ("5", "1", ["CDC-D3-F-5", "DENV3-R1"]),
    ("7", "1", ["CDC-D3-F-7", "CDC-D3-R-7"]),
    ("9", "1", ["CDC-D3-F-9", "CDC-D3-R-9"]),
    ("15", "1", ["CDC-D3-F-15", "CDC-D3-R-15"]),
    ("17", "1", ["CDC-D3-F-17", "CDC-D3-R-17"]),
    ("21", "1", ["CDC-D3-F-21", "CDC-D3-R-21"]),
    ("23", "1", ["CDC-D3-F-23", "CDC-D3-R-23"]),
    ("33", "1", ["CDC-D3-F-33", "CDC-D3-R-33"]),
    ("36", "1", ["CNR-D3-F-36", "CNR-D3-R-36Bis", "DENV3-R-8"]),  # trio, tailles 559/581
    ("2", "2", ["CDC-D3-F-2", "CDC-D3-R-2"]),
    ("4", "2", ["CDC-D3-F-4", "CDC-D3-R-4"]),
    ("8", "2", ["CDC-D3-F-8", "CDC-D3-R-8"]),
    ("10", "2", ["CDC-D3-F-10", "CDC-D3-R-10"]),
    ("12", "2", ["CDC-D3-F-12", "CDC-D3-R-12"]),
    ("14", "2", ["CDC-D3-F-14", "CDC-D3-R-14"]),
    ("16", "2", ["CDC-D3-F-16", "CDC-D3-R-16"]),
    ("18", "2", ["CDC-D3-F-18", "CDC-D3-R-18"]),
    ("20", "2", ["CDC-D3-F-20", "CDC-D3-R-20"]),
    ("22", "2", ["CDC-D3-F-22", "CDC-D3-R-22"]),
    ("24", "2", ["CDC-D3-F-24", "CDC-D3-R-24"]),
    ("26", "2", ["CDC-D3-F-26", "CDC-D3-R-26"]),
    ("28", "2", ["CDC-D3-F-28", "CDC-D3-R-28"]),
    ("30", "2", ["CDC-D3-F-30", "CDC-D3-R-30"]),
    ("32", "2", ["CDC-D3-F-32", "CDC-D3-R-32"]),
    ("34", "2", ["CDC-D3-F-34", "CDC-D3-R-34"]),
    ("6", "3", ["CDC-D3-F-6", "CDC-D3-R-6"]),
    ("11", "3", ["CDC-D3-F-11", "CDC-D3-R-11"]),
    ("13", "3", ["CDC-D3-F-13", "CDC-D3-R-13"]),
    ("19", "3", ["CDC-D3-F-19", "CDC-D3-R-19"]),
    ("25", "3", ["CDC-D3-F-25", "CDC-D3-R-25"]),
    ("27", "3", ["CDC-D3-F-27", "CDC-D3-R-27"]),
    ("29", "3", ["CDC-D3-F-29", "CDC-D3-R-29"]),
    ("31", "3", ["CDC-D3-F-31", "CDC-D3-R-31"]),
    ("35", "3", ["CDC-D3-F-35", "CDC-D3-R-35"]),
]

def to_regex(seq):
    return "".join(IUPAC_REGEX.get(b.upper(), b) for b in seq)

def revcomp(seq):
    return "".join(COMP.get(b.upper(), 'N') for b in reversed(seq))

def read_fasta(path):
    seqs = {}
    name, buf = None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf)
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
        if name:
            seqs[name] = "".join(buf)
    return seqs

def read_primers(path):
    primers = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            name, seq, pool, direction = parts[0], parts[1], parts[2], parts[3]
            primers.append({"name": name, "seq": seq.upper(), "pool": pool, "dir": direction.upper()})
    return primers

def find_exact(genome_seq, primer_seq):
    m = re.search(to_regex(primer_seq), genome_seq, re.IGNORECASE)
    return (m.start(), m.end()) if m else None

def count_mismatches(primer_seq, window):
    return sum(1 for p, g in zip(primer_seq, window) if g.upper() not in IUPAC_SET.get(p.upper(), set()))

def mismatch_positions(primer_seq, window):
    """Retourne la liste des indices (0-based, depuis le 5' de l'amorce TELLE QUE COMPARÉE) qui different."""
    return [i for i, (p, g) in enumerate(zip(primer_seq, window)) if g.upper() not in IUPAC_SET.get(p.upper(), set())]

def best_site_with_mismatches(genome_seq, primer_seq):
    L = len(primer_seq)
    best = (None, L + 1, "")
    for i in range(0, len(genome_seq) - L + 1):
        window = genome_seq[i:i+L]
        mm = count_mismatches(primer_seq, window)
        if mm < best[1]:
            best = (i, mm, window)
            if mm == 0:
                break
    return best

def find_primer(genome_seq, primer):
    hit = find_exact(genome_seq, primer["seq"])
    if hit:
        return {"start": hit[0], "end": hit[1], "strand": "+", "mismatches": 0, "found": True, "mm_positions_5to3": []}
    hit = find_exact(genome_seq, revcomp(primer["seq"]))
    if hit:
        return {"start": hit[0], "end": hit[1], "strand": "-", "mismatches": 0, "found": True, "mm_positions_5to3": []}

    L = len(primer["seq"])
    start_f, mm_f, window_f = best_site_with_mismatches(genome_seq, primer["seq"])
    start_r, mm_r, window_r = best_site_with_mismatches(genome_seq, revcomp(primer["seq"]))

    if mm_f <= mm_r:
        # comparaison faite dans l'orientation F (= 5'->3' d'origine) : positions directement utilisables
        positions = mismatch_positions(primer["seq"], window_f)
        return {"start": start_f, "end": start_f + L, "strand": "+", "mismatches": mm_f,
                "found": False, "mm_positions_5to3": [p + 1 for p in positions]}  # +1 pour affichage 1-based
    else:
        # comparaison faite sur revcomp(amorce) : il faut retourner les indices dans le repere 5'->3' d'origine
        positions_in_revcomp = mismatch_positions(revcomp(primer["seq"]), window_r)
        positions_5to3 = [L - 1 - p for p in positions_in_revcomp]  # miroir
        return {"start": start_r, "end": start_r + L, "strand": "-", "mismatches": mm_r,
                "found": False, "mm_positions_5to3": [p + 1 for p in sorted(positions_5to3)]}

def map_primers_to_genome(primers, genome_seq, genome_name):
    results = []
    for p in primers:
        r = find_primer(genome_seq, p)
        results.append({**p, **r})
    n_exact = sum(1 for r in results if r["found"])
    print(f"\n=== {genome_name} ===")
    print(f"Amorces avec correspondance exacte (IUPAC) : {n_exact} / {len(primers)}")
    approx = [r for r in results if not r["found"]]
    for r in sorted(approx, key=lambda x: x["mismatches"]):
        L = len(next(p["seq"] for p in primers if p["name"] == r["name"]))
        pos_str = ", ".join(f"{p}/{L} depuis le 5'" for p in r["mm_positions_5to3"])
        proche_3prime = any(p >= L - 4 for p in r["mm_positions_5to3"])
        alerte = "  <-- PROCHE DU 3' (plus penalisant)" if proche_3prime else ""
        print(f"  {r['name']:20s} pool {r['pool']} {r['dir']} — {r['mismatches']} mismatch(es), position dans l'amorce : {pos_str}{alerte}")
    return results


def sanitize(name):
    """Nettoie un nom de sequence pour en faire un nom de fichier valide."""
    return re.sub(r'[^A-Za-z0-9_.-]', '_', name)

def draw_scheme(ref_results, ref_name, ref_len, outfile="primer_scheme_diagram.png"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    from matplotlib.lines import Line2D

    by_name = {r["name"]: r for r in ref_results}
    pool_colors = POOL_COLORS
    y_by_pool = {"1": 2, "2": 1, "3": 0}

    plt.rcParams["font.family"] = "DejaVu Sans"
    fig, ax = plt.subplots(figsize=(FIGURE_WIDTH, FIGURE_HEIGHT))
    fig.patch.set_facecolor("white")

    bar_h = 0.46
    any_missing = []
    all_cross_x = []
    all_cross_y = []
    for key, pool, names in AMPLICON_GROUPS:
        hits = [by_name[n] for n in names if n in by_name]
        if len(hits) < 2:
            any_missing.append(key)
            continue
        x0 = min(h["start"] for h in hits)
        x1 = max(h["end"] for h in hits)
        y = y_by_pool.get(pool, 3)
        total_mismatches = sum(h["mismatches"] for h in hits)  # somme (ex: trio F-36/R-36Bis/R-8)
        face = pool_colors.get(pool, "#AAAAAA")

        ax.add_patch(Rectangle((x0, y - bar_h/2), x1 - x0, bar_h,
                                facecolor=face, edgecolor="none", zorder=3))
        if (x1 - x0) > ref_len * 0.02:
            ax.text((x0 + x1) / 2, y, key, ha="center", va="center",
                     fontsize=7.5, color=AMPLICON_TEXT_COLOR, fontweight="normal", zorder=4)
        amplicon_bp = x1 - x0
        ax.text((x0 + x1) / 2, y - bar_h/2 - 0.10, f"{amplicon_bp} pb", ha="center", va="top",
                 fontsize=6, color=LENGTH_LABEL_COLOR, fontstyle="italic", zorder=4)

        if total_mismatches > 0:
            # une croix par mismatch, reparties le long de la largeur de l'amplicon (avec marge)
            margin = (x1 - x0) * 0.12
            xs_start, xs_end = x0 + margin, x1 - margin
            if total_mismatches == 1:
                xs = [(xs_start + xs_end) / 2]
            else:
                xs = [xs_start + i * (xs_end - xs_start) / (total_mismatches - 1) for i in range(total_mismatches)]
            all_cross_x.extend(xs)
            all_cross_y.extend([y + bar_h/2 + 0.11] * len(xs))

    if any_missing:
        print(f"AVERTISSEMENT — amplicons non dessines (amorce manquante dans primers.tsv) : {any_missing}")

    # une petite croix rouge par mismatch (au lieu d'un point unique), pour visualiser la severite
    if all_cross_x:
        ax.scatter(all_cross_x, all_cross_y, color=MISMATCH_COLOR, marker="x", s=MISMATCH_MARKER_SIZE,
                   linewidths=0.8, zorder=5)

    # ligne de fond (backbone) unique, discrete, reliant les 3 rangees visuellement
    for y in [0, 1, 2]:
        ax.plot([0, ref_len], [y, y], color=BACKBONE_COLOR, lw=6, zorder=1, solid_capstyle="round")

    # regle du genome, epuree (graduation tous les 2000 pb seulement)
    ax.plot([0, ref_len], [-0.65, -0.65], color=BACKBONE_COLOR, lw=1, zorder=1)
    step = 2000
    for tick in range(0, ref_len + 1, step):
        ax.text(tick, -0.85, f"{tick:,}", ha="center", va="top", fontsize=7.5, color=AXIS_LABEL_COLOR)

    ax.set_xlim(-50, ref_len + 50)
    ax.set_ylim(-1.05, 2.45)
    ax.set_yticks([0, 1, 2])
    ax.set_yticklabels(["Pool 3", "Pool 2", "Pool 1"], fontsize=9.5, fontweight="normal", color=TITLE_COLOR)
    ax.tick_params(axis="y", length=0)
    ax.set_xticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

    ax.set_title(f"{ref_name} — schéma d'amplicons tuilés ({ref_len:,} pb)",
                 fontsize=11.5, fontweight="normal", color=TITLE_COLOR, pad=10, loc="left")

    legend_elems = [
        Line2D([0], [0], marker="s", color="none", markerfacecolor=POOL_COLORS["1"], markersize=10, label="Pool 1"),
        Line2D([0], [0], marker="s", color="none", markerfacecolor=POOL_COLORS["2"], markersize=10, label="Pool 2"),
        Line2D([0], [0], marker="s", color="none", markerfacecolor=POOL_COLORS["3"], markersize=10, label="Pool 3"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=MISMATCH_COLOR, markersize=6, label="≥1 mismatch"),
    ]
    ax.legend(handles=legend_elems, loc="upper center", bbox_to_anchor=(0.5, -0.14),
              ncol=4, frameon=False, fontsize=8.5, handletextpad=0.4, columnspacing=1.5)

    plt.tight_layout()
    plt.savefig(outfile, dpi=200, facecolor="white", bbox_inches="tight")
    print(f"Diagramme sauvegarde : {outfile}")


def main():
    if len(sys.argv) < 3:
        print(__doc__); sys.exit(1)
    primers = read_primers(sys.argv[1])
    print(f"{len(primers)} amorces chargees depuis {sys.argv[1]}")
    ref_seqs = read_fasta(sys.argv[2])
    ref_name, ref_seq = next(iter(ref_seqs.items()))
    ref_len = len(ref_seq)
    print(f"Reference : {ref_name} ({ref_len} pb)")
    ref_results = map_primers_to_genome(primers, ref_seq, ref_name)

    with open("primer_positions_reference_v4.tsv", "w") as f:
        f.write("name\tpool\tdirection\tstart\tend\tstrand\tmismatches\texact_match\n")
        for r in ref_results:
            f.write(f"{r['name']}\t{r['pool']}\t{r['dir']}\t{r['start']}\t{r['end']}\t{r['strand']}\t{r['mismatches']}\t{r['found']}\n")
    print("Coordonnees sauvegardees dans primer_positions_reference_v4.tsv")

    try:
        draw_scheme(ref_results, ref_name, ref_len, outfile=f"primer_scheme_diagram_{sanitize(ref_name)}.png")
    except ImportError as e:
        print(f"matplotlib non disponible ({e}) — diagramme non genere")

    for path in sys.argv[3:]:
        seqs = read_fasta(path)
        for name, seq in seqs.items():
            results = map_primers_to_genome(primers, seq, name)
            genome_len = len(seq)
            try:
                draw_scheme(results, name, genome_len, outfile=f"primer_scheme_diagram_{sanitize(name)}.png")
            except ImportError as e:
                print(f"matplotlib non disponible ({e}) — diagramme non genere pour {name}")

if __name__ == "__main__":
    main()
