#!/usr/bin/env python3
"""
Trace la profondeur de couverture le long du genome, avec
surlignage des zones fragiles identifiees.

Usage dans le tableau de commande : python3 plot_depth_igv_style_v2.py barcode01.depth.txt [barcode02.depth.txt ...]
Genere un PNG par fichier fourni.
"""
import sys

# =========================================================================
# CONFIGURATION — modifie librement les valeurs ci-dessous
# =========================================================================
FIGURE_WIDTH = 14
FIGURE_HEIGHT = 4.2
DPI = 200

TITLE_PREFIX = "Profondeur de couverture"   # suivi automatiquement de " — <nom echantillon>"
XLABEL = "Position sur MH544651.1 (pb)"
YLABEL = "Profondeur (X)"

# Couleurs
DEPTH_FILL_COLOR = "#B0B0AC"      # remplissage de l'aire de couverture
DEPTH_LINE_COLOR = "#8A8A85"      # contour fin de la courbe
THRESHOLD_COLOR  = "#5A8FB0"      # ligne du seuil min-depth
TITLE_COLOR      = "#3A3A38"
AXIS_LABEL_COLOR = "#3A3A38"
TICK_COLOR       = "#3A3A38"
GRID_COLOR       = "#EDEDEA"

# Seuil de profondeur minimale du pipeline (ligne pointillee horizontale)
MIN_DEPTH_THRESHOLD = 20
MIN_DEPTH_LABEL = "seuil min-depth=20"

# Zones fragiles a surligner : name, start, end, color, note (note = legende courte)
FRAGILE_ZONES = [
    {"name": "Amplicon 25", "start": 6770, "end": 7386, "color": "#E8A33D",
     "note": "0 mismatch sur  amorces ref. et Guyanaise"},
    {"name": "Amplicon 36", "start": 10124, "end": 10706, "color": "#C0392B",
     "note": "mismatch des 2 amorces reverse du trio sur la souche Guyanaise et ref. "},
]
# =========================================================================

def read_depth(path):
    positions, depths = [], []
    with open(path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            try:
                pos = int(parts[1])
                d = int(parts[2])
            except ValueError:
                continue
            positions.append(pos)
            depths.append(d)
    return positions, depths

def plot_one(depth_path, outfile):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    positions, depths = read_depth(depth_path)
    if not positions:
        print(f"AVERTISSEMENT : aucune donnee lue dans {depth_path}")
        return
    genome_len = max(positions)
    max_depth = max(depths) if depths else 1

    fig, ax = plt.subplots(figsize=(FIGURE_WIDTH, FIGURE_HEIGHT))
    fig.patch.set_facecolor("white")

    ax.fill_between(positions, depths, color=DEPTH_FILL_COLOR, linewidth=0, zorder=2)
    ax.plot(positions, depths, color=DEPTH_LINE_COLOR, linewidth=0.4, zorder=3)

    ax.axhline(MIN_DEPTH_THRESHOLD, color=THRESHOLD_COLOR, linewidth=1, linestyle="--", zorder=4, alpha=0.8)
    ax.text(genome_len * 0.995, MIN_DEPTH_THRESHOLD, f" {MIN_DEPTH_LABEL}", va="bottom", ha="right",
            fontsize=7.5, color=THRESHOLD_COLOR)

    # Surlignage des zones fragiles : bandeau + nom + note, TOUT regroupe en haut
    # (evite toute collision avec le label de l'axe des abscisses, place lui tout en bas)
    for zone in FRAGILE_ZONES:
        ax.axvspan(zone["start"], zone["end"], color=zone["color"], alpha=0.15, zorder=1)
        ax.add_patch(Rectangle((zone["start"], max_depth * 1.03), zone["end"] - zone["start"],
                                max_depth * 0.05, facecolor=zone["color"], edgecolor="none", zorder=5))
        mid = (zone["start"] + zone["end"]) / 2
        ax.text(mid, max_depth * 1.11, zone["name"], ha="center", va="bottom",
                fontsize=9, fontweight="bold", color=zone["color"], zorder=6)
        ax.text(mid, max_depth * 1.20, zone["note"], ha="center", va="bottom",
                fontsize=6.8, color=zone["color"], zorder=6)

    sample_name = depth_path.split("/")[-1].replace(".depth.txt", "")
    ax.set_title(f"{TITLE_PREFIX} — {sample_name}", fontsize=12.5,
                 fontweight="normal", color=TITLE_COLOR, loc="left", pad=14)
    ax.set_xlabel(XLABEL, fontsize=9.5, color=AXIS_LABEL_COLOR, labelpad=8)
    ax.set_ylabel(YLABEL, fontsize=9.5, color=AXIS_LABEL_COLOR)
    ax.set_xlim(0, genome_len)
    ax.set_ylim(-max_depth * 0.06, max_depth * 1.32)   # marge en haut pour les annotations, pas en bas
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    ax.tick_params(labelsize=8, colors=TICK_COLOR)
    ax.grid(axis="y", color=GRID_COLOR, linewidth=0.6, zorder=0)

    plt.tight_layout()
    plt.savefig(outfile, dpi=DPI, facecolor="white", bbox_inches="tight")
    print(f"Figure sauvegardee : {outfile}")

def main():
    if len(sys.argv) < 2:
        print(__doc__); sys.exit(1)
    for path in sys.argv[1:]:
        outfile = path.replace(".depth.txt", "") + "_depth_plot.png"
        plot_one(path, outfile)

if __name__ == "__main__":
    main()
