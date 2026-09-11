#!/usr/bin/env python3
# =============================================================================
# figure_couverture.py —  DENV-3
#
# Génère deux panneaux :
#   A. carte de chaleur  profondeur médiane par amplicon x barcode
#   B. profil de profondeur position par position
#
# PRÉREQUIS
#   pip install matplotlib numpy
#
# ENTRÉES (2 fichiers seulement)
#   1) le BED du schéma d'amorces
#   2) un TSV à 3 colonnes : barcode <TAB> position <TAB> profondeur
#      produit par :
#        for i in $(seq -w 1 23); do
#          samtools depth -a -d 0 barcode${i}.primertrimmed.rg.sorted.bam \
#            | awk -v s="barcode${i}" '{print s"\t"$2"\t"$3}'
#        done > depth_all_v5.tsv
#
# USAGE
#   python3 figure_couverture.py
# =============================================================================

import re
import statistics
import numpy as np
import matplotlib
matplotlib.use('Agg')          # pas d'affichage interactif : on écrit un fichier
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Rectangle

# ---- À ADAPTER --------------------------------------------------------------
BED    = 'DENV3_scheme.bed'
DEPTH  = 'depth_all_v5.tsv'
SORTIE = 'figure_couverture_v5'     # sans extension : .png et .pdf seront écrits
SEUIL  = 20                         # seuil de masquage artic (--min-depth)
REF_LABEL = 'MH544651.1'
# -----------------------------------------------------------------------------


# =============================================================================
# ÉTAPE 1 — lire le BED et en déduire les bornes d'insert de chaque amplicon
# =============================================================================
# Un BED ARTIC contient une ligne par AMORCE, pas par amplicon. Pour chaque
# amplicon n il y a au moins DENV3_n_LEFT et DENV3_n_RIGHT (parfois des
# variantes _alt). L'INSERT — la partie réellement séquencée après primer
# trimming — va de la FIN de l'amorce LEFT au DÉBUT de l'amorce RIGHT.
# Quand plusieurs variantes existent on prend la borne la plus interne :
#   max() des fins LEFT, min() des débuts RIGHT.

amp = {}
for ligne in open(BED):
    champs = ligne.rstrip('\n').split('\t')
    if len(champs) < 6:
        continue
    debut, fin, nom = int(champs[1]), int(champs[2]), champs[3]

    # adapter cette regex si tes amorces sont nommées autrement
    # (ex. CNR-D3-F-25 / CNR-D3-R-25  ->  r'CNR-D3-([FR])-(\d+)')
    m = re.match(r'DENV3_(\d+)_(LEFT|RIGHT)', nom)
    if not m:
        print(f"  [!] nom d'amorce non reconnu, ignoré : {nom}")
        continue

    n, cote = int(m.group(1)), m.group(2)
    a = amp.setdefault(n, {'L': [], 'R': []})
    a['L' if cote == 'LEFT' else 'R'].append((debut, fin))

ins = {n: (max(f for d, f in a['L']), min(d for d, f in a['R']))
       for n, a in amp.items() if a['L'] and a['R']}
amps = sorted(ins)
print(f"Schéma : {len(amps)} amplicons, de {min(amps)} à {max(amps)}")


# =============================================================================
# ÉTAPE 2 — charger les profondeurs
# =============================================================================
# Structure : depth['barcode01'][position] = profondeur
# ~250 000 lignes, quelques secondes.

depth = {}
for ligne in open(DEPTH):
    bc, pos, d = ligne.split('\t')
    depth.setdefault(bc, {})[int(pos)] = int(d)

bcs = sorted(depth)
GLEN = max(max(v) for v in depth.values())
print(f"Profondeurs : {len(bcs)} échantillons, génome de {GLEN} nt")


# =============================================================================
# ÉTAPE 3 — matrice barcode x amplicon
# =============================================================================
# M[i,j]   = profondeur MÉDIANE de l'amplicon j chez le barcode i.
#            La médiane est préférable à la moyenne : elle n'est pas tirée
#            vers le bas par les quelques bases de bord d'amplicon.
# LOW[i,j] = True si plus de 5 % des positions passent sous le seuil, donc
#            si des N apparaîtront dans le consensus.
# Le max(..., 0.5) évite un zéro, impossible à représenter en échelle log.

M   = np.zeros((len(bcs), len(amps)))
LOW = np.zeros_like(M, dtype=bool)

for i, bc in enumerate(bcs):
    D = depth[bc]
    for j, n in enumerate(amps):
        s, e = ins[n]
        v = [D.get(p, 0) for p in range(s + 1, e + 1)]   # BED est 0-based
        M[i, j]   = max(statistics.median(v), 0.5)
        LOW[i, j] = sum(1 for x in v if x < SEUIL) / len(v) > 0.05


# =============================================================================
# ÉTAPE 4 — tracer
# =============================================================================
fig = plt.figure(figsize=(15, 11))
gs  = fig.add_gridspec(2, 1, height_ratios=[1.35, 1], hspace=0.30)

# ---- Panneau A : carte de chaleur -------------------------------------------
# LogNorm est indispensable : la profondeur va de 0 à ~800X. En échelle
# linéaire, tout ce qui est sous 100X apparaît uniformément sombre et les
# chutes deviennent invisibles.

ax = fig.add_subplot(gs[0])
im = ax.imshow(M, aspect='auto', cmap='YlGnBu',
               norm=LogNorm(vmin=1, vmax=max(1000, M.max())))

ax.set_xticks(range(len(amps)))
ax.set_xticklabels(amps, fontsize=8)
ax.set_yticks(range(len(bcs)))
ax.set_yticklabels([b.replace('barcode', 'bc') for b in bcs], fontsize=8)
ax.set_xlabel("Amplicon", fontsize=11)
ax.set_title(f"A.  Profondeur médiane par amplicon  —  {len(bcs)} isolats DENV-3, "
             f"{len(amps)} amplicons",
             fontsize=12, loc='left', pad=12)

# encadrer en rouge les couples (échantillon, amplicon) qui produiront des N
for i in range(len(bcs)):
    for j in range(len(amps)):
        if LOW[i, j]:
            ax.add_patch(Rectangle((j - .5, i - .5), 1, 1, fill=False,
                                   edgecolor='#c1272d', lw=1.8))

cb = fig.colorbar(im, ax=ax, pad=0.012, fraction=0.025)
cb.set_label("Profondeur (X, échelle log)", fontsize=9)

# entrée de légende factice : imshow ne génère pas de légende tout seul
ax.plot([], [], 's', mfc='none', mec='#c1272d', mew=1.8, ms=9,
        label=f"> 5 % des positions sous {SEUIL}X")
ax.legend(loc='upper left', bbox_to_anchor=(0, -0.09), frameon=False, fontsize=9)

# ---- Panneau B : profil positionnel -----------------------------------------
# 23 courbes fines et claires en arrière-plan (la dispersion entre échantillons),
# puis la médiane en trait épais foncé par-dessus.

ax2 = fig.add_subplot(gs[1])

pos = np.arange(1, GLEN + 1)
mat = np.zeros((len(bcs), GLEN))
for i, bc in enumerate(bcs):
    D = depth[bc]
    mat[i] = [D.get(p, 0) for p in pos]

for i in range(len(bcs)):
    ax2.plot(pos, np.maximum(mat[i], .3), color='#9fb8c8', lw=.4, alpha=.55)

med = np.median(mat, axis=0)
ax2.plot(pos, np.maximum(med, .3), color='#14415c', lw=1.3,
         label=f'Médiane des {len(bcs)} isolats')
ax2.axhline(SEUIL, color='#c1272d', ls='--', lw=1.2,
            label=f'Seuil de masquage ({SEUIL}X)')

ax2.set_yscale('log')
ax2.set_ylim(.3, max(2000, mat.max() * 1.5))
ax2.set_xlim(0, GLEN)
ax2.set_xlabel(f"Position sur le génome ({REF_LABEL}, nt)", fontsize=11)
ax2.set_ylabel("Profondeur (X)", fontsize=11)
ax2.set_title("B.  Profil de profondeur par position", fontsize=12, loc='left', pad=10)

# zones grises : hors couverture du schéma, avant le 1er amplicon et après le dernier.
# Ce sont des N structurels, pas des échecs de séquençage — à distinguer clairement.
a1s  = ins[amps[0]][0]
a36e = ins[amps[-1]][1]
for (x0, x1) in [(0, a1s), (a36e, GLEN)]:
    ax2.axvspan(x0, x1, color='#bbbbbb', alpha=.45, zorder=0)
ax2.text(a36e + (GLEN - a36e) / 2, 3.5, "3′ hors schéma", ha='center',
         fontsize=7.5, color='#555555', rotation=90, va='bottom')
ax2.text(a1s / 2, 3.5, "5′", ha='center', fontsize=7.5, color='#555555',
         rotation=90, va='bottom')

# amplicons à commenter — À ADAPTER selon ce que tu veux mettre en avant
A_SIGNALER = [(36, '#c1272d', "ampl. 36\nchute résiduelle")]
for n, couleur, libelle in A_SIGNALER:
    if n not in ins:
        continue
    s, e = ins[n]
    ax2.axvspan(s, e, color=couleur, alpha=.13, zorder=0)
    ax2.text((s + e) / 2, 1100, libelle, ha='center', fontsize=8,
             color=couleur, weight='bold')

ax2.legend(loc='lower left', fontsize=9, framealpha=.92)
ax2.grid(axis='y', ls=':', lw=.5, alpha=.55)

fig.suptitle("Couverture de séquençage par amplicon",
             fontsize=13.5, y=0.965, weight='bold')


# =============================================================================
# ÉTAPE 5 — écrire figure et tableau
# =============================================================================
# PNG pour Word, PDF vectoriel pour l'impression et la soumission.
for ext in ('png', 'pdf'):
    fig.savefig(f'{SORTIE}.{ext}', dpi=300, bbox_inches='tight', facecolor='white')
    print(f"écrit : {SORTIE}.{ext}")

# le tableau brut, à mettre en annexe
with open(f'{SORTIE}.tsv', 'w') as fh:
    fh.write("amplicon\tdebut\tfin\tlongueur\t" +
             "\t".join(b.replace('barcode', 'bc') for b in bcs) + "\n")
    for j, n in enumerate(amps):
        s, e = ins[n]
        fh.write(f"{n}\t{s}\t{e}\t{e-s}\t" +
                 "\t".join(f"{M[i,j]:.0f}" for i in range(len(bcs))) + "\n")
print(f"écrit : {SORTIE}.tsv")


# =============================================================================
# ÉTAPE 6 — résumé texte, à recopier dans le mémoire
# =============================================================================
print("\n" + "=" * 64)
print("AMPLICONS À SURVEILLER")
print("=" * 64)
for j, n in enumerate(amps):
    s, e = ins[n]
    nlow = int(LOW[:, j].sum())
    if nlow == 0:
        continue
    concernes = [bcs[i].replace('barcode', 'bc') for i in range(len(bcs)) if LOW[i, j]]
    print(f"amplicon {n:>2} ({e-s} pb, {s}-{e}) : médiane {np.median(M[:, j]):.0f}X"
          f" — {nlow} échantillon(s) sous seuil : {', '.join(concernes)}")
print("=" * 64)
