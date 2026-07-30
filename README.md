# Dengue3S-quen-age

Vérification d'appariement d'amorces — Schéma DENV-3 (CDC / CNR Guyane)

Script Python de vérification informatique de l'appariement d'un schéma d'amorces PCR multiplexées (amplicons tuilés) sur un génome de référence

Objectif

Le protocole ARTIC de séquençage MinION en amplicons tuilés repose sur un jeu d'amorces PCR réparties en plusieurs pools (ici : 3 pools, 73 amorces, 36 amplicons). Certaines de ces amorces contiennent des positions dégénérées (codes IUPAC : R, Y, S, W, K, M, B, D, H, V, N) pour couvrir la diversité génétique du virus. Ce script vérifie, pour chaque amorce :

si elle s'apparie exactement (compte tenu des dégénérescences) sur la séquence de référence donnée ;
si non, la position du meilleur site possible et le nombre de mésappariements (mismatches) à cet endroit ;
la position exacte du ou des mismatch(es) au sein de l'amorce, avec repérage de leur proximité à l'extrémité 3' (la plus déterminante pour l'efficacité d'extension par la polymérase) ;
génère un diagramme visuel du schéma d'amplicons tuilés, coloré par pool, avec signalement des amplicons présentant un mismatch.

Le script peut aussi être utilisé pour tester l'appariement du même jeu d'amorces sur d'autres génomes (autres souches, autres génotypes), afin d'évaluer le risque de "dropout" d'amplicon sur des séquences divergentes.

Origine du schéma d'amorces

Les séquences d'amorces vérifiées ici sont celles décrites dans :

Lagrave A, Enfissi A, Tirera S, et al. Re-Emergence of DENV-3 in French Guiana: Retrospective Analysis of Cases That Circulated in the French Territories of the Americas from the 2000s to the 2023–2024 Outbreak. Viruses. 2024;16(8):1298. https://doi.org/10.3390/v16081298

Schéma dérivé du protocole CDC (Santiago et al., PLoS Negl Trop Dis 2013) et adapté par le CNR Guyane pour le DENV-3.

Prérequis
Python ≥ 3.9
matplotlib (pour la génération du diagramme)
bash
pip install matplotlib --break-system-packages
Fichiers d'entrée
1. Fichier des amorces (primers.tsv)

Fichier texte à séparateur tabulation, 4 colonnes, sans en-tête obligatoire (les lignes commençant par # sont ignorées) :

nom_amorce	sequence	pool	sens
DENV3-F-1	AGTTGTTAGTCTACGTGGACCGAC	1	F
CDC-D3-R-1	ATCATSCGYGGCTCTCCAT	1	R
sens : F (forward) ou R (reverse)
Les séquences peuvent contenir des codes IUPAC dégénérés (R, Y, S, W, K, M, B, D, H, V, N)
2. Génome(s) de référence (FASTA)

Séquence(s) génomique(s) contre laquelle/lesquelles tester l'appariement, au format FASTA standard.

Utilisation
bash
python3 primer_mapping.py primers.tsv reference.fasta [autre_genome1.fasta autre_genome2.fasta ...]
Le premier génome fourni est traité comme la référence principale (utilisée pour générer le diagramme).
Tout génome supplémentaire est testé avec le même jeu d'amorces.
Exemple
bash
efetch -db nucleotide -id "MH544651.1" -format fasta > MH544651_ref.fasta
python3 primer_mapping.py primers.tsv MH544651_ref.fasta
Sorties générées
Console : pour chaque génome testé, nombre d'amorces avec correspondance exacte, et détail des amorces en mismatch (position dans l'amorce, distance à l'extrémité 3').
primer_positions_reference_vX.tsv : coordonnées, brin, nombre de mismatches et statut (correspondance exacte ou non) pour chaque amorce sur le génome de référence.
primer_scheme_diagram_vX.png : représentation visuelle du schéma d'amplicons tuilés, un rectangle par amplicon, coloré par pool, avec un marqueur pour les amplicons comportant un mismatch.

Les couleurs et dimensions du diagramme sont modifiables dans le bloc CONFIGURATION en tête de script.

Principe de l'algorithme
Traduction IUPAC → expression régulière : chaque code dégénéré de l'amorce est converti en classe de caractères équivalente (ex. S → [GC]).
Recherche exacte : recherche par expression régulière sur le brin direct, puis sur le brin complémentaire.
Recherche approchée (si l'étape 2 échoue) : comparaison glissante de l'amorce le long de l'intégralité du génome, position par position, avec décompte des mésappariements à chaque position testée. Le site présentant le nombre minimal de mésappariements est retenu comme "meilleur site".
Localisation du mismatch dans l'amorce : les positions divergentes sont recalculées dans le référentiel 5'→3' d'origine de l'amorce (indépendamment du brin sur lequel la correspondance a été trouvée), pour évaluer leur proximité à l'extrémité 3'.
Limites connues
La recherche approchée (étape 3) a une complexité en O(longueur du génome × longueur de l'amorce) par amorce : adaptée à des génomes de quelques dizaines de kb, non optimisée pour des génomes de grande taille.
Le regroupement des amorces en amplicons repose sur une liste explicite (AMPLICON_GROUPS) propre à ce schéma précis (73 amorces, 36 amplicons, 3 pools) ; toute utilisation avec un autre schéma d'amorces nécessite d'adapter cette liste.
Un mismatch proche de l'extrémité 3' est signalé comme potentiellement plus pénalisant pour l'efficacité PCR, mais le script ne modélise pas quantitativement l'impact réel sur l'efficacité d'amplification (pas de calcul de Tm ni d'énergie libre d'hybridation).
Contexte du projet

Ce script a été développé dans le cadre du travail de mémoire de DU (les agents infectieux à l'ère de la génomique) portant sur la surveillance génomique du DENV-3 en Guadeloupe


Marwan
