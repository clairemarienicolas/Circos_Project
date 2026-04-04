#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -M claire-marie.nicolas@etudiant.univ-rennes.fr
#$ -m bea

# 1. Charger Python
module load python/3.9.5

# 2. Lancer les scripts
# L'option "set -e" fait stopper le script dès qu'une commande échoue,
# au lieu de continuer silencieusement avec des fichiers manquants ou corrompus.
set -e

echo "Génération du karyotype..."
python3 scripts/script1_chr.py Sb_data/Sbicolor_313_v3.0.fa results/karyotype.txt

echo "Calcul de la densité des TE..."
python3 scripts/script2_TEdensity.py Sb_data/Sbicolor_313_v3.0.hardmasked.fa results/TE_density.txt

echo "Calcul de la densité des gènes..."
python3 scripts/script3_genedensity.py Sb_data/Sbicolor_313_v3.1.gene_exons.gff3 results/gene_density.txt

echo "Calcul du nombre moyen d'exons par gène et par fenêtre..."
python3 scripts/script4_nbexons.py Sb_data/Sbicolor_313_v3.1.gene_exons.gff3 results/nb_exons.txt

echo "Calcul de la longueur moyenne des exons par fenêtre..."
python3 scripts/script5_lenexons.py Sb_data/Sbicolor_313_v3.1.gene_exons.gff3 results/len_exons.txt

echo "Traitement des blocs de synténie..."
python3 scripts/script6_syntenie.py Sb_data/Sb_Sb.aligncoords.gcoords.txt results/links.txt

echo "Terminé !"
