import sys
from collections import defaultdict
""" Densité des éléments transposables

NOUS ON TRAVAIL AVEC LE FICHIER FASTA CAR ON A SORGHUM BICOLOR
> un tableau avec (cf diapo):
    - une colonne avec les id des chromosomes
    - une colonne avec la position du début de la fenêtre
    - une colonne avec la position de fin de la fenêtre 
    - une colonne avec la valeur de densité -> nombre de N dans l'intervalle de la fenêtre """


# Paramètres
WINDOW_SIZE = 100000  # 100 kb

# Entrée
gff_file = sys.argv[1]
output_file = sys.argv[2]

# Dictionnaire : chr -> liste des positions TE
te_positions = defaultdict(list)

# Lecture du GFF
with open(gff_file) as f:
    for line in f:
        if line.startswith(">"): #check si le symbol est bon !!
            continue
        
        cols = line.strip().split("\t")
        
        chrom = cols[0]
        feature_type = cols[2]
        start = int(cols[3])
        end = int(cols[4])

        # !!!! Adapter si besoin selon le fichier
        if "transpos" in feature_type.lower() or "TE" in feature_type:
            te_positions[chrom].append((start, end))

# Calcul densité
with open(output_file, "w") as out:
    
    for chrom in te_positions:
        # Taille du chromosome = max position observée
        max_pos = max(end for start, end in te_positions[chrom])
        
        # Création des fenêtres
        for window_start in range(0, max_pos, WINDOW_SIZE):
            window_end = window_start + WINDOW_SIZE
            count = 0
            
            for start, end in te_positions[chrom]:
                # chevauchement avec la fenêtre
                if not (end < window_start or start > window_end):
                    count += 1
            
            out.write(f"{chrom}\t{window_start}\t{window_end}\t{count}\n")


