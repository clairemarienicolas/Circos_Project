"""
Script : script5_syntenie.py

Description :
Ce script convertit un fichier de synténie SynMap (.gcoords.txt) en un fichier 
de liens (links) compatible avec Circos. Au lieu de tracer chaque gène 
individuellement, il regroupe les gènes par blocs de synténie (délimités 
par des '#' dans le fichier source) et calcule les coordonnées globales (min/max).

Le script :
- Lit le fichier de coordonnées d'alignement par blocs.
- Extrait les positions de début et de fin de chaque gène via le séparateur "||".
- Calcule l'enveloppe globale (start/end) de chaque bloc de synténie.
- Filtre pour ne garder que les chromosomes principaux (Chr01, etc.).
- Filtre par taille de bloc (distance > 50kb) pour assurer la lisibilité.
- Écrit les 6 colonnes tabulées requises par Circos.

Format de sortie :
Chr01    48286845    52228681    Chr01    71614120    71801583

Utilisation :
python scripts/script5_syntenie.py input.gcoords.txt results/output.txt

Exemple :
python scripts/script5_syntenie.py data/Sb_Sb.aligncoords.gcoords.txt results/links.txt
"""

import sys
import module as md

def format_chr_name(raw_name):
    """
    Transforme un identifiant (ex: '1', 'a31607_1' ou 'super_12') en 'ChrXX'.
    """
    # On cherche le dernier chiffre après un underscore ou on prend le tout
    parts = raw_name.replace('||', '_').split('_')
    last_part = parts[-1]
    try:
        return f"Chr{int(last_part):02d}"
    except ValueError:
        return None

def write_link(out_file, chr_a, block_a, chr_b, block_b, min_size):
    """
    Calcule les bornes d'un bloc et l'écrit s'il dépasse la taille minimale.
    Le filtre est appliqué sur les DEUX côtés (A et B) pour éviter les blocs
    asymétriques aberrants : un bloc doit être suffisamment grand sur les deux
    génomes pour être biologiquement significatif.
    """
    if not block_a or not block_b:
        return 0

    name_a = format_chr_name(chr_a)
    name_b = format_chr_name(chr_b)
    
    start_a, end_a = min(block_a), max(block_a)
    start_b, end_b = min(block_b), max(block_b)

    size_a = end_a - start_a
    size_b = end_b - start_b

    # Filtre de taille sur A ET B : les deux régions doivent dépasser le seuil
    if size_a >= min_size and size_b >= min_size and name_a and name_b:
        out_file.write(f"{name_a}\t{start_a}\t{end_a}\t{name_b}\t{start_b}\t{end_b}\n")
        return 1
    return 0

def transform(input_file, output_file):
    min_dist_size = 200_000  # Filtre à 200kb pour meilleure lisibilité Circos
    links_count = 0
    
    current_block_a = []
    current_block_b = []
    current_chr_a = None
    current_chr_b = None

    with open(input_file, 'r') as f, open(output_file, 'w') as out:
        for line in f:
            line = line.strip()
            if not line or line.startswith("# "): 
                continue

            if line.startswith("#"):
                # On traite le bloc qui vient de se terminer
                links_count += write_link(out, current_chr_a, current_block_a, 
                                          current_chr_b, current_block_b, min_dist_size)
                # Reset
                current_block_a, current_block_b = [], []
                current_chr_a, current_chr_b = None, None
            else:
                parts = line.split('\t')
                if len(parts) < 5: continue

                # Extraction selon le format
                data_a = parts[1].split('||')
                data_b = parts[5].split('||')

                if len(data_a) >= 3 and len(data_b) >= 3:
                    current_chr_a = data_a[0]
                    current_chr_b = data_b[0]
                    current_block_a.extend([int(data_a[1]), int(data_a[2])])
                    current_block_b.extend([int(data_b[1]), int(data_b[2])])

        # Ne pas oublier le dernier bloc après la boucle
        links_count += write_link(out, current_chr_a, current_block_a, 
                                  current_chr_b, current_block_b, min_dist_size)

    print(f"Succès : {links_count} blocs de synténie écrits dans {output_file}")

def main():
    if len(sys.argv) != 3:
        print("Usage : python script5_syntenie.py <input> <output>")
        sys.exit(1)
    
    input_file = md.get_input_filename()
    output_file = md.get_output_filename()
    transform(input_file, output_file)

if __name__ == "__main__":
    main()