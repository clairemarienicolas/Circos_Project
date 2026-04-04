"""
Script : script3_genedensity.py

Description :
Ce script génère un fichier de densité des gènes à partir d’un fichier
d'annotation GFF3. Il compte le nombre de gènes présents dans chaque fenêtre 
de taille fixe le long du génome.

Le script :
- lit le fichier GFF3 (via md.parse_gff3_genes)
- filtre uniquement les entités de type 'gene'
- ignore les contigs et scaffolds pour ne garder que les chromosomes (ex: Chr01)
- découpe chaque chromosome en fenêtres de 100 kb
- compte le nombre de gènes dont la position est incluse dans chaque fenêtre
- écrit un fichier texte tabulé compatible avec Circos

Format de sortie :
Chr01    1       100000    12
Chr01    100001  200000    5

Utilisation :
python scripts/script3_genedensity.py input.gff output.txt

Exemple :
python scripts/script3_genedensity.py data/Sbicolor_313_v3.1.gene_exons.gff3 results/gene_density.txt
"""

import sys
import module as md

def calculate_gene_density(input_file, output_file):
    """
    Lit un fichier GFF3 et calcule, pour chaque fenêtre de 100 kb,
    le nombre de gènes présents.

    Paramètres :
    - input_file  : fichier GFF3 en entrée
    - output_file : fichier texte tabulé en sortie
    """
    window_size = 100000

    # Lecture du GFF3 : md.parse_gff3_genes retourne :
    #   gene_positions : {gene_id: (chrom, start)}
    #   chr_lengths : {chrom: longueur_max}
    gene_positions, chr_lengths = md.parse_gff3_genes(input_file)

    with open(output_file, 'w') as out:
        for chrom in sorted(chr_lengths.keys()):

            # Génération des fenêtres : md.get_windows retourne [(start, end), ...]
            for start, end in md.get_windows(chr_lengths[chrom], window_size):

                # Gènes dans la fenêtre : md.get_genes_in_window retourne une liste de gene_id
                genes_in_window = md.get_genes_in_window(gene_positions, chrom, start, end)
                count = len(genes_in_window)

                out.write(f"{chrom}\t{start + 1}\t{end}\t{count}\n")

    print(f"Fichier {output_file} créé avec succès.")


def main():
    if len(sys.argv) != 3:
        print("Usage : script3_genedensity.py input_filename output_filename")
        sys.exit(1)

    input_file  = md.get_input_filename()
    output_file = md.get_output_filename()

    calculate_gene_density(input_file, output_file)


if __name__ == "__main__":
    main()
