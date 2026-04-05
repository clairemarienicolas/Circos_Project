"""
Script : script4_nbexons.py

Description :
Ce script génère un fichier de données à partir d'un fichier d'annotation GFF3.
Il calcule, pour chaque fenêtre de 500 kb le long du génome,
le nombre moyen d'exons par gène.

Le script :
- lit le fichier GFF3 (via md.parse_gff3_genes et md.parse_gff3_exons)
- filtre uniquement les entités de type 'gene', 'mRNA' et 'exon'
- ignore les contigs et scaffolds pour ne garder que les chromosomes (ex: Chr01)
- associe chaque exon à son gène parent via les lignes mRNA
- découpe chaque chromosome en fenêtres de 500 kb
- calcule le nombre moyen d'exons par gène dans chaque fenêtre
- écrit un fichier texte tabulé compatible avec Circos

Format de sortie :
Chr01    1       500000    5.30
Chr01    500001  1000000   4.80

Utilisation :
python scripts/script4_nbexons.py input.gff3 output.txt

Exemple :
python scripts/script4_nbexons.py data/Sbicolor_313_v3.1.gene_exons.gff3 results/nb_exons.txt
"""

import sys
import module as md

def calculate_mean_nb_exons_per_gene(input_file, output_file):
    """
    Lit un fichier GFF3 et calcule, pour chaque fenêtre de 500 kb,
    le nombre moyen d'exons par gène.

    Paramètres :
    - input_file  : fichier GFF3 en entrée
    - output_file : fichier texte tabulé en sortie
    """
    window_size = 500000 #Fenêtre de 500000 pour avoir des moyennes plus "représentatives" 
    # Car peu de gènes par fenêtre de 100000

    # Lecture des gènes : md.parse_gff3_genes retourne :
    #   gene_positions : {gene_id: (chrom, start)}
    #   chr_lengths : {chrom: longueur_max}
    gene_positions, chr_lengths = md.parse_gff3_genes(input_file)

    # Lecture des exons : md.parse_gff3_exons retourne :
    #   gene_exon_count : {gene_id: nombre d'exons}
    #   gene_exon_lengths : {gene_id: [longueurs des exons]} (non utilisé ici)
    gene_exon_count, _ = md.parse_gff3_exons(input_file)

    with open(output_file, 'w') as out:
        for chrom in sorted(chr_lengths.keys()):

            for start, end in md.get_windows(chr_lengths[chrom], window_size):
                genes_in_window = md.get_genes_in_window(gene_positions, chrom, start, end)

                # On ignore les fenêtres sans gènes (bords de chromosome)
                if not genes_in_window:
                    continue

                nb_exons_list = [gene_exon_count[g] for g in genes_in_window
                                 if g in gene_exon_count and gene_exon_count[g] > 0]
                if nb_exons_list:
                    mean_nb_exons = sum(nb_exons_list) / len(nb_exons_list)
                    out.write(f"{chrom}\t{start + 1}\t{end}\t{mean_nb_exons:.2f}\n")

    print(f"Fichier {output_file} créé avec succès.")


def main():
    if len(sys.argv) != 3:
        print("Usage : script4_nbexons.py input.gff3 output.txt")
        sys.exit(1)

    input_file = md.get_input_filename()
    output_file = md.get_output_filename()

    calculate_mean_nb_exons_per_gene(input_file, output_file)


if __name__ == "__main__":
    main()
