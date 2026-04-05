"""
Script : script5_lenexons.py

Description :
Ce script génère un fichier de données à partir d'un fichier d'annotation GFF3.
Il calcule, pour chaque fenêtre de 500 kb le long du génome,
la longueur moyenne des exons.

Le script :
- lit le fichier GFF3 (via md.parse_gff3_genes et md.parse_gff3_exons)
- filtre uniquement les entités de type 'gene', 'mRNA' et 'exon'
- ignore les contigs et scaffolds pour ne garder que les chromosomes (ex: Chr01)
- associe chaque exon à son gène parent via les lignes mRNA
- découpe chaque chromosome en fenêtres de 500 kb
- calcule la longueur moyenne des exons dans chaque fenêtre
- écrit un fichier texte tabulé compatible avec Circos

Format de sortie :
Chr01    1       500000    312.45
Chr01    500001  1000000   298.10

Utilisation :
python scripts/script5_lenexons.py input.gff3 output.txt

Exemple :
python scripts/script5_lenexons.py data/Sbicolor_313_v3.1.gene_exons.gff3 results/len_exons.txt
"""

import sys
import module as md

def calculate_mean_len_exons(input_file, output_file):
    """
    Lit un fichier GFF3 et calcule, pour chaque fenêtre de 500 kb,
    la longueur moyenne des exons.

    Paramètres :
    - input_file  : fichier GFF3 en entrée
    - output_file : fichier texte tabulé en sortie
    """
    window_size = 500000

    # Lecture des gènes : md.parse_gff3_genes retourne :
    #   gene_positions : {gene_id: (chrom, start)}
    #   chr_lengths : {chrom: longueur_max}
    gene_positions, chr_lengths = md.parse_gff3_genes(input_file)

    # Lecture des exons : md.parse_gff3_exons retourne :
    #   gene_exon_count : {gene_id: nombre d'exons} (non utilisé ici)
    #   gene_exon_lengths : {gene_id: [longueurs des exons]}
    _, gene_exon_lengths = md.parse_gff3_exons(input_file)

    with open(output_file, 'w') as out:
        for chrom in sorted(chr_lengths.keys()):

            for start, end in md.get_windows(chr_lengths[chrom], window_size):
                genes_in_window = md.get_genes_in_window(gene_positions, chrom, start, end)

                # On ignore les fenêtres sans gènes (bords de chromosome)
                if not genes_in_window:
                    continue

                # On rassemble toutes les longueurs d'exons des gènes de la fenêtre
                # en une seule liste à plat : pour chaque gène de la fenêtre,
                # on récupère sa liste de longueurs d'exons et on les ajoute toutes.
                # Le "if g in gene_exon_lengths" évite une erreur si un gène
                # n'a pas d'exons enregistrés (ex : gène sans transcrit dans le GFF3).
                all_exon_lengths = [l for g in genes_in_window
                                    if g in gene_exon_lengths
                                    for l in gene_exon_lengths[g]]
                if all_exon_lengths:
                    mean_len = sum(all_exon_lengths) / len(all_exon_lengths)
                    out.write(f"{chrom}\t{start + 1}\t{end}\t{mean_len:.2f}\n")

    print(f"Fichier {output_file} créé avec succès.")


def main():
    if len(sys.argv) != 3:
        print("Usage : script5_lenexons.py input.gff3 output.txt")
        sys.exit(1)

    input_file = md.get_input_filename()
    output_file = md.get_output_filename()

    calculate_mean_len_exons(input_file, output_file)


if __name__ == "__main__":
    main()
