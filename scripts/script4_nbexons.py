"""
Script : script4_nbexons.py

Description :
Ce script génère un fichier de données à partir d'un fichier d'annotation GFF3.
Il calcule, pour chaque fenêtre de 500 kb le long du génome,
le nombre moyen d'exons par gène.

Le script :
- lit le fichier GFF3
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
import script1_chr as sc1

def get_attribute(attrs, key):
    """
    Extrait la valeur d'un attribut dans la colonne 9 d'un GFF3.
    Exemple : get_attribute("ID=Sobic.001G000100.v3.1;Name=...", "ID")
              -> "Sobic.001G000100.v3.1"
    """
    for attr in attrs.split(';'):
        attr = attr.strip()
        if attr.startswith(key + "="):
            return attr[len(key)+1:]
    return ""

def transform(input_file, output_file):
    """
    Lit un fichier GFF3 et calcule, pour chaque fenêtre de 500 kb,
    le nombre moyen d'exons par gène.

    Paramètres :
    - input_file  : fichier GFF3 en entrée
    - output_file : fichier texte tabulé en sortie
    """

    window_size = 500000

    gene_positions  = {}  # {gene_id: (chrom, start)}
    gene_exon_count = {}  # {gene_id: nombre d'exons}
    chr_lengths     = {}  # {chrom: position_max}

    # Table de correspondance transcrit -> gène
    # Nécessaire car les exons pointent vers le transcrit (mRNA), pas vers le gène directement
    transcript_to_gene = {}  # {transcript_id: gene_id}

    # Lecture du GFF3
    with open(input_file, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue

            parts = line.split('\t')
            if len(parts) < 9:
                continue

            # On ne garde que les chromosomes principaux
            if not parts[0].startswith("Chr"):
                continue

            chrom     = parts[0]
            feat_type = parts[2]
            start     = int(parts[3])
            end       = int(parts[4])
            attrs     = parts[8]

            # Mise à jour de la longueur du chromosome
            if chrom not in chr_lengths or end > chr_lengths[chrom]:
                chr_lengths[chrom] = end

            # Récupération des gènes
            if feat_type == "gene":
                gene_id = get_attribute(attrs, "ID")
                if gene_id:
                    gene_positions[gene_id]  = (chrom, start)
                    gene_exon_count[gene_id] = 0

            # Récupération des transcrits pour la table de correspondance
            elif feat_type == "mRNA":
                transcript_id = get_attribute(attrs, "ID")
                gene_id       = get_attribute(attrs, "Parent")
                if transcript_id and gene_id:
                    transcript_to_gene[transcript_id] = gene_id

            # Comptage des exons par gène
            elif feat_type == "exon":
                parent_transcript = get_attribute(attrs, "Parent")
                gene_id = transcript_to_gene.get(parent_transcript, "")
                if gene_id in gene_exon_count:
                    gene_exon_count[gene_id] += 1

    # Calcul par fenêtre et écriture du fichier de sortie
    with open(output_file, 'w') as out:
        for chrom in sorted(chr_lengths.keys()):
            max_len = chr_lengths[chrom]

            for win_start in range(0, max_len, window_size):
                win_end = win_start + window_size
                if win_end > max_len:
                    win_end = max_len

                # Gènes dont le début est dans cette fenêtre
                genes_in_window = [
                    gene_id for gene_id, (c, s) in gene_positions.items()
                    if c == chrom and win_start <= s < win_end
                ]

                # On ignore les fenêtres sans gènes (bords de chromosome)
                if not genes_in_window:
                    continue

                nb_exons_list = [gene_exon_count[g] for g in genes_in_window if gene_exon_count[g] > 0]
                if nb_exons_list:
                    mean_nb_exons = sum(nb_exons_list) / len(nb_exons_list)
                    out.write(f"{chrom}\t{win_start + 1}\t{win_end}\t{mean_nb_exons:.2f}\n")

    print(f"Fichier {output_file} créé avec succès.")


def main():
    if len(sys.argv) != 3:
        print("Usage : script6_nb_exons.py input.gff3 output.txt")
        sys.exit(1)

    input_file  = sc1.get_input_filename()
    output_file = sc1.get_output_filename()

    transform(input_file, output_file)


if __name__ == "__main__":
    main()
