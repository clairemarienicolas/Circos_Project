"""
Script : script7_len_exons.py

Description :
Ce script génère un fichier de données à partir d'un fichier d'annotation GFF3.
Il calcule, pour chaque fenêtre de 500 kb le long du génome,
la longueur moyenne des exons.

Le script :
- lit le fichier GFF3
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
python scripts/script7_len_exons.py input.gff3 output.txt

Exemple :
python scripts/script7_len_exons.py data/Sbicolor_313_v3.1.gene_exons.gff3 results/len_exons.txt
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
    la longueur moyenne des exons.

    Paramètres :
    - input_file  : fichier GFF3 en entrée
    - output_file : fichier texte tabulé en sortie
    """

    window_size = 500000

    gene_positions = {}  # {gene_id: (chrom, start)}
    gene_exons     = {}  # {gene_id: [longueur_exon1, longueur_exon2, ...]}
    chr_lengths    = {}  # {chrom: position_max}

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
                    gene_positions[gene_id] = (chrom, start)
                    gene_exons[gene_id]     = []

            # Récupération des transcrits pour la table de correspondance
            elif feat_type == "mRNA":
                transcript_id = get_attribute(attrs, "ID")
                gene_id       = get_attribute(attrs, "Parent")
                if transcript_id and gene_id:
                    transcript_to_gene[transcript_id] = gene_id

            # Récupération des longueurs d'exons par gène
            elif feat_type == "exon":
                parent_transcript = get_attribute(attrs, "Parent")
                gene_id = transcript_to_gene.get(parent_transcript, "")
                exon_length = end - start + 1
                if gene_id in gene_exons:
                    gene_exons[gene_id].append(exon_length)

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

                # On rassemble toutes les longueurs d'exons des gènes de la fenêtre
                all_exon_lengths = [l for g in genes_in_window for l in gene_exons[g]]
                if all_exon_lengths:
                    mean_len = sum(all_exon_lengths) / len(all_exon_lengths)
                    out.write(f"{chrom}\t{win_start + 1}\t{win_end}\t{mean_len:.2f}\n")

    print(f"Fichier {output_file} créé avec succès.")


def main():
    if len(sys.argv) != 3:
        print("Usage : script7_len_exons.py input.gff3 output.txt")
        sys.exit(1)

    input_file  = sc1.get_input_filename()
    output_file = sc1.get_output_filename()

    transform(input_file, output_file)


if __name__ == "__main__":
    main()

