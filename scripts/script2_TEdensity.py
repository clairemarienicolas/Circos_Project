"""
Script : script2_TEdensity.py

Description :
Ce script génère un fichier de densité des répétitions à partir d’un fichier
FASTA hardmasked. Dans ce type de fichier, les bases masquées sont notées "N".

Le script :
- lit le fichier FASTA
- garde uniquement les chromosomes (et pas les contigs)
- découpe chaque chromosome en fenêtres de 100 kb
- compte le nombre de "N" dans chaque fenêtre
- écrit un fichier texte tabulé compatible avec Circos

Format de sortie :
Chr1    1       100000   10832
Chr1    100001  200000   1683

Utilisation :
python scripts/density_N.py input.fasta output.txt

Exemple :
python scripts/script2_TEdensity.py data/Sbicolor_313_v3.0.hardmasked.fa config/TE_density.txt
"""

import sys
import script1_chr as sc1

def transform(input_file, output_file):
    """
    Lit un fichier FASTA hardmasked et calcule, pour chaque fenêtre de 100 kb,
    le nombre de N présents dans la séquence.

    Paramètres :
    - input_file : fichier FASTA en entrée
    - output_file : fichier texte tabulé en sortie
    """

    chromosomes = {}
    current_chr = None
    window_size = 100000

    # Lecture du fichier FASTA
    with open(input_file, "r") as fh:
        for line in fh:
            line = line.strip()

            # Si c'est un header FASTA
            if line.startswith(">"):
                header = line[1:].split()[0]

                # on garde seulement les chromosomes
                if header.startswith("Chr"):
                    current_chr = header
                    chromosomes[current_chr] = []
                else:
                    current_chr = None

            # Si c'est une ligne de séquence
            elif current_chr is not None:
                chromosomes[current_chr].append(line)

    # Écriture du fichier de sortie
    with open(output_file, "w") as out:
        for chr_name in chromosomes:

            sequence = "".join(chromosomes[chr_name])

            chr_length = len(sequence)

            # découpage en fenêtres de 100 kb
            for start in range(0, chr_length, window_size):
                end = start + window_size
                
                # Si on dépasse la fin du chromosome, on s'arrête à la vraie longueur
                if end > chr_length:
                    end = chr_length

                window_seq = sequence[start:end]
                n_count = window_seq.count("N")
                density = n_count / window_size #pour avoir un ratio entre 0 et 1 au lieu du nombre de bases
                out.write(f"{chr_name}\t{start + 1}\t{end}\t{density:.4f}\n") #start+1 car Circos commence à 1 et pas à 0 comme Python

    print(f"Fichier {output_file} créé avec succès.")

def main():
    if len(sys.argv) != 3:
        print("Usage : script2_TEdensity.py input_filename output_filename")
        sys.exit(1)

    input_file = sc1.get_input_filename()
    output_file = sc1.get_output_filename()

    transform(input_file, output_file)


if __name__ == "__main__":
    main()