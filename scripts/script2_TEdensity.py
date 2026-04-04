"""
Script : script2_TEdensity.py

Description :
Ce script génère un fichier de densité des répétitions à partir d’un fichier
FASTA hardmasked. Dans ce type de fichier, les bases masquées sont notées "N".

Le script :
- lit le fichier FASTA (via md.parse_fasta_sequences)
- garde uniquement les chromosomes (et pas les contigs)
- découpe chaque chromosome en fenêtres de 100 kb
- compte le nombre de "N" dans chaque fenêtre
- écrit un fichier texte tabulé compatible avec Circos

Format de sortie :
Chr1    1       100000   10832
Chr1    100001  200000   1683

Utilisation :
python scripts/script2_TEdensity.py input.fasta output.txt

Exemple :
python scripts/script2_TEdensity.py data/Sbicolor_313_v3.0.hardmasked.fa results/TE_density.txt
"""

import sys
import module as md

def calculate_TE_density(input_file, output_file):
    """
    Lit un fichier FASTA hardmasked et calcule, pour chaque fenêtre de 100 kb,
    le nombre de N présents dans la séquence.

    Paramètres :
    - input_file : fichier FASTA en entrée
    - output_file : fichier texte tabulé en sortie
    """

    window_size = 100000

    # Lecture du FASTA : md.parse_fasta_sequences retourne {chr_name: séquence}
    chromosomes = md.parse_fasta_sequences(input_file)

    with open(output_file, "w") as out:
        for chr_name in sorted(chromosomes.keys()):
            sequence   = chromosomes[chr_name]
            chr_length = len(sequence)

            # Découpage en fenêtres : md.get_windows retourne [(start, end), ...]
            for start, end in md.get_windows(chr_length, window_size):
                window_seq = sequence[start:end]
                n_count    = window_seq.count("N")
                out.write(f"{chr_name}\t{start + 1}\t{end}\t{n_count}\n")

    print(f"Fichier {output_file} créé avec succès.")

def main():
    if len(sys.argv) != 3:
        print("Usage : script2_TEdensity.py input_filename output_filename")
        sys.exit(1)

    input_file = md.get_input_filename()
    output_file = md.get_output_filename()

    calculate_TE_density(input_file, output_file)


if __name__ == "__main__":
    main()