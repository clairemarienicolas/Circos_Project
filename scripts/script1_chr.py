"""
Script : script1_chr.py

Description :
Ce script permet de générer un fichier karyotype.txt à partir d’un fichier FASTA
contenant un génome (ici le sorgho), qui est compatible avec une utilisation pour Circos.

Le script :
- lit le fichier FASTA
- sélectionne uniquement les chromosomes (ex : Chr01, Chr02…)
- ignore les contigs
- calcule la longueur de chaque chromosome
- écrit un fichier karyotype au format Circos

Format de sortie :
chr - Chr01 Chr01 début fin couleur

Exemple :
chr - Chr01 Chr01 1 80884392 88,114,107

Utilisation :
python scripts/script1_chr.py input.fasta output.txt

Exemple :
python scripts/script1_chr.py data/Sbicolor_313_v3.0.hardmasked.fa results/karyotype.txt
"""

import sys
import module as md


def main():
    """Fonction principale qui gère les arguments et lance la génération du karyotype"""

    if len(sys.argv) != 3:
        print("Usage : script1_chr.py input_filename output_filename")
        sys.exit(1)

    input_file = md.get_input_filename()
    output_file = md.get_output_filename()

    md.read_fasta(input_file, output_file)


if __name__ == "__main__":
    main()
