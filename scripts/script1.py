"""
Script : script1.py

Description :
Ce script permet de générer un fichier karyotype.txt à partir d’un fichier FASTA
contenant un génome (ici le sorgho). Le fichier produit est compatible avec Circos.

Le script :
- lit le fichier FASTA
- sélectionne uniquement les chromosomes (ex : Chr01, Chr02…)
- ignore les contigs/non mappés
- calcule la longueur de chaque chromosome
- écrit un fichier karyotype au format Circos

Format de sortie :
chr - Chr01 Chr01 1 longueur couleur

Exemple :
chr - Chr01 Chr01 1 80884392 88,114,107

Utilisation :
python scripts/script1.py input.fasta output.txt

Exemple concret :
python scripts/script1.py data/Sbicolor_313_v3.0.hardmasked.fa config/karyotype.txt
"""

import sys


def get_input_filename():
    """Récupère le nom du fichier FASTA en entrée depuis la ligne de commande"""
    return sys.argv[1]


def get_output_filename():
    """Récupère le nom du fichier de sortie depuis la ligne de commande"""
    return sys.argv[2]


def transform(input_file, output_file):
    """
    Lit un fichier FASTA et génère un fichier karyotype pour Circos.

    Paramètres :
    - input_file : chemin vers le fichier FASTA
    - output_file : chemin vers le fichier karyotype à créer
    """

    # dictionnaire pour stocker les longueurs des chromosomes
    chromosomes = {}

    # variable pour savoir sur quel chromosome on est en train de lire
    current_chr = None

    # couleur RGB utilisée pour tous les chromosomes
    couleur = "88,114,107"

    # ouverture du fichier FASTA
    with open(input_file, "r") as fh:
        for line in fh:
            line = line.strip()

            # si la ligne commence par ">" → header FASTA
            if line.startswith(">"):
                header = line[1:].split()[0]  # enlève ">" et garde le premier mot

                # on garde uniquement les chromosomes (on exclut les contigs)
                if header.startswith("Chr"):
                    current_chr = header
                    chromosomes[current_chr] = 0
                else:
                    current_chr = None

            # si c'est une ligne de séquence et qu'on est sur un chromosome
            elif current_chr is not None:
                chromosomes[current_chr] += len(line)

    # écriture du fichier karyotype
    with open(output_file, "w") as out:
        for chr_name in chromosomes:
            out.write(
                f"chr - {chr_name} {chr_name} 1 {chromosomes[chr_name]} {couleur}\n"
            )

    print(f"Fichier {output_file} créé avec succès.")


def main():
    """Fonction principale qui gère les arguments et lance le traitement"""

    if len(sys.argv) != 3:
        print("Usage : script1.py input_filename output_filename")
        sys.exit(1)

    input_file = get_input_filename()
    output_file = get_output_filename()

    transform(input_file, output_file)


if __name__ == "__main__":
    main()