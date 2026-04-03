"""
Script : script3_genedensity.py

Description :
Ce script génère un fichier de densité des gènes à partir d’un fichier
d'annotation GFF3. Il compte le nombre de gènes présents dans chaque fenêtre 
de taille fixe le long du génome.

Le script :
- lit le fichier GFF3
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
import script1_chr as sc1

def transform(input_file, output_file):
    """
    Lit un fichier GFF3 et calcule, pour chaque fenêtre de 100 kb,
    le nombre de gènes présents.

    Paramètres :
    - input_file : fichier GFF3 en entrée
    - output_file : fichier texte tabulé en sortie
    """
    
    gene_positions = {}
    chr_lengths = {}
    window_size = 100000

    # Lecture du GFF pour extraire les gènes
    with open(input_file, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            
            parts = line.split('\t')
            # On filtre : type "gene" et seulement les chromosomes principaux (on ne prend pas en compte les contigs)
            if parts[2] == "gene" and parts[0].startswith("Chr"):
                chrom = parts[0]
                start_gene = int(parts[3])
                end_gene = int(parts[4])
                
                # On utilise le début du gène pour déterminer sa fenêtre
                if chrom not in gene_positions:
                    gene_positions[chrom] = []
                gene_positions[chrom].append(start_gene)
                
                # On suit la longueur max pour définir la fin du chromosome
                if chrom not in chr_lengths or end_gene > chr_lengths[chrom]:
                    chr_lengths[chrom] = end_gene

    # Calcul par fenêtre et écriture du fichier de sortie
    with open(output_file, 'w') as out:
        for chrom in sorted(gene_positions.keys()):
            positions = sorted(gene_positions[chrom])
            max_len = chr_lengths[chrom]
            
            for start in range(0, max_len, window_size):
                end = start + window_size
                if end > max_len:
                    end = max_len
                
                # Compte des gènes dans la fenêtre 
                count = sum(1 for p in positions if start <= p < end)
                
                # Écriture adaptée à Circos (le start commence à 1)
                out.write(f"{chrom}\t{start + 1}\t{end}\t{count}\n")

    print(f"Fichier {output_file} créé avec succès.")

def main():
    if len(sys.argv) != 3:
        print("Usage : script3_genedensity.py input_filename output_filename")
        sys.exit(1)

    input_file = sc1.get_input_filename()
    output_file = sc1.get_output_filename()

    transform(input_file, output_file)

if __name__ == "__main__":
    main()