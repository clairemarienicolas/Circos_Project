"""
Script : module.py

Description :
Ce script apporte les fonctions utilitaires réutilisables par tous les autres scripts 
du projet (lecture FASTA, lecture GFF3, découpage en fenêtres...).

Les fonctions suivantes sont importées par les autres scripts via :
    import module as md

Fonctions disponibles :
    md.get_input_filename() -> lit sys.argv[1]
    md.get_output_filename() -> lit sys.argv[2]
    md.get_attribute(attrs, key) -> extrait un champ de la colonne 9 GFF3
    md.parse_fasta_sequences(file) -> lit un FASTA, retourne {chr: séquence}
    md.parse_gff3_genes(file) -> lit un GFF3, retourne (gene_positions, chr_lengths)
    md.parse_gff3_exons(file) -> lit un GFF3, retourne (gene_exon_count, gene_exon_lengths)
    md.get_windows(chr_length, window_size) -> génère les fenêtres (start, end)
    md.get_genes_in_window(gene_positions, chrom, win_start, win_end) -> liste des gènes

"""

import sys

# Fonctions de gestion des arguments (input, output)

def get_input_filename():
    """Récupère le nom du fichier en entrée depuis la ligne de commande (sys.argv[1])"""
    return sys.argv[1]


def get_output_filename():
    """Récupère le nom du fichier de sortie depuis la ligne de commande (sys.argv[2])"""
    return sys.argv[2]


# Fonctions GFF3

def get_attribute(attrs, key):
    """
    Extrait la valeur d'un attribut dans la colonne 9 d'un GFF3.

    La colonne 9 ressemble à : ID=Sobic.001G000100.v3.1;Name=Sobic.001G000100;Parent=...
    On découpe sur les ';', puis on cherche le champ qui commence par 'clé='.

    Paramètres :
    - attrs : la chaîne de caractères de la colonne 9
    - key   : le nom de l'attribut à extraire (ex: "ID", "Parent")

    Retourne la valeur trouvée, ou une chaîne vide "" si l'attribut est absent.

    Exemple :
        get_attribute("ID=Sobic.001G000100.v3.1;Name=Sobic.001G000100", "ID")
        -> "Sobic.001G000100.v3.1"
    """
    for attr in attrs.split(';'):
        attr = attr.strip()
        if attr.startswith(key + "="):
            return attr[len(key)+1:]
    return ""


def parse_gff3_genes(input_file):
    """
    Lit un fichier GFF3 et extrait les informations sur les gènes.

    Pour chaque ligne de type 'gene' sur un chromosome principal (Chr...),
    on récupère l'identifiant, le chromosome et la position de début.
    On suit aussi la position maximale de chaque chromosome.

    Paramètres :
    - input_file : chemin vers le fichier GFF3

    Retourne :
    - gene_positions : dict {gene_id: (chrom, start)}
      Ex: {"Sobic.001G000100.v3.1": ("Chr01", 48286)}
    - chr_lengths : dict {chrom: longueur_max}
      Ex: {"Chr01": 80884392}
    """
    gene_positions = {} # {gene_id: (chrom, start)}
    chr_lengths = {} # {chrom: position_max}

    with open(input_file, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue

            parts = line.split('\t')
            if len(parts) < 9:
                continue

            # On ne garde que les chromosomes et on exclus les contigs
            if not parts[0].startswith("Chr"):
                continue

            chrom = parts[0]
            feat_type = parts[2]
            start = int(parts[3])
            end = int(parts[4])
            attrs = parts[8]

            # Mise à jour de la longueur du chromosome
            if chrom not in chr_lengths or end > chr_lengths[chrom]:
                chr_lengths[chrom] = end

            # Récupération des gènes uniquement
            if feat_type == "gene":
                gene_id = get_attribute(attrs, "ID")
                if gene_id:
                    gene_positions[gene_id] = (chrom, start)

    return gene_positions, chr_lengths


def parse_gff3_exons(input_file):
    """
    Lit un fichier GFF3 et extrait les informations sur les exons, par gène.

    On passe par les lignes 'mRNA' pour construire une table de correspondance
    transcrit -> gène, car les exons pointent vers le transcrit et non vers
    le gène directement.

    Paramètres :
    - input_file : chemin vers le fichier GFF3

    Retourne :
    - gene_exon_count : dict {gene_id: nombre d'exons}
      Ex: {"Sobic.001G000100.v3.1": 5}
    - gene_exon_lengths : dict {gene_id: [longueur_exon1, longueur_exon2, ...]}
      Ex: {"Sobic.001G000100.v3.1": [205, 120, 85]}
    """
    transcript_to_gene = {}  # {transcript_id: gene_id}
    gene_exon_count = {}  # {gene_id: nombre d'exons}
    gene_exon_lengths = {}  # {gene_id: [longueurs des exons]}

    with open(input_file, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue

            parts = line.split('\t')
            if len(parts) < 9:
                continue

            if not parts[0].startswith("Chr"):
                continue

            feat_type = parts[2]
            start = int(parts[3])
            end = int(parts[4])
            attrs = parts[8]

            # Les lignes mRNA permettent de relier transcrit -> gène
            if feat_type == "mRNA":
                transcript_id = get_attribute(attrs, "ID")
                gene_id = get_attribute(attrs, "Parent")
                if transcript_id and gene_id:
                    transcript_to_gene[transcript_id] = gene_id
                    # Initialisation des compteurs si premier transcrit du gène
                    if gene_id not in gene_exon_count:
                        gene_exon_count[gene_id]   = 0
                        gene_exon_lengths[gene_id] = []

            # Les lignes exon sont rattachées au gène via la table ci-dessus
            elif feat_type == "exon":
                parent_transcript = get_attribute(attrs, "Parent")
                gene_id = transcript_to_gene.get(parent_transcript, "")
                if gene_id in gene_exon_count:
                    gene_exon_count[gene_id] += 1
                    gene_exon_lengths[gene_id].append(end - start + 1)

    return gene_exon_count, gene_exon_lengths


# Fonctions FASTA

def parse_fasta_sequences(input_file):
    """
    Lit un fichier FASTA et retourne la séquence complète de chaque chromosome.

    Ignore les contigs (seules les entrées dont le header commence par "Chr"
    sont conservées).

    Paramètres :
    - input_file : chemin vers le fichier FASTA

    Retourne :
    - chromosomes : dict {chr_name: séquence_complète}
      Ex: {"Chr01": "ATGCNNATGC..."}
    """
    chromosomes = {} # dictionnaire pour stocker les longueurs des chromosomes
    current_chr = None # variable pour savoir sur quel chromosome on est en train de lire

    # ouverture du fichier FASTA
    with open(input_file, 'r') as fh:
        for line in fh:
            line = line.strip()
            
            # si la ligne commence par ">", c'est un header FASTA
            if line.startswith(">"):
                header = line[1:].split()[0] # enlève ">" et garde le premier mot
                if header.startswith("Chr"): # on garde uniquement les chromosomes (on exclut les contigs)
                    current_chr = header
                    chromosomes[current_chr] = []
                else:
                    current_chr = None

            # si c'est une ligne de séquence et qu'on est sur un chromosome
            elif current_chr is not None:
                chromosomes[current_chr].append(line)

    # On joint les morceaux de séquence en une seule chaîne par chromosome
    return {chr_name: "".join(seq) for chr_name, seq in chromosomes.items()}


def read_fasta(input_file, output_file):
    """
    Génère un fichier karyotype pour Circos à partir d'un fichier FASTA.

    Paramètres :
    - input_file  : chemin vers le fichier FASTA
    - output_file : chemin vers le fichier karyotype à créer
    """
    couleur = "100,149,237"
    chromosomes = parse_fasta_sequences(input_file)

    # écriture du fichier
    with open(output_file, "w") as out:
        for chr_name, sequence in chromosomes.items():
            out.write(f"chr\t-\t{chr_name}\t{chr_name}\t1\t{len(sequence)}\t{couleur}\n")

    print(f"Fichier {output_file} créé avec succès.")


# Fonctions fenêtres

def get_windows(chr_length, window_size):
    """
    Génère toutes les fenêtres (start, end) pour un chromosome donné.

    La dernière fenêtre est ajustée pour ne pas dépasser la longueur du chromosome.
    Les valeurs retournées sont en coordonnées Python (start à 0).
    Pour Circos, il faudra écrire start+1.

    Paramètres :
    - chr_length  : longueur totale du chromosome
    - window_size : taille de chaque fenêtre (ex: 100000 pour 100 kb)

    Retourne une liste de tuples (start, end).

    Exemple :
        get_windows(250000, 100000)
        -> [(0, 100000), (100000, 200000), (200000, 250000)]
    """
    windows = []
    for start in range(0, chr_length, window_size):
        end = min(start + window_size, chr_length)
        windows.append((start, end))
    return windows


def get_genes_in_window(gene_positions, chrom, win_start, win_end):
    """
    Retourne la liste des identifiants de gènes dont le début se trouve
    dans la fenêtre [win_start, win_end[ sur le chromosome donné.

    Paramètres :
    - gene_positions : dict {gene_id: (chrom, start)} (retourné par parse_gff3_genes)
    - chrom          : nom du chromosome (ex: "Chr01")
    - win_start      : début de la fenêtre (inclus)
    - win_end        : fin de la fenêtre (exclus)

    Retourne une liste de gene_id.
    """
    return [
        gene_id for gene_id, (c, s) in gene_positions.items()
        if c == chrom and win_start <= s < win_end
    ]




