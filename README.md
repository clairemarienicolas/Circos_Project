# Projet Circos - Analyse du génome de Sorghum bicolor

Ce projet vise à visualiser les caractéristiques génomiques du sorgho (*Sorghum bicolor*) à l'aide de l'outil **Circos**. Il permet de représenter simultanément la densité en éléments transposables, la densité en gènes, la complexité des gènes (nombre et longueur des exons) et les blocs de synténie intra-génomique.

---

## 1. Structure du projet

### Dossier Sb_data/ 
Contient toutes les données brutes
* **Génome hardmasked** : Sbicolor_313_v3.0.fa
* **Génome hardmasked** : Sbicolor_313_v3.0.hardmasked.fa
* **Annotation des gènes et exons** : Sbicolor_313_v3.1.gene_exons.gff3
* **Synténie intra-génomique** : Sb_Sb.aligncoords.gcoords.txt 

### Dossier scripts/
Contient tous les scripts permettant de préprocesser les données pour être utilisées par Circos
* **Génération du karyotype** : script1_chr.py 
* **Densité des éléments transposables par fenêtre de 100kb** : script2_TEdensity.py 
* **Densité des gènes par fenêtre de 100kb** : script3_genedensity.py
* **Nombre moyen d'exons par gène par fenêtre de 500kb** : script4_nbexons.py  
* **Longueur moyenne des exons par fenêtre de 500kb** : script5_lenexons.py  
* **Blocs de synténie intragénomique** : script6_syntenie.py 
* **Module contenant toutes les fonctions utilisées dans les autres scripts** : module.py
* **Script Shell d'automatisation pour run tous les scripts facilement** : run_all_scripts.sh

### Dossier results/
Contient tous les résultats des scripts ci-dessus, prêts à être utilisés par Circos
* **Résultat du script 1** : karyotype.txt
* **Résultat du script 2** : TE_density.txt
* **Résultat du script 3** : gene_density.txt
* **Résultat du script 4** : nb_exons.txt
* **Résultat du script 5** : len_exons.txt
* **Résultat du script 6** : links.txt

### Dossier circos/
* **Script de soumission du job Circos** : circos_bash.sh 
* **Fichier principal pour configurer notre Circos (à modifier)** : circos.conf 
* **Configuration de l'idéogramme** : ideogram.conf
* **Labels des chromosomes** : ideogram.label.conf 
* **Position de l'idéogramme** : ideogram.position.conf
* **Graduations** : ticks.conf 
* **Bandes chromosomiques** : bands.conf 

---

## 2. Pré-requis

Le projet est configuré pour fonctionner sur le cluster **GenOuest** avec :

- **Python 3.9.5**
- **Circos 0.69.8**
- **Perl 5.26.2**

---

## 3. Données

| Fichier | Description |
| :--- | :--- |
| `Sbicolor_313_v3.0.fa` | Génome du sorgho |
| `Sbicolor_313_v3.0.hardmasked.fa` | Génome du sorgho hardmasked : les éléments transposables sont remplacés par des "N" |
| `Sbicolor_313_v3.1.gene_exons.gff3` | Annotation des gènes et exons (~30 000 gènes) |
| `Sb_Sb.aligncoords.gcoords.txt` | Coordonnées génomiques des blocs de synténie intra-génomique détectés par DAGChainer |

---

## 4. Prétraitement des données

### 4.1 Module utilitaire (`module.py`)

Toutes les fonctions communes aux scripts sont centralisées dans `module.py`, importé par chaque script via :

```python
import module as md
```

Fonctions disponibles :

| Fonction | Description |
| :--- | :--- |
| `get_input_filename()` | Lit le fichier d'entrée depuis `sys.argv[1]` |
| `get_output_filename()` | Lit le fichier de sortie depuis `sys.argv[2]` |
| `get_attribute(attrs, key)` | Extrait un attribut de la colonne 9 d'un GFF3 |
| `parse_fasta_sequences(file)` | Lit un FASTA, retourne `{chr: séquence}` |
| `read_fasta(input_file, output_file)` | Lit un FASTA, crée le fichier `karyotype.txt` |
| `parse_gff3_genes(file)` | Lit un GFF3, retourne positions des gènes et longueurs des chromosomes |
| `parse_gff3_exons(file)` | Lit un GFF3, retourne comptages et longueurs des exons par gène |
| `get_windows(chr_length, window_size)` | Génère les fenêtres `(start, end)` |
| `get_genes_in_window(gene_positions, chrom, start, end)` | Retourne les gènes d'une fenêtre |

### 4.2 Scripts de prétraitement

| Script | Entrée | Sortie | Description |
| :--- | :--- | :--- | :--- |
| `script1_chr.py` | `.fa` | `karyotype.txt` | Longueur de chaque chromosome pour Circos |
| `script2_TEdensity.py` | `.fa` hardmasked | `TE_density.txt` | Comptage des "N" par fenêtre de 100 kb |
| `script3_genedensity.py` | `.gff3` | `gene_density.txt` | Nombre de gènes par fenêtre de 100 kb |
| `script4_nbexons.py` | `.gff3` | `nb_exons.txt` | Nombre moyen d'exons par gène par fenêtre de 500 kb |
| `script5_lenexons.py` | `.gff3` | `len_exons.txt` | Longueur moyenne des exons par fenêtre de 500 kb |
| `script6_syntenie.py` | `.gcoords.txt` | `links.txt` | Blocs de synténie filtrés (> 200 kb) |

## 5. Utilisation sur le cluster GenOuest

### 5.1 Génération des fichiers résultats

Depuis la racine du projet :

```bash
sbatch scripts/run_all_scripts.sh
```

Ce script soumet un job SLURM qui exécute successivement les 6 scripts Python.

### 5.2 Génération de la figure Circos

Depuis le dossier `circos/` :

```bash
cd circos/
sbatch circos_bash.sh
```

> **Important** : Circos doit être lancé **depuis le dossier `circos/`** et non depuis la racine du projet, car `circos.conf` utilise des chemins relatifs (`../results/`, `etc/`).

Les fichiers générés apparaissent dans `circos/` :

- `circos.png` — image bitmap (3000×3000 px)
- `circos.svg` — image vectorielle

---

## 6. Configuration Circos (`circos.conf`)

Le fichier `circos.conf` est le seul fichier de configuration à modifier. Les autres fichiers (ideogram, ticks, bands) ont été préconfigurés.

### Représentations incluses (de l'extérieur vers l'intérieur)

| Anneau | Type | Fichier | Couleur | Description |
| :--- | :--- | :--- | :--- | :--- |
| 1 (le plus extérieur)| Heatmap | `TE_density.txt` | `orrd-9-seq` (blanc -> rouge) | Densité en éléments transposables |
| 2 | Heatmap | `gene_density.txt` | `blues-9-seq` (blanc -> bleu) | Densité en gènes |
| 3 | Histogramme | `nb_exons.txt` | Bleu acier | Nombre moyen d'exons par gène |
| 4 | Histogramme | `len_exons.txt` | Bleu ciel | Longueur moyenne des exons |
| 5 (le plus intérieur) | Rubans | `links.txt` | Rouge transparent | Blocs de synténie intra-génomique |

> Les palettes de couleurs utilisées (ColorBrewer `orrd`, `blues`) sont adaptées aux personnes daltoniennes.

---

## 7. Résultat

La figure finale représente les 10 chromosomes du sorgho (*Sorghum bicolor* v3.1) avec :

- Une **anticoïncidence visible** entre les régions riches en TE (rouge) et les régions riches en gènes (bleu)
- Les **blocs de synténie intra-génomique** reflétant l'histoire polyploïde du sorgho
- La **complexité structurale des gènes** (nombre et taille des exons) variable selon les régions chromosomiques
