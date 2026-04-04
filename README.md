# Projet Circos - Analyse du génome de Sorghum bicolor

Ce projet vise à visualiser les caractéristiques génomiques du sorgho (densité de gènes, répétitions, synténie) à l'aide de l'outil **Circos**.

## 1. Structure du projet

* `Sb_data/` : Données brutes (Fasta, GFF3, fichiers de synténie).
* `scripts/` : Scripts Python de prétraitement des données.
* `results/` : Fichiers formatés générés pour Circos (Karyotype, densités, liens).
* `circos/` : Fichiers de configuration de Circos (.conf).
* `run_all_scripts.sh` : Script de soumission pour le cluster (Slurm/SGE).

## 2. Pré-requis

Le projet a été développé pour fonctionner sur le cluster **GenOuest** avec :
* **Python 3.9.5**
* **Circos 0.69.8**
* Modules Python standards (sys, os) et fichiers utilitaires (`module.py`).

## 3. Prétraitement des données

Pour générer les fichiers nécessaires à Circos, nous utilisons une chaîne de scripts Python. Chaque script a un rôle spécifique :

| Script | Fonction | Sortie |
| :--- | :--- | :--- |
| `script1_chr.py` | Génère le karyotype à partir du FASTA | `karyotype.txt` |
| `script2_TEdensity.py` | Calcule la densité des TE (fenêtres 100kb) | `TE_density.txt` |
| `script3_genedensity.py`| Calcule la densité des gènes (GFF3) | `gene_density.txt` |
| `script4_nbexons.py` | Moyenne d'exons par gène | `nb_exons.txt` |
| `script5_lenexons.py` | Longueur moyenne des exons | `len_exons.txt` |
| `script6_syntenie.py` | Transforme la synténie en liens Circos | `links.txt` |

## 4. Utilisation

### Étape 1 : Régénération des résultats
Pour lancer tous les calculs sur le cluster, utilisez le script de soumission :
```bash
sbatch run_all_scripts.sh