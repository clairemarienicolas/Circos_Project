#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -M claire-marie.nicolas@etudiant.univ-rennes.fr
#$ -m bea

. /local/env/envperl-5.26.2.sh
. /local/env/envcircos-0.69.8.sh

circos -conf circos.conf
