#! /bin/bash
#$ -S /bin/bash
#$ -M armel.salmon@univ-rennes1.fr
#$ -m bea
#$ -cwd


. /local/env/envperl-5.26.2.sh
. /local/env/envcircos-0.69.8.sh

circos -conf circos.conf
