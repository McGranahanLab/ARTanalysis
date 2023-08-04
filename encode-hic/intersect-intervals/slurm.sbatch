#!/usr/bin/bash

#SBATCH --job-name=bedtools
#SBATCH --output=slurm/%x-%A_%a.out
#SBATCH --partition=cpu
#SBATCH --time=0-00:10:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#       --array=1-40:1
#SBATCH --array=41-60:1

IFS=',' read COMPARTMENTS_BEDGRAPH REPLISEQ_BED OUPTUT_FILE <<< $(sed '/^$/d' parameters.csv | awk --assign LN=${SLURM_ARRAY_TASK_ID} 'NR==LN+1{print}')

conda activate nk479

OUPTUT_FILE=outputs/${OUPTUT_FILE}
mkdir --parents `dirname ${OUPTUT_FILE}`

bedtools intersect \
	-a ${REPLISEQ_BED} -wa \
	-b ${COMPARTMENTS_BEDGRAPH} -wb > ${OUPTUT_FILE}
