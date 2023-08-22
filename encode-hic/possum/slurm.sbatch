#!/usr/bin/bash

#SBATCH --job-name=possum
#SBATCH --output=slurm/%x-%A_%a.out
#SBATCH --partition=cpu
#SBATCH --time=0-00:30:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=60G
#       --array=1-40:1
#SBATCH --array=41-70:1

INPUTS_ROOT=inputs
OUTPUTS_ROOT=outputs

# ---- no more user settings

IFS=',' read INPUT_HIC RESOLUTION OUTPUT_SUBDIR OUTPUT_FILE PARAMETERS <<< $(sed '/^$/d' parameters.csv | awk --assign LN=${SLURM_ARRAY_TASK_ID} 'NR==LN+1{print}')

echo 'INPUT_HIC:' ${INPUT_HIC}
echo 'RESOLUTION:' ${RESOLUTION}
echo 'OUTPUT_SUBDIR:' ${OUTPUT_SUBDIR}
echo 'OUTPUT_FILE:' ${OUTPUT_FILE}
echo 'PARAMETERS:' ${PARAMETERS}

conda activate nk479
module load GCC/7.3.0-2.30

OUTPUT_PATH=${OUTPUTS_ROOT}/$(dirname ${INPUT_HIC})/${OUTPUT_SUBDIR}
mkdir --parents ${OUTPUT_PATH}

OUTPUT_WIG=${OUTPUT_PATH}/${OUTPUT_FILE}.wig
OUTPUT_BIGWIG=${OUTPUT_PATH}/${OUTPUT_FILE}.bigWig
OUTPUT_BEDGRAPH=${OUTPUT_PATH}/${OUTPUT_FILE}.bedGraph

GWevIntra.exe ${PARAMETERS} -T ${SLURM_CPUS_PER_TASK} ${INPUTS_ROOT}/${INPUT_HIC} ${OUTPUT_PATH}/${OUTPUT_FILE}.wig ${RESOLUTION}
wigToBigWig ${OUTPUT_WIG} hg38.chrom.sizes ${OUTPUT_BIGWIG}
bigWigToBedGraph ${OUTPUT_BIGWIG} ${OUTPUT_BEDGRAPH}
