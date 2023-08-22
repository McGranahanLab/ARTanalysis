#!/usr/bin/bash

#SBATCH --job-name=wget
#SBATCH --output=slurm/%x-%A_%a.out
#SBATCH --partition=cpu
#SBATCH --time='0-00:30:00'
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#       --array=1-13:2
#SBATCH --array=14-20:2

OUTPUT_ROOT=outputs

# download file using core #1
IFS=',' read SAMPLE_ACCESSION OUTPUT_FILE FILE_ACCESSION EXTENSION VERSION <<< $(awk -v LN=${SLURM_ARRAY_TASK_ID} 'NR==LN+1{print}' data_files.csv)
OUTPUT_DOC=${OUTPUT_ROOT}/${OUTPUT_FILE}.${EXTENSION}
SOURCE_URL=https://www.encodeproject.org/files/${FILE_ACCESSION}/@@download/${FILE_ACCESSION}.${EXTENSION}
mkdir --parents `dirname ${OUTPUT_DOC}`
wget --quiet --output-document ${OUTPUT_DOC} ${SOURCE_URL} &

# download file using core #2
if (( `cat data_files.csv | wc -l` > $((${SLURM_ARRAY_TASK_ID}+1)) )); then
	IFS=',' read SAMPLE_ACCESSION OUTPUT_FILE FILE_ACCESSION EXTENSION VERSION <<< $(awk -v LN=${SLURM_ARRAY_TASK_ID} 'NR==LN+2{print}' data_files.csv)
	OUTPUT_DOC=${OUTPUT_ROOT}/${OUTPUT_FILE}.${EXTENSION}
	SOURCE_URL=https://www.encodeproject.org/files/${FILE_ACCESSION}/@@download/${FILE_ACCESSION}.${EXTENSION}
	mkdir --parents `dirname ${OUTPUT_DOC}`
	wget --quiet --output-document ${OUTPUT_DOC} ${SOURCE_URL} &
fi

wait # for the downloads to finish
