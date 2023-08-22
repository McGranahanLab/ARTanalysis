# Identify compartments

This script uses POSSUM to identify compartments in the ENCODE-processed contact matrices for the samples described by `parameters.csv`. For each sample and resolution specified, POSUM will run with default parameters and with `-n SCALE`. The tracks are converted into bedGraph format to intersect with genomic intervals.

The 50kb `-n SCALE` results are used in the manuscript.

```bash
sbatch slurm.sbatch
```
