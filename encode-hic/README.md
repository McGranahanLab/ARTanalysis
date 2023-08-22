# Installing software

## Conda environment

A conda (4.12.0) environment is first created to keep all software for the project:

```bash
conda create --name nk479 \
bedtools \
ucsc-bigwigsummary \
ucsc-bigwiginfo \
ucsc-bigwigtowig \
ucsc-fetchchromsizes \
ucsc-wigtobigwig \
ucsc-bigwigtobedgraph

conda activate nk419
```
## POSSUM

POSSUM was used to define compartments, and is installed into the conda environment:

```bash
github clone git@github.com:moshe-olshansky/EigenVector.git #a1df194
github clone git@github.com:aidenlab/straw.git #2525edc

ln --symbolic --relative straw EigenVector/C++/PowerMethod/STRAW

module load GCC/7.3.0-2.30

cd EigenVector/C++/PowerMethod

BIN=${CONDA_PREFIX}/bin
g++ -O2 -o ${BIN}/ev.exe theEigenVector_flip_new.cpp theBestEigen.c thdMul.c hgFlipSign.c STRAW/C++/straw.cpp -I . -I STRAW/C++ -lz -lcurl -lpthread
g++ -O2 -o ${BIN}/GWevIntra.exe GWevIntra_new.cpp theBestEigen.c thdMul.c hgFlipSign.c STRAW/C++/straw.cpp -I . -I STRAW/C++ -lz -lcurl -lpthread
```
