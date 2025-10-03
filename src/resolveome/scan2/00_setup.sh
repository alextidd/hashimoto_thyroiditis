module load ISG/conda
conda create -n scan2 -c conda-forge -c bioconda -c jluquette -c dranew -c soil scan2
conda activate scan2
conda install bioconductor-bsgenome.hsapiens.1000genomes.hs37d5
conda install bioconductor-bsgenome.hsapiens.ucsc.hg38
# add library(rlang) to integrate_tables.R and digest_depth.R
R -e 'install.packages("rlang")'
