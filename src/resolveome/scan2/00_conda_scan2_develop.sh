#!/bin/bash

# modules
module load ISG/conda

# build conda env
conda create -v -n scan2_dev -c conda-forge -c bioconda -c dranew -c soil \
    snakemake-executor-plugin-slurm \
    gatk4 \
    snakemake python sigprofilermatrixgenerator shyaml \
    shapeit eagle-phase \
    'bcftools>=1.19' bedtools samtools pysam htslib \
    bioconductor-bsgenome bioconductor-genomeinfodb bioconductor-genomeinfodbdata bioconductor-rsamtools bioconductor-genomicranges \
    r-fastghquad r-argparse r-r.utils r-bedtoolsr r-future r-future.apply r-viridislite r-viridis \
    r-devtools r-yaml r-pracma r-progressr r-rhpcblasctl r-reticulate r-qs r-digest r-upsetr r-kernsmooth \
    bioconductor-annotationdbi bioconductor-annotationfilter bioconductor-biovizbase \
    bioconductor-bsgenome.hsapiens.1000genomes.hs37d5 \
    bioconductor-bsgenome.hsapiens.ucsc.hg19 bioconductor-bsgenome.hsapiens.ucsc.hg38 \
    bioconductor-genomicalignments bioconductor-genomicfeatures \
    bioconductor-gviz bioconductor-txdb.hsapiens.ucsc.hg19.knowngene \
    bioconductor-variantannotation bioconductor-dnacopy bioconductor-ctc \
    bx-python pandas \
    r-rcpp r-rlang r-curl r-harmonicmeanp r-r.matlab \
    r-lme4 r-lmertest r-nnls r-openxlsx r-pheatmap r-rcolorbrewer r-readxl \
    r-inline r-fastcluster r-heatmap3 r-rio r-mclust \
    r-ggpmisc r-ggplot2 r-forcats r-ggsignif r-patchwork r-broom r-ggseqlogo r-geomtextpath r-ggrepel r-gghighlight \
    r-ggforce r-concaveman r-umap r-scales r-ggbreak r-ggh4x \
    s3transfer \
    ucsc-bedgraphtobigwig ucsc-bigwigaverageoverbed \
    vcf2maf vcflib bigtools wiggletools \
    sigprofilerassignment

# add generic cluster plugin
conda install -y -c bioconda -c conda-forge snakemake-executor-plugin-cluster-generic

# r libraries
export LD_LIBRARY_PATH=/software/conda/users/at31/scan2_dev/lib:$LD_LIBRARY_PATH
/software/conda/users/at31/scan2_dev/bin/Rscript -e "devtools::install_github('AlexandrovLab/SigProfilerMatrixGeneratorR')"
/software/conda/users/at31/scan2_dev/bin/Rscript -e "devtools::install_github('parklab/r-scan2')"
/software/conda/users/at31/scan2_dev/bin/Rscript -e "install.packages('stringfish', repos='https://cloud.r-project.org')"

# sigprofiler data - GRCh37/GRCh38
python -c "from SigProfilerMatrixGenerator import install as genInstall; genInstall.install('GRCh37', rsync=False, bash=True)"
python -c "from SigProfilerMatrixGenerator import install as genInstall; genInstall.install('GRCh38', rsync=False, bash=True)"
