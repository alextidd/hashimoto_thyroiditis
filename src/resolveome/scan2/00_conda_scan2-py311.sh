# create env with compatible python and snakemake
module load ISG/conda
conda create -n scan2-py311 python=3.11 -c conda-forge -c bioconda -y
conda activate scan2-py311
conda install -c conda-forge -c bioconda snakemake=9.13.4 -y

# install R + r-scan2 dependencies in same env (conda R or system R; conda recommended)
conda install -c conda-forge r-base r-devtools r-reticulate r-data.table r-dplyr -y

# # point reticulate to this python (export in your job/script)
# export RETICULATE_PYTHON=$(which python)   # this will be the python in scan2-py311