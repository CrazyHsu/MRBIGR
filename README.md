# MRBIGR: <ins>M</ins>endelian <ins>R</ins>andomization-<ins>B</ins>ased <ins>I</ins>nference of <ins>G</ins>enetic <ins>R</ins>egulation

MRBIGR is a multifunctional toolkit for pre-GWAS, GWAS and post-GWAS of both traditional and multi-omics data. MRBIGR provides all the components needed to build a complete GWAS pipeline, and integrates with rich post-GWAS analysis tools such as QTL annotation and haplotype analysis. In particular, Mendelian randomization (MR) analysis, MR-based network construction, module identification and gene ontology analysis are proposed for further genetic regulation studies. Additionally, it also produces rich plots for visualization of the analysis results and other formatted data.

## Installation
### Installation using docker (recommended)
```bash
git clone https://github.com/CrazyHsu/MRBIGR.git
cd MRBIGR
docker build -t mrbigr_test .
curr_dir=$(pwd)
docker run -it -d --name mrbigr_env -v $curr_dir/MRBIGR_data:/root/MRBIGR/MRBIGR_data \
    -v $curr_dir/MRBIGR_output:/root/MRBIGR/MRBIGR_output \
    mrbigr_test:latest
```

### Installation using conda
```
git clone https://github.com/CrazyHsu/MRBIGR.git
cd MRBIGR
conda create -n mrbigr python=3.7 -y
conda activate mrbigr
python setup.py build
python setup.py install

pip install pyranges
conda install -y -c conda-forge r-rcppeigen r-xml r-rsqlite r-europepmc r=3.6 rpy2 vcftools
Rscript -e 'install.packages(c("data.table", "ggplot2", "ggsignif", "ggrepel","Matrix", "igraph", "network", "GGally", "sna","tidyr","ggraph","lme4","ggpubr","pheatmap","factoextra", "R.utils"), repos="https://cloud.r-project.org")'
Rscript -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org");BiocManager::install(c("AnnotationForge","clusterProfiler","ggtree"))'
Rscript -e 'install.packages("bigsnpr", dependence=T, repos="https://cloud.r-project.org")'
Rscript -e 'install.packages("https://cran.r-project.org/src/contrib/Archive/FactoMineR/FactoMineR_1.42.tar.gz", repos=NULL)'

echo "export PATH=`pwd`/utils:\$PATH" >> ~/.bashrc
echo "export LD_LIBRARY_PATH=`pwd`/utils/libs:\$LD_LIBRARY_PATH" >> ~/.bashrc
source ~/.bashrc
```

## Usage
To reproduce the images in the MRBIGR paper, you can follow the instruments in [`reproduce/run_from_docker.sh`](reproduce/run_from_docker.sh)

You can also refer to the [online website](https://mrbigr.github.io) for more usage about MRBIGR.


