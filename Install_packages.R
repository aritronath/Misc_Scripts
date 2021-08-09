install.packages('Seurat')
install.packages('ggplot2')
install.packages('randomForest')
install.packages('glmnet')
install.packages('GSVA')
install.packages('caret')
install.packages('ggpubr')
install.packages('rms')
install.packages('survminer')
install.packages('data.table')

remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GSVA")
BiocManager::install("GSEABase")
BiocManager::install("biomaRt")
BiocManager::install("impute")
BiocManager::install("ComplexHeatmap")
BiocManager::install("SingleR")
BiocManager::install("celldex")
BiocManager::install("celldex")

# InferCNV and optional (requires JAGS to be installed also) 
install.packages("rjags")
BiocManager::install("infercnv")

# Optional for InferCNV
devtools::install_github("bmbroom/tsvio")
devtools::install_github("bmbroom/NGCHMR", ref="stable")
devtools::install_github("broadinstitute/inferCNV_NGCHM")

