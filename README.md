# EnsembleFS: Ensemble feature selection methods for analysis of molecular data
## Description
EnsembleFS is a R package for single and ensemble feature selection (FS) of molecular data or clinical data (numeric data formats).
This tool is based on several feature filters, such as the test Manna-Whitneya (U-test), the Monte Carlo Feature Selection (MCFS), the MultiDimensional Feature Selection (two variants: MDFS-1D and MDFS-2D), and the Minimum Redundancy Maximum Relevance (MRMR) for discovering the most important biomarkers and used the machine learning algorithms to evaluate the quality of feature sets. Predictive models are built using the Random Forest algorithm. It can be applied to two-class problems.

Moreover, EnsembleFS support users in analysis and interpretation of molecular data. The information about each of top biomarkers is extracted from diverse biological databases, namely the Gene Ontology, the Kyoto Encyclopedia of Genes and Genomes, the Reactome, the WikiPathways, the Transfac, the miRNA targets, the Human Protein Atlas, the CORUM, and the Human Phenotype Ontology.
The proposed tool accept molecular data includes different types of gene identifiers, such as Ensembl, NCBI Entrez gene ID, Refseq, Illumina, and Uniprot.

EnsembleFS allows the user to:
- filter the most informative features (biomarkers) by using up to five FS methods from molecular data generated from high-throughput molecular biology experiments;
- filter the most informative features by using up to five FS methods from clinical data (numeric data);
- add any other FS methods to default list of basic filters (U-test, MCFS, MRMR, MDFS-1D, and MDFS-2D);
- establish the selected parameters for predictive models, such as the number of top N informative features;
- remove redundant features by building a the Spearman correlation matrix that identifies highly correlated features;
- evaluate the stability of feature subsets and performance of predictive models;
- find information about selected molecular markers (gene ontology, pathways, tissue specificity, miRNA targets, regulatory motif, protein complexes, disease phenotypes) in nine biological databases.

## Install the development version from GitHub:

```r
# install.packages("devtools")
devtools::install_github("biocsuwb/EnsembleFS-package")
```
#### Note: To install EnsembleFS package in your R environment make sure you have Java installed (rJava R package).
## Examples 

EnsembleFS makes it trivial to run many algorithms and use the best one or an ensemble.

```r

data <- read.csv2('examplesDataTest.csv')
class <- data$class
data$class <- NULL

# showing list methods
list.methods()

# run end-to-end Ensemble for comparison of feature selection methods.
result <- ensembleFS(x = data,
                     y = class,
                     methods = c("fs.utest", "fs.mcfs", "fs.mrmr", "fs.mdfs.1D", "fs.mdfs.2D"),
                     method.cv = "kfoldcv",
                     params.cv = list(k = 3, niter = 5),
                     level.cor = 1,
                     params = list(adjust = "holm", feature.number = 10, alpha = 0.05),
                     asm = c("fs.utest", "fs.mcfs", "fs.mrmr", "fs.mdfs.1D", "fs.mdfs.2D"),
                     model = c("fs.utest", "fs.mcfs", "fs.mrmr", "fs.mdfs.1D", "fs.mdfs.2D"))
                     

# showing result ensemble
graph.result(result$stability, "stability")
graph.result(result$model, "auc")

# getting information about genes from gprogiler2
gene.top <- get.top.gene(result$selected.feature, 15 , 20)
info.gene <- get.info.top.gene(gene.top, condition.methods = 'union')

# feature selection U-test
var.utest <- fs.utest(x = data, y = class, params = list(adjust = "holm", alpha = 0.05))

# feature selection MCFS-ID
var.mcfs <- fs.mcfs(x = data, y = class)

#create subset indexes cross-validation
list.index.cross <- cross.val(x = data,
                              y = class,
                              method = 'kfoldcv',
                              params.cv = list(niter = 10, k = 3)
                              
# feature selection for one method in cross-validation
list.selected.var <- feature.selection.cv(x = data,
                                          y = class,
                                          method = 'fs.mdfs.2D',
                                          list.index.cross = list.index.cross,
                                          params = list(adjust = 'holm', alpha = 0.05)
 
# Compute Lustgarten’s stability measure for one method
asm <- stabilty.selection(list.selected.var, list.index.cross, 100)

# Train model Random Forest for one method
model.result <- build.model.crossval(x = data,
                                     y = class,
                                     list.selected.var = list.selected.var,
                                     list.index.cross = list.index.cross,
                                     nvar = 100)



```
