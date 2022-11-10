# EnsembleFS: Ensemble feature selection methods for analysis of molecular data
## Description
EnsembleFS is a R package for single feature selection (FS) and ensemble feature selection of molecular data or clinical data (numeric data formats).
This tool is based on several feature filters, such as the test Manna-Whitneya (U-test), the Monte Carlo Feature Selection (MCFS) ([Dramiński & Koronacki 2018](https://www.jstatsoft.org/article/view/v085i12)) , the MultiDimensional Feature Selection (two variants: MDFS-1D and MDFS-2D) ([Mnich & Rudnicki 2020](https://www.sciencedirect.com/science/article/abs/pii/S0020025520302048)), and the Minimum Redundancy Maximum Relevance (MRMR) 
([Ding 2005](https://pubmed.ncbi.nlm.nih.gov/15852500/)) for discovering the most important biomarkers and used the machine learning algorithms to evaluate the quality of feature sets. Predictive models are built using the Random Forest algorithm ([Breiman 2001](https://link.springer.com/article/10.1023/A:1010933404324)). It can be applied to two-class problems.

Moreover, EnsembleFS support users in analysis and interpretation of molecular data. The information about each of top biomarkers is extracted from diverse biological databases, namely the Gene Ontology [GO](https://pubmed.ncbi.nlm.nih.gov/33290552/), the Kyoto Encyclopedia of Genes and Genomes ([KEGG](https://pubmed.ncbi.nlm.nih.gov/18477636/)), the Reactome ([React](https://pubmed.ncbi.nlm.nih.gov/32907876/)), the WikiPathways ([WP](https://pubmed.ncbi.nlm.nih.gov/33211851/)), the Transfac ([TF](https://pubmed.ncbi.nlm.nih.gov/8594589/)), the miRTarBase ([MIRNA](https://academic.oup.com/nar/article/48/D1/D148/5606625)), the Human Protein Atlas ([HPA](https://pubmed.ncbi.nlm.nih.gov/25613900/)), the [CORUM](https://pubmed.ncbi.nlm.nih.gov/30357367/), and the Human Phenotype Ontology ([HPO](https://pubmed.ncbi.nlm.nih.gov/33264411/)).
The proposed tool accept molecular data includes different types of gene identifiers, such as Ensembl, NCBI Entrez gene ID, Refseq, Illumina, and Uniprot.

EnsembleFS allows the user to:
- filter the most informative features by using up to five FS methods (U-test, MCFS, MRMR, MDFS-1D, and MDFS-2D) from molecular data generated from high-throughput molecular biology experiments
and also clinical data (numeric formats);
- create and add their own feature filters to default list of basic feature filters (U-test, MCFS, MRMR, MDFS-1D, and MDFS-2D) for more complex FS task;
- build a Random Forest classifiers using selected top N features (Fig.1); 
- evaluate the stability of feature subsets (the Lustgarten adjusted stability measure ([ASM](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2815476/)) ) and the performance of predictive models (the area under receiver operator curve  (AUC), accuracy (ACC), and the Matthews correlation coefficient(MCC);
- compare the predictive performance of models and the stability of selected feature sets for selected FS algorithms; 
- establish the selected parameters for predictive models, such as the number of top N informative features;
- remove redundant features by building a the Spearman correlation matrix that identifies highly correlated features;
- product plots to visualize the model results;
- create detailed raport with feature selection and modeling results;
- find information about selected top molecular markers (gene ontology, pathways, tissue specificity, miRNA targets, regulatory motif, protein complexes, disease phenotypes) in nine biological databases (GO, KEGG, React, WP, TF, MIRNA, HPA, CORUM, and HPO).

![Fig.1](https://github.com/biocsuwb/Images/blob/main/Scheme1.png?raw=true)
Fig.1 The scheme of ensemble feature selection and supervised classification.
![Fig.2](https://github.com/biocsuwb/Images/blob/main/Scheme2.png?raw=true)
Fig.2 The scheme of construction of the combined set of N-top relevant features.
![Fig.3](https://github.com/biocsuwb/Images/blob/main/Scheme3.png?raw=true)
Fig.3 The scheme for biological information collection and integration about biomarkers.

## Example data sets
The RNA-sequencing data of tumor-adjacent normal tissues of lung adenocarcinoma cancer patients from The Cancer Genome Atlas database ([TCGA](https://www.cancer.gov/tcga)) was used. The preprocessing of data involved standard steps for RNA-Seq data. The log2 transformation was performed. Features with zero and near-zero (1%) variance across patients were removed. After the preprocessing procedure the primary dataset contains 574 samples (59 normal and 515 tumor) described with 20172 differentially expressed genes (DEGs). This dataset includes highly correlated features and the number of cancer samples is roughly ten times more than normal samples. For testing purposes, the number of molecular markers was limited to random 2000 DEGs ([exampleData.csv](https://github.com/biocsuwb/EnsembleFS-package/tree/main/data)) and 500 DEGs with the
highest difference in the gene expression level between tumor and normal tissues ([exampleData_500.csv](https://github.com/biocsuwb/EnsembleFS-package/tree/main/data)). 

## Install the development version from GitHub:

```r
# install.packages("devtools")
devtools::install_github("biocsuwb/EnsembleFS-package")
```
#### Note: To install EnsembleFS package in your R environment make sure you have Java installed (rJava R package).
## Examples 

#### Loading data
```r
data <- read.csv2('exampleData.csv')
class <- data$class
data$class <- NULL
```

#### Showing the list of available feature selection methods
```r
list.methods()
```

#### Model configuration parameters
EnsembleFS allows user to set some parameter values, such as:
- feature selection methods: ***methods = c("fs.utest", "fs.mcfs", "fs.mrmr", "fs.mdfs.1D", "fs.mdfs.2D")***;
- U-test and MDFS parameter, ***multitest correction: adjust = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")***;
- U-test and MDFS parameter, significance level: ***alpha = 0.05***;
- MRMR parameter, number of significant features: ***feature.number = 100***;
- MCFS parameter, cut-off method: ***cutoff.method = c("permutations", "criticalAngle", "kmeans")***;
- correlation coefficient: ***level.cor = 0.75***;
- validation methods: ***method.cv = c('kfoldcv','rsampling')***;
- number of repetitions: ***niter = 5***;
- train-test-split the data: ***k = 3***.



#### Run end-to-end EnsembleFS for ensemble feature selection and comparison of feature filters
Selected feature filters: U-test, MCFS, MRMR, and MDFS-1D.
```r
result <- ensembleFS(x = data,
                     y = class,
                     methods = c("fs.utest", "fs.mcfs", "fs.mrmr", "fs.mdfs.1D"),
                     method.cv = "kfoldcv",
                     params.cv = list(k = 3, niter = 5),
                     level.cor = 1,
                     params = list(adjust = "holm", cutoff.method = "kmeans", feature.number = 10, alpha = 0.05),
                     asm = c("fs.utest", "fs.mcfs", "fs.mrmr", "fs.mdfs.1D"),
                     model = c("fs.utest", "fs.mcfs", "fs.mrmr", "fs.mdfs.1D"))
 ```
                     
#### Visualizing the model results;
```r
graph.result(result$stability, "stability")
graph.result(result$model, "auc")
```
![Fig.4](https://github.com/biocsuwb/Images/blob/main/ASM&ACC.png?raw=true)
Fig.4 The average values for accurancy (ACC) vs N top features for various features filters and the ASM similarity measure between 15 feature subsets vs N top features.
#### Getting information about biomarkers from databases:
#### the Gene Ontology, the KEGG, the Reactome, the WikiPathways, the Transfac, the miRTarBase, the Human Protein Atlas, the CORUM, and the Human Phenotype Ontology.
```r
gene.top <- get.top.gene(result$selected.feature, 15 , 20)
info.gene <- get.info.top.gene(gene.top, condition.methods = 'union')
```
#### Feature selection U-test
```r
var.utest <- fs.utest(x = data, y = class, params = list(adjust = "holm", alpha = 0.05))
```

#### Feature selection MCFS
```r
var.mcfs <- fs.mcfs(x = data, y = class)
```

#### Create subset indexes cross-validation
```r
list.index.cross <- cross.val(x = data,
                              y = class,
                              method = 'kfoldcv',
                              params.cv = list(niter = 10, k = 3)
```                              
#### Feature selection for one method in cross-validation
```r
list.selected.var <- feature.selection.cv(x = data,
                                          y = class,
                                          method = 'fs.mdfs.2D',
                                          list.index.cross = list.index.cross,
                                          params = list(adjust = 'holm', alpha = 0.05)
 ```
#### Compute Lustgarten’s stability measure for one method
```r
asm <- stabilty.selection(list.selected.var, list.index.cross, 100)
```

#### Train model Random Forest for one method
```r
model.result <- build.model.crossval(x = data,
                                     y = class,
                                     list.selected.var = list.selected.var,
                                     list.index.cross = list.index.cross,
                                     nvar = 100)
```
