# EnsembleFS: Ensemble feature selection methods for analysis of molecular data
## Description
EnsembleFS is a R package for single feature selection (FS) and ensemble feature selection of molecular data or clinical data (numeric data formats).
This tool is based on several feature filters, such as the test Manna-Whitneya (U-test), the Monte Carlo Feature Selection (MCFS) ([Dramiński & Koronacki 2018](https://www.jstatsoft.org/article/view/v085i12)) , the MultiDimensional Feature Selection (two variants: MDFS-1D and MDFS-2D) ([Mnich & Rudnicki 2020](https://www.sciencedirect.com/science/article/abs/pii/S0020025520302048)), and the Minimum Redundancy Maximum Relevance (MRMR) 
([Ding 2005](https://pubmed.ncbi.nlm.nih.gov/15852500/)) for discovering the most important biomarkers and used the machine learning algorithms (ML) to evaluate the quality of feature sets. Predictive models are built using the Random Forest algorithm ([Breiman 2001](https://link.springer.com/article/10.1023/A:1010933404324)). It can be applied to two-class problems.

Moreover, EnsembleFS support users in analysis and interpretation of molecular data. The information about each of top biomarkers is extracted from diverse biological databases, namely the Gene Ontology ([GO](https://pubmed.ncbi.nlm.nih.gov/33290552/)), the Kyoto Encyclopedia of Genes and Genomes ([KEGG](https://pubmed.ncbi.nlm.nih.gov/18477636/)), the Reactome ([React](https://pubmed.ncbi.nlm.nih.gov/32907876/)), the WikiPathways ([WP](https://pubmed.ncbi.nlm.nih.gov/33211851/)), the Transfac ([TF](https://pubmed.ncbi.nlm.nih.gov/8594589/)), the miRTarBase ([MIRNA](https://academic.oup.com/nar/article/48/D1/D148/5606625)), the Human Protein Atlas ([HPA](https://pubmed.ncbi.nlm.nih.gov/25613900/)), the [CORUM](https://pubmed.ncbi.nlm.nih.gov/30357367/), and the Human Phenotype Ontology ([HPO](https://pubmed.ncbi.nlm.nih.gov/33264411/)).
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
The RNA-sequencing data of tumor-adjacent normal tissues of lung adenocarcinoma cancer patients from The Cancer Genome Atlas database ([TCGA](https://www.cancer.gov/tcga)) was used. The preprocessing of data involved standard steps for RNA-Seq data. The log2 transformation was performed. Features with zero and near-zero (1%) variance across patients were removed. After the preprocessing procedure the primary dataset contains 574 samples (59 normal and 515 tumor) described with 20172 differentially expressed genes (DEGs). This dataset includes highly correlated features and the number of cancer samples is roughly ten times more than normal samples. For testing purposes, the number of molecular markers was limited to random 2000 DEGs ([exampleData.csv](https://github.com/biocsuwb/EnsembleFS-package/tree/main/data)) and 5000 DEGs ranked by
highest difference in the gene expression level between tumor and normal tissues ([exampleData_5000.csv](https://github.com/biocsuwb/EnsembleFS-package/tree/main/data)). 

## Install the development version from GitHub:

```r
# install.packages("devtools")
devtools::install_github("biocsuwb/EnsembleFS-package")
```
## Notes: 
To install EnsembleFS package in your R environment make sure you have Java installed (rJava R package).

To accelerate processing by using a CUDA GPU the EnsembleFS package must be compiled with CUDA (set the MDFS-2D parameter: use.cuda = TRUE) 

## Example 1 - individual feature selection

#### Loading data
```r
data <- read.csv2('exampleData.csv')
class <- data$class
data$class <- NULL
```

#### Model configuration parameters
- U-test and MDFS-1D, and MDFS-2D parameter, ***multitest correction: adjust = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")***;
- U-test and MDFS-1D, and MDFS-2D parameter, significance level: ***alpha = 0.05***;
- MDFS-2D parameter: use.cuda = FALSE 
- whether to use CUDA acceleration (must be compiled with CUDA) for MDFS-2D method
- MCFS parameter, cut-off method: ***cutoff.method = c("permutations", "criticalAngle", "kmeans")***;
- correlation coefficient: ***level.cor = 0.75***;
- validation methods: ***method.cv = c('kfoldcv','rsampling')***;
- number of repetitions: ***niter = 5***;
- train-test-split the data: ***k = 3***.


#### Feature selection U-test
```r
var.utest <- fs.utest(x = data, y = class, params = list(adjust = "holm", alpha = 0.05))
```
#### Feature selection MDFS-1D
```r
var.mdfs1D <- fs.mdfs.1D(x = data, y = class, params = list(adjust = "holm", alpha = 0.05))
```
#### Feature selection MDFS-2D
```r
var.mdfs2D <- fs.mdfs.2D(x = data, y = class, params = list(adjust = "holm", alpha = 0.05, use.cuda = FALSE))
```
#### Feature selection MCFS
```r
var.mcfs <- fs.mcfs(x = data, y = class,  params = list(cutoff.method = "kmeans"))
```
#### Feature selection MRMR
```r
var.mrmr <- fs.mrmr(x = data, y = class,  params = list(feature.number = 100))
```

#### Creating the cross-validation index array  (3-fold cross-validation repeated 10 times)
```r
list.index.cross <- cross.val(x = data,
                              y = class,
                              method = 'kfoldcv',
                              params.cv = list(niter = 10, k = 3)
```                              
#### Feature selection in a cross-validation scenario for individual FS method (eg. MDFS-1D)
```r
list.selected.var <- feature.selection.cv(x = data,
                                          y = class,
                                          method = 'fs.mdfs.1D',
                                          list.index.cross = list.index.cross,
                                          params = list(adjust = 'holm', alpha = 0.05)
 
# show 10 top features from the first list of most informative features (list 1 of 30)
print(list.selected.var[[1]][1:10,])

       name       Pvalue adjustPvalue
834   ADPRH 2.780479e-26 5.560957e-23
117   SRPK1 1.919463e-23 3.837007e-20
680  RASIP1 4.470157e-23 8.931374e-20
663   PDIA4 7.100476e-23 1.417965e-19
353   PRKCE 8.564362e-23 1.709447e-19
119  NCKAP5 1.785219e-22 3.561512e-19
464    FGF2 3.127923e-22 6.237079e-19
538   IL3RA 5.110742e-22 1.018571e-18
636    EMP2 7.996179e-22 1.592839e-18
1952  TRPV2 3.432695e-21 6.834496e-18

```

#### Building and testing ML models on top N = 100 features with the individual FS algorithm (eg. MDFS-1D); 30 models in total.
```r
model.result <- build.model.crossval(x = data,
                                     y = class,
                                     list.selected.var = list.selected.var,
                                     list.index.cross = list.index.cross,
                                     nvar = 100)
# show a predictive model result (list 1 of 30)
print(model.result[[1]]) 

 Accuracy       AUC       MCC 
0.9807692 0.9736842 0.9589080 
```
#### Computing ASM value for top 100 features
```r
asm <- stabilty.selection(list.selected.var, list.index.cross,  nvar = 100)
# show ASM value 
print(asm)
0.6915541
```
## Example 2 - ensemble feature selection

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
- correlation coefficient: ***level.cor = 1***;
- validation methods: ***method.cv = c('kfoldcv','rsampling')***;
- number of repetitions: ***niter = 5***;
- train-test-split the data: ***k = 3***.

#### Runing end-to-end EnsembleFS for ensemble feature selection and classification model construction 
Number of top features: N = 5, 10, 15, 20, ..., 50, 75, 100

Model validation technique: 3-fold cross-validation repeated 10 times

Selected feature filters: U-test, MCFS, MRMR, and MDFS-1D.

The cut off value of the Spearman correlation coefficient: 0.75
```r
result <- ensembleFS(x = data,
                     y = class,
                     methods = c("fs.utest", "fs.mcfs", "fs.mrmr", "fs.mdfs.1D"),
                     method.cv = "kfoldcv",
                     params.cv = list(k = 3, niter = 5),
                     level.cor = 0.75,
                     params = list(adjust = "holm", cutoff.method = "kmeans", feature.number = 100, alpha = 0.05),
                     asm = c("fs.utest", "fs.mcfs", "fs.mrmr", "fs.mdfs.1D"),
                     model = c("fs.utest", "fs.mcfs", "fs.mrmr", "fs.mdfs.1D"))
 ```
                     
#### Visualizing the model results
```r
graph.result(result$model, "acc")
graph.result(result$stability, "stability")
# graph.result(result$model, "auc")
# graph.result(result$model, "mcc")
```
![Fig.4](https://github.com/biocsuwb/Images/blob/main/ASM&ACC.png?raw=true)
Fig.4 The average values for accurancy (ACC) vs N top features (5, 10, 15, 20, ..., 50, 75, 100) for various features filters and the ASM similarity measure between m = 15 feature subsets vs N top features.
#### Showing m-list of top biomarkers for each of filter FS methods
```r
result$selected.feature
```

#### Showing the combined list of top biomarkers.
How many times a biomarker has occurred in m feature subsets: level.freq = 7

Number of top N biomarkers for each of filter FS methods: number.gene = 100
```r
gene.top <- get.top.gene(list.imp.var.cv = result$selected.feature, level.freq = 7, number.gene = 100)
```

#### Getting information about biomarkers from nine biological databases

Combination of a set of biomarkers: union
```r
info.gene <- get.info.top.gene(gene.top, condition.methods = 'union')
```
## Example 3 - create and add their own feature filters to default list of basic feature filters 
```r
data <- read.csv2('exampleData.csv')
class <- data$class
data$class <- NULL
```
####   Installing the required package
```r
install.packages('mt')
library(mt)
```
####  Creating new feature filter eg. ReliefF ([Kononenko 1994](https://link.springer.com/chapter/10.1007/3-540-57868-4_57))
```r
feature.number = 100
```
#### Rules for adding a new FS method to the benchmark procedure:
1. the function name of the new FS method should have the prefix 'fs.';
2. the created function should take as input arguments:
- x tabular input, numeric type, where columns are variables, and the lines are observations;
- y decision variable as a binary vector of length equal to a number observation;
- params is a list of hyperparameters of the new FS method added, if such are required, and in their absence, the params argument of the above-mentioned
  functions should be omitted;
3. function should return a data frame that consists of two columns : name (names of relevant biomarkers) and score (validity metric
variables for an individual FS method, e.g. p-value for U-test); 
Example:
```r
fs.relieff <- function(x, y, params = list(feature.number = 100)){
  result <- fs.relief(x, y)
  stats <- as.data.frame(result$stats)
  var.names <- row.names(stats)
  scores <- stats[,1]
  var.imp <- as.data.frame(cbind(var.names, scores))
  names(var.imp) <- c('name', 'score')
  var.imp <- var.imp[order(var.imp$score, decreasing=T),][1:params$feature.number,]
  return(var.imp)
}
```
#### Run end-to-end EnsembleFS for ensemble feature selection and comparison of feature filters 
Feature filters: U-test, MCFS, MRMR, MDFS-1D, and ReliefF.
```r
result2 <- ensembleFS(x = data,
                     y = class,
                     methods = c("fs.utest", 'fs.mrmr', 'fs.mcfs' , 'fs.mdfs.1D' ,'fs.relieff'),
                     method.cv = "kfoldcv",
                     params.cv = list(k = 3, niter = 5),
                     level.cor = 0.75,
                     params = list(adjust = "fdr", feature.number = 100, alpha = 0.05, use.cuda = FALSE, cutoff.method = 'kmeans'),
                     asm = c("fs.utest", 'fs.mrmr' , 'fs.mcfs' , 'fs.mdfs.1D' ,'fs.relieff'),
                     model = c("fs.utest",  'fs.mrmr',  'fs.mcfs' , 'fs.mdfs.1D' ,'fs.relieff'))
 ```
 #### Visualizing the model results
```r
graph.result2(result2$model, "acc")
```
![ensembleFS_relieff_auc](https://user-images.githubusercontent.com/36896714/201480835-8af6ba21-10a0-45ef-86c2-353711520d4b.png)
