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
- ***to install EnsembleFS package in your R environment make sure you have Java installed (rJava R package);***
- ***to accelerate processing by using a CUDA GPU the EnsembleFS package must be compiled with CUDA (set the MDFS-2D parameter: use.cuda = TRUE).*** 

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
- MDFS-2D parameter: ***use.cuda = FALSE*** 
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
... etc.
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
print(list.methods())
"fs.mcfs"  "fs.mdfs.1D"  "fs.mdfs.2D"  "fs.mrmr"  "fs.utest"  
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

#### Building and testing ML models on top N features with each of selected feature filters (eg. U-test, MCFS, MRMR, and MDFS-1D); 
- number of top features: ***N = 5, 10, 15, 20, ..., 50, 75, 100;***
- model validation technique: ***3-fold cross-validation repeated 5 times;***
- selected feature filters: ***U-test, MCFS, MRMR, and MDFS-1D;***
- the cut off value of the Spearman correlation coefficient: ***0.75.***
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
#### Showing m-list of top biomarkers for each of feature filters (m = 3 x 5 in this case)
```r
print(result$selected.feature)

$fs.utest[[1]]
             name       Pvalue adjustPvalue
1            EMP2 1.440258e-17 2.880515e-14
88        ONECUT1 2.889262e-14 5.527158e-11
89       BAIAP2L2 3.195846e-14 6.110458e-11
99          CILP2 5.888033e-14 1.119904e-10
108         FOLR3 1.097770e-13 2.078079e-10
115          PVT1 1.317063e-13 2.486614e-10
120      MGC42105 1.607918e-13 3.026102e-10
123         LYPD1 1.689979e-13 3.173781e-10
125           XDH 1.866637e-13 3.501812e-10  
... etc.
```

#### Showing the combined list of top biomarkers.
- number of top N biomarkers for each of feature filters: ***number.gene = 100***
- how many times a biomarker has occurred in m feature subsets: ***level.freq = 7***
```r
gene.top <- get.top.gene(list.imp.var.cv = result$selected.feature, level.freq = 7, number.gene = 100)
print(gene.top)

$utest
 [1] C1orf69   C2orf15   CHMP4C    L2HGDH    MCCC2     MGC42105  NET1      ONECUT1   S1PR5     SEZ6      SIX4      SLC2A3   
[13] TMEM62    TNFSF11   XDH       ZDHHC9    ABCC3     ALG6      CD2AP     GJA9      GTPBP5    IGFL2     LOC81691  ZNF107   
[25] BAIAP2L2  CCDC54    DCPS      FAM155B   KRT86     PRNP      PROC      TOM1L1    TTC39C    ZC3HAV1L  ABHD5     APH1B    
[37] CILP2     DHFR      GJA1      HOXA5     KIAA0406  NETO1     ONECUT2   RARRES2   ADAM28    KIR3DL1   LYPD1     NRXN1    
[49] PTGFRN    RAE1      UNC5CL    AOX1      C1orf126  C1orf198  C7orf16   CP        FOLR3     RALGPS2   SEMA3B    C11orf67 
[61] FAM105A   FAM66C    FUNDC1    IL11RA    LOC148709 PTK6      RMST      RUSC2     SLC2A5    TLE4      DNM3      ERBB2    
[73] FRK       KLHL4     LPHN3     MMP13     NUDT5     PPARG     SEPHS2    SLC24A3   TCTEX1D1  ZNF786    BMP5      C19orf52 
[85] C8orf85   GPER      LOC387646 NELF      SCN8A     SLC25A25  SPP2      UST       ZNF643   

$mcfs
[1] EMP2 GPT2

$mrmr
[1] RS1    GPT2   B3GNT3

$mdfs.1D
 [1] C2orf15      MCCC2        MGC42105     NET1         ONECUT1      ONECUT2      TMEM62       TOM1L1       XDH         
[10] ABCC3        C1orf126     SEZ6         SLC2A3       TNFSF11      BAIAP2L2     CD2AP        CHMP4C       CILP2       
[19] DCPS         FAM155B      IGFL2        KIAA0406     LOC81691     ZC3HAV1L     ALG6         APH1B        CCDC54      
[28] HOXA5        KRT86        L2HGDH       NUDT5        PROC         SCN8A        SLC12A8      TTC39C       UNC5CL      
[37] UST          ZDHHC9       ADAM28       C1orf69      DHFR         DNM3         GJA1         KIR3DL1      LOC148709   
[46] LYPD1        PTK6         SIX4         TSPYL5       ALPL         CP           FOLR3        GTPBP5       RALGPS2     
[55] S1PR5        SEMA3B       TLE4         C11orf9      C1orf198     CCDC43       HIST1H2BE    LOC387646    PRNP        
[64] RAE1         ZNF107       ABHD5        ERBB2        FAM105A      GJA9         GRIP1        IL11RA       LGALS3BP    
[73] MAL          NELF         PLAC8        PTGFRN       QPRT         BMP5         C11orf67     C8orf85      CLCN4       
[82] KATNA1       KBTBD12      LOC100125556 NRXN1        RARRES2      SLC24A3      SPP2         UBFD1        ZNF643      
[91] ZNF786

```

#### Getting information about biomarkers from nine biological databases

- combination of a set of biomarkers: ***union***
```r
info.gene <- get.info.top.gene(gene.top, condition.methods = 'union')
print(info.gene)

        term source               term.ID                                                                    term.name
1     L2HGDH  GO:MF            GO:0047545                                    2-hydroxyglutarate dehydrogenase activity
2     L2HGDH     HP            HP:0040147                                                 L-2-hydroxyglutaric acidemia
3     L2HGDH   REAC     REAC:R-HSA-880009                     Interconversion of 2-oxoglutarate and 2-hydroxyglutarate
4     L2HGDH     WP             WP:WP4519                               Cerebral organic acidurias, including diseases
5      MCCC2  CORUM            CORUM:7190                                             3-methylcrotonyl-CoA carboxylase
6      MCCC2  GO:CC            GO:0002169                      3-methylcrotonyl-CoA carboxylase complex, mitochondrial
7      MCCC2  GO:CC            GO:1905202                                      methylcrotonoyl-CoA carboxylase complex
8      MCCC2  GO:MF            GO:0004485                                     methylcrotonoyl-CoA carboxylase activity
9      MCCC2     WP             WP:WP5031                                            Biotin metabolism, including IMDs
10   ONECUT1  CORUM             CORUM:746                                                      C/EBPalpha-HNF6 complex
11   ONECUT1   KEGG            KEGG:04950                                         Maturity onset diabetes of the young
12      SIX4     WP             WP:WP3595          mir-124 predicted interactions with cell cycle and differentiation 
13   TNFSF11  GO:BP            GO:0051466             positive regulation of corticotropin-releasing hormone secretion
14   TNFSF11     HP            HP:0004499                                  Chronic rhinitis due to narrow nasal airway
15   TNFSF11     TF             TF:M05659                                          Factor: ZNF879; motif: NGGTTTATAAKM
... etc.
```
## Example 3 - creating and adding a new feature filter to default list of basic FS methods
```r
data <- read.csv2('exampleData.csv')
class <- data$class
data$class <- NULL
```
####   Installing the required packages 
```r
install.packages('mt')
library(mt)
```
####  Creating a new feature filter eg. ReliefF ([Kononenko 1994](https://link.springer.com/chapter/10.1007/3-540-57868-4_57))

#### Rules for adding a new FS method to the benchmark procedure:
1. the function name of the new FS method must have the prefix 'fs.';
2. the created function should take as input arguments:
- ***x*** tabular input data, numeric type, where columns are features, and the rows are observations;
- ***y*** decision variable as a binary vector, y length equals the total number of observations;
- ***params*** list of hyperparameters for a new FS method, if such are required, and in their absence, the params argument of the above-mentioned
  function should be omitted: ***params = list(feature.number = 100)***
3. function should return a data frame that includes two columns: ***name*** (names of relevant biomarkers) and ***score*** (value of variable importance metric, e.g. p-value for U-test); 
**Code example:**
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
graph.result(result2$model, "acc")
graph.result(result2$stability, "stability")
```
![Fig.4](https://github.com/biocsuwb/Images/blob/main/ACC&ASM.png?raw=true)
Fig.4 The average values for accurancy (ACC) vs N top features (5, 10, 15, 20, ..., 50, 75, 100) for various features filters and the ASM similarity measure between m = 15 feature subsets vs N top features.
