# EnsembleFS: Ensemble feature selection methods for the analysis of molecular omics data
## Description
EnsembleFS is an R package for single feature selection (FS) and ensemble feature selection of molecular data or clinical data (numeric data formats).
This tool is based on several feature filters, such as the test Manna-Whitney'a (U-test), the Monte Carlo Feature Selection (MCFS) ([Dramiński & Koronacki 2018](https://www.jstatsoft.org/article/view/v085i12)), the MultiDimensional Feature Selection (two variants: MDFS-1D and MDFS-2D) ([Mnich & Rudnicki 2020](https://www.sciencedirect.com/science/article/abs/pii/S0020025520302048)), and the Minimum Redundancy Maximum Relevance (MRMR) 
([Ding 2005](https://pubmed.ncbi.nlm.nih.gov/15852500/)) for discovering the most important biomarkers and using machine learning algorithms (ML) to evaluate the quality of feature sets. Predictive models are built using the Random Forest algorithm ([Breiman 2001](https://link.springer.com/article/10.1023/A:1010933404324)). It can be applied to two-class problems.

Moreover, EnsembleFS supports users in the analysis and interpretation of molecular data. The information about each of the top biomarkers is extracted from diverse biological databases, namely the Gene Ontology ([GO](https://pubmed.ncbi.nlm.nih.gov/33290552/)), the Kyoto Encyclopedia of Genes and Genomes ([KEGG](https://pubmed.ncbi.nlm.nih.gov/18477636/)), the Reactome ([React](https://pubmed.ncbi.nlm.nih.gov/32907876/)), the WikiPathways ([WP](https://pubmed.ncbi.nlm.nih.gov/33211851/)), the Transfac ([TF](https://pubmed.ncbi.nlm.nih.gov/8594589/)), the miRTarBase ([MIRNA](https://academic.oup.com/nar/article/48/D1/D148/5606625)), the Human Protein Atlas ([HPA](https://pubmed.ncbi.nlm.nih.gov/25613900/)), the [CORUM](https://pubmed.ncbi.nlm.nih.gov/30357367/), and the Human Phenotype Ontology ([HPO](https://pubmed.ncbi.nlm.nih.gov/33264411/)).
The proposed tool accepts molecular data that includes different gene identifiers, such as Ensembl, NCBI Entrez gene ID, Refseq, Illumina, and Uniprot.

EnsembleFS allows the user to:
- filter the most informative features by using up to five FS methods (U-test, MCFS, MRMR, MDFS-1D, and MDFS-2D) from molecular data generated from high-throughput molecular biology experiments
and also clinical data (numeric formats);
- create and add their own feature filters to the default list of basic feature filters (U-test, MCFS, MRMR, MDFS-1D, and MDFS-2D) for more complex FS tasks;
- build Random Forest classifiers using selected top N features (Fig.1); 
- evaluate the stability of feature subsets (the Lustgarten adjusted stability measure ([ASM](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2815476/)) ) and the performance of predictive models (the area under receiver operator curve  (AUC), accuracy (ACC), and the Matthews correlation coefficient(MCC);
- compare the predictive performance of models and the stability of selected feature sets for selected FS algorithms; 
- establish the selected parameters for predictive models, such as the number of top N informative features;
- remove redundant features by building the Spearman correlation matrix that identifies highly correlated features;
- product plots to visualize the model results;
- find information about selected top molecular markers (gene ontology, pathways, tissue specificity, miRNA targets, regulatory motif, protein complexes, disease phenotypes) in nine biological databases (GO, KEGG, React, WP, TF, MIRNA, HPA, CORUM, and HPO).

![Fig.1](https://github.com/biocsuwb/Images/blob/main/Scheme1.png?raw=true)
Fig.1 The scheme of ensemble feature selection and supervised classification.
![Fig.2](https://github.com/biocsuwb/Images/blob/main/Scheme2.png?raw=true)
Fig.2 The construction scheme of the combined set of N-top relevant features.
![Fig.3](https://github.com/biocsuwb/Images/blob/main/Scheme3.png?raw=true)
Fig.3 The scheme for collection and integration of biological information about biomarkers.


## Functional modules of EnsembleFS R package 
Table of main functions
[Functional modules of EnsembleFS R package.pdf](https://github.com/biocsuwb/Images/blob/9b5d8bb15d0a0b2e0cc6ec730cdb14e7ecd99314/Functional%20modules%20of%20EnsembleFS%20R%20package.pdf)

## Install the development version from GitHub:
```r
install.packages("devtools")
devtools::install_github("biocsuwb/EnsembleFS-package")
```
## Notes: 
- ***to install the EnsembleFS package in your R environment, make sure that you have Java installed (rJava R package);***
- ***to accelerate processing by using a CUDA GPU, the EnsembleFS package must be compiled with CUDA (set the MDFS-2D parameter: use.cuda = TRUE).*** 

## Example data sets
The RNA-sequencing data of tumor-adjacent normal tissues of lung adenocarcinoma cancer patients from The Cancer Genome Atlas database ([TCGA](https://www.cancer.gov/tcga)) was used. The preprocessing of data involved standard steps for RNA-Seq data. The log2 transformation was performed. Features with zero and near-zero (1%) variance across patients were removed. After preprocessing, the primary dataset contains 574 samples (59 normal and 515 tumors) described with 20172 differentially expressed genes (DEGs). This dataset includes highly correlated features, and the number of cancer samples is roughly ten times more than normal samples. For testing purposes, the number of molecular markers was limited to 2000 DEGs ranked by the highest difference in the gene expression level between tumor and normal tissues ([exampleData_TCGA_LUAD_2000.csv](https://github.com/biocsuwb/EnsembleFS-package/tree/main/data)). 

The clinical data of 1394 breast cancer patients from Molecular Taxonomy of Breast Cancer International Consortium (METABRIC) project ([Pereira et. al 2016](https://www.jstatsoft.org/article/view/v085i12](https://www.nature.com/articles/ncomms11479))). This data is available at ([exampleData_METABRIC_clinical.csv](https://github.com/biocsuwb/EnsembleFS-package/tree/main/data)). First column ("class") includes the clinical endpoints of patients. 

## Example 1 - individual feature selection

#### Loading data
```r
download.file("https://raw.githubusercontent.com/biocsuwb/EnsembleFS-package/main/data/exampleData_TCGA_LUAD_2000.csv", 
              destfile = "exampleData_TCGA_LUAD_2000.csv", method = "curl")

data <- read.csv2('exampleData_TCGA_LUAD_2000.csv')
decisions <- data$class
data$class <- NULL
```

#### Model (optional) configuration parameters
- U-test and MDFS-1D, and MDFS-2D parameter, multitest correction:  ***adjust = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")***;
- U-test and MDFS-1D, and MDFS-2D parameter, significance level: ***alpha = 0.05***;
- MDFS-2D parameter, CPU/GPU architecture: ***use.cuda = FALSE*** 
- MRMR parameter, number of significant features: ***feature.number = 10***;
- MCFS parameter, cut-off method: ***cutoff.method = c("permutations", "criticalAngle", "kmeans", "mean", "contrast")***;
- validation methods: ***method.cv = c('kfoldcv','rsampling')***;
- number of repetitions: ***niter = 10***;
- train-test-split the data: ***k = 3*** for stratified k-fold cross-validation and ***test.size = 0.3*** for random sampling.


#### Feature selection U-test
```r
library (ensembleFS)
var.utest <- fs.utest(x = data, y = decisions, params = list(adjust = "holm", alpha = 0.05))
```
#### Feature selection MDFS-1D
```r
var.mdfs1D <- fs.mdfs.1D(x = data, y = decisions, params = list(adjust = "holm", alpha = 0.05))
```
#### Feature selection MDFS-2D
```r
var.mdfs2D <- fs.mdfs.2D(x = data, y = decisions, params = list(adjust = "holm", alpha = 0.05, use.cuda = FALSE))
```
#### Feature selection MCFS
```r
var.mcfs <- fs.mcfs(x = data, y = decisions,  params = list(cutoff.method = "kmeans"))
```
#### Feature selection MRMR
```r
var.mrmr <- fs.mrmr(x = data, y = decisions,  params = list(feature.number = 10))
```

#### Creating the cross-validation index array  (3-fold cross-validation repeated 10 times)
```r
list.index.cross <- cross.val(x = data,
                              y = decisions,
                              method = 'kfoldcv',
                              params.cv = list(k = 3, niter = 10)
                              )
```                              
#### Feature selection in a cross-validation scenario for individual FS method eg. MDFS-1D
```r
list.selected.var <- feature.selection.cv(x = data,
                                          y = decisions,
                                          method = 'fs.mdfs.1D',
                                          list.index.cross = list.index.cross,
                                          params = list(adjust = 'holm', alpha = 0.05)
                                          )
 
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

#### Building and testing ML models on top 100 features with the MDFS-1D method; 30 models in total.
```r
model.result <- build.model.crossval(x = data,
                                     y = decisions,
                                     list.selected.var = list.selected.var,
                                     list.index.cross = list.index.cross,
                                     nvar = 100)
# show a predictive model result (list 1 of 30)
print(model.result[[1]]) 

 Accuracy       AUC       MCC 
0.9807692 0.9736842 0.9589080 
```
#### Computing ASM value with 30 sets of top 100 features
```r
asm <- stabilty.selection(list.selected.var, list.index.cross,  nvar = 100)
# show ASM value 
print(asm)
0.6915541
```
## Example 2 - ensemble feature selection

#### Loading data
```r
data <- read.csv2('exampleData_TCGA_LUAD_2000.csv')
decisions <- data$class
data$class <- NULL
```

#### Showing the list of available feature selection methods
```r
print(list.methods())
"fs.mcfs"  "fs.mdfs.1D"  "fs.mdfs.2D"  "fs.mrmr"  "fs.utest"  
```

#### Model (optional) configuration parameters
EnsembleFS allows the user to set some parameter values, such as:
- feature selection methods: ***methods = c("fs.utest", "fs.mcfs", "fs.mrmr", "fs.mdfs.1D", "fs.mdfs.2D")***;
- U-test and MDFS parameter, multitest correction: ***adjust = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")***;
- U-test and MDFS parameter, significance level: ***alpha = 0.05***;
- MRMR parameter, number of significant features: ***feature.number = 100***;
- MCFS parameter, cut-off method: ***cutoff.method = c("permutations", "criticalAngle", "kmeans", "mean", "contrast")***;
- correlation coefficient: ***level.cor = 1***;
- validation methods: ***method.cv = c('kfoldcv','rsampling')***;
- number of repetitions: ***niter = 10***;
- train-test-split the data: ***k = 3*** for stratified k-fold cross-validation and ***test.size = 0.3*** for random sampling.


#### Building and testing ML models on top N features (N = 5, 10, 15, 20, 30, 40, 50, 75, 100) with each of the selected feature filters
- selected feature filters: ***methods = c("fs.utest", "fs.mcfs", "fs.mrmr", "fs.mdfs.1D")***;
- U-test and MDFS parameter: ***adjust = "holm"***;
- MRMR parameter: ***feature.number = 150;***;
- MCFS parameter: ***cutoff.method = "kmeans"***;
- model validation technique: ***rsampling, test set 30%, repeated 10 times;***
- the cut off value of the Spearman correlation coefficient: ***0.75;***
```r
result <- ensembleFS(x = data,
                     y = decisions,
                     methods = c("fs.utest", "fs.mcfs", "fs.mrmr", "fs.mdfs.1D"),
                     method.cv = "rsampling",
                     params.cv = list(test.size = 0.3, niter = 10),
                     level.cor = 0.75,
                     params = list(adjust = "holm", cutoff.method = "kmeans",feature.number = 150, alpha = 0.05),
                     asm = c("fs.utest", "fs.mcfs", "fs.mrmr", "fs.mdfs.1D"),
                     model = c("fs.utest", "fs.mcfs", "fs.mrmr", "fs.mdfs.1D")
                     )
 ```
                     
#### Visualizing the model results
```r
graph.result(result$model, "acc")
graph.result(result$stability, "stability")
# graph.result(result$model, "auc")
# graph.result(result$model, "mcc")
```
![Fig.4](https://github.com/biocsuwb/Images/blob/main/ACC&ASM.png?raw=true)
Fig.4 The average values for accurancy (ACC) vs N top features (5, 10, 15, 20, ..., 50, 75, 100) for various features filters and the ASM similarity measure between m = 10 feature subsets vs N top features.
#### Showing m-list of top biomarkers for each of feature filters (m = 10 in this case)
```r
print(result$selected.feature)

$fs.utest[[1]]
           name       Pvalue adjustPvalue
1         OR6K3 8.581769e-32 1.716354e-28
2        CELA2B 1.483157e-29 2.964831e-26
3         CHRM2 1.862059e-28 3.720393e-25
4         STX11 7.275286e-28 1.452875e-24
5          EMP2 9.019200e-28 1.800232e-24
6      C13orf36 9.737433e-28 1.942618e-24
7         TNNC1 9.968668e-28 1.987752e-24
8          AGER 1.040613e-27 2.073941e-24
9          FHL1 1.166631e-27 2.323929e-24
10         SDPR 1.508215e-27 3.002855e-24
12         WWC2 1.666373e-27 3.314415e-24
14       CHRNA2 1.768252e-27 3.513517e-24
15         SGCG 1.771934e-27 3.519061e-24
17  PALM2.AKAP2 1.789309e-27 3.551778e-24
... etc.
```

#### Showing the combined list of top biomarkers.
- number of top N biomarkers for each of feature filters: ***number.gene = 100***
- how many times a biomarker has occurred in m feature subsets: ***level.freq = 5***
```r
gene.top <- get.top.gene(list.imp.var.cv = result$selected.feature, level.freq = 5, number.gene = 100)
print(gene.top)

$utest
  [1] ADCY8       ADRB2       ANGPTL7     C10orf67    C13orf36    C9orf140    CAT         CAV1        CCDC85A    
 [10] CELA2B      CGNL1       CHRM2       CHRNA2      CYP1A2      EFNA3       EMP2        EPAS1       FAM189A2   
 [19] FAM46B      GCOM1       GOLM1       GPT2        GYPE        KAL1        MGAT3       NCKAP5      NECAB1     
 [28] ODAM        OTC         OTUD1       OVCH2       PALM2.AKAP2 PDLIM2      PRKCE       PTPN21      PTPRQ      
 [37] PYCR1       RCC1        RTKN2       SCUBE1      SGEF        SH3GL3      SLC25A10    SLC39A8     STX11      
 [46] TACC1       TNNC1       UPK3B       WNT3A       WWC2        ARHGAP31    C14orf132   C20orf202   CD5L       
 [55] DPYSL2      FAM83A      KCNT2       LIMS2       OR6K3       SGCG        TMEM184A    AGER        C10orf116  
 [64] CCDC141     DOCK4       EFNA4       FABP4       FGF10       GPM6A       LOC401093   PAFAH1B3    PAICS      
 [73] RNF144B     SDPR        STX1A       STXBP6      ACADL       AMOTL1      ANGPT4      C16orf59    C19orf59   
 [82] FGF11       FHL1        MEX3A       OCIAD2      OVCH1       PECAM1      PHACTR1     SASH1       SPAG4      
 [91] ADM2        ADRA1A      ETV4        GRK5        KIAA1324L   PPP1R14B    PRKG2       SPOCK2      AGBL1      
[100] DENND3      LOC158376   MYOC        THBD        TNS1       

$mcfs
 [1] ADRB2       C10orf67    C14orf132   CAT         CBLC        EMP2        ETV4        FAM83A      GYPE       
[10] PYCR1       SGCG        SH3GL3      STX11       STX1A       WWC2        CLIC5       EFNA3       FABP4      
[19] GOLM1       GPM6A       PALM2.AKAP2 PDLIM2      SLC25A10    TNNC1       CAV1        EFNA4       FAM46B     
[28] TMEM184A    ACADL       SDPR        SLC39A8     C16orf59    KCNT2       LIMS2       ODAM        PTPRQ      
[37] LOC401093   SFTPC       TGFBR2     

$mrmr
  [1] ADCY8       ADM2        ALAS2       ANGPT4      ANGPTL7     ATP10B      B3GNT3      C10orf67    C13orf36   
 [10] C2orf71     CD101       CD5L        CELA2B      CHRM2       CHRNA2      COL10A1     CYP1A2      EFNA4      
 [19] EMP2        EPAS1       ETV4        FAM83A      FERMT1      FGF10       FUT2        GOLM1       GPT2       
 [28] HTR3C       KIF26B      LGR4        LOC158376   MYO7A       MYOC        OR6K3       OTC         OTUD1      
 [37] OTX1        OVCH1       OVCH2       PALM2.AKAP2 PLEKHN1     PPP1R15A    PTPN21      PYCR1       RASAL1     
 [46] RCC1        RS1         RTKN2       SGCG        SGEF        SH3GL3      SPTBN2      STX11       TMEM184A   
 [55] TMPRSS4     WWC2        C20orf202   CACNA1S     CAV1        CYP3A7      FABP4       FGF11       IGSF9      
 [64] KCNA4       NCKAP5      SLC4A1      SPOCK2      WNT3A       C16orf59    CMTM2       E2F8        GP9        
 [73] HS6ST2      OCIAD2      ODAM        SPAG4       SPP1        COLEC10     EPN3        MEX3A       OR2W3      
 [82] PECAM1      PPAP2C      PRKG2       PTPN5       PTPRQ       SEMA3G      WFS1        CLEC4M      CTHRC1     
 [91] EFNA3       FAM46B      GYPE        KCNT2       PRKCE       RASIP1      TGM1        C13orf15    MND1       
[100] SLC39A8     SLC6A4      TUBB3      

$mdfs.1D
 [1] RTKN2     SGEF      STX1A     C8orf84   DPYSL2    GPT2      LOC401093 MEX3A     PIP5K1B   PPP1R14B  PRKG2    
[12] ADM2      ANKRD29   C6orf155  CAT       CAV1      CCDC141   CELSR3    CYP1A2    ESAM      FAM46B    FAM83A   
[23] KAL1      LIN7A     PLA2G4F   PTPRQ     RADIL     RCC1      RNF144B   SLC39A8   TNNC1     TUBB3     ADRB1    
[34] AMOTL1    ANXA3     B3GNT3    C10orf116 C10orf67  CCDC85A   DNAH14    EMP2      ETV4      FABP4     GCOM1    
[45] HSD17B6   IGSF10    KCNT2     KL        NCKAP5    NECAB1    NKAPL     NMUR1     OTUD1     PAFAH1B3  PLEKHN1  
[56] PRKCE     PYCR1     SEMA5A    STXBP6    TACC1     TFR2      WFS1     
```

#### Showing the predictive performance of the random forest model
```r
print(result$model)
     N  mean.acc  mean.auc  mean.mcc      sd.acc     sd.auc     sd.mcc  method
1    5 0.9847001 0.9476866 0.9118262 0.011247365 0.04758166 0.05924979   utest
2   10 0.9899995 0.9618989 0.9408039 0.004710309 0.03878885 0.03258368   utest
3   15 0.9926742 0.9739264 0.9581730 0.003655541 0.02114019 0.02050122   utest
4   20 0.9926742 0.9739104 0.9589725 0.003655541 0.02112748 0.01785854   utest
5   30 0.9926144 0.9738867 0.9578773 0.003885302 0.02171403 0.02105812   utest
6   40 0.9926144 0.9738867 0.9578773 0.003885302 0.02171403 0.02105812   utest
7   50 0.9926742 0.9739104 0.9589725 0.003655541 0.02112748 0.01785854   utest
8   75 0.9926742 0.9739104 0.9589725 0.003655541 0.02112748 0.01785854   utest
9  100 0.9926742 0.9739104 0.9589725 0.003655541 0.02112748 0.01785854   utest
10   5 0.9926910 0.9692675 0.9558952 0.006565496 0.03992986 0.04390362    mcfs
11  10 0.9919811 0.9665294 0.9538377 0.005186949 0.03098983 0.03094842    mcfs
12  15 0.9926144 0.9709342 0.9588457 0.004946283 0.02515812 0.02533537    mcfs
13  20 0.9926144 0.9709342 0.9588457 0.004946283 0.02515812 0.02533537    mcfs
14  30 0.9939729 0.9746406 0.9653022 0.003727046 0.02195584 0.02148852    mcfs
15  40 0.9933236 0.9742675 0.9627090 0.004333648 0.02157585 0.02211227    mcfs
16  50 0.9939729 0.9746406 0.9653022 0.003727046 0.02195584 0.02148852    mcfs
17  75 0.9939729 0.9746406 0.9653022 0.003727046 0.02195584 0.02148852    mcfs
18 100 0.9933236 0.9742675 0.9627090 0.004333648 0.02157585 0.02211227    mcfs
19   5 0.9865872 0.9428229 0.9225840 0.009476662 0.04228822 0.05398912    mrmr
20  10 0.9905535 0.9603092 0.9464660 0.006657340 0.03148771 0.03645462    mrmr
21  15 0.9906906 0.9662615 0.9470864 0.005460468 0.02668610 0.03061786    mrmr
22  20 0.9920076 0.9707854 0.9553583 0.004102557 0.02408033 0.02009801    mrmr
23  30 0.9939729 0.9746406 0.9653022 0.003727046 0.02195584 0.02148852    mrmr
24  40 0.9939729 0.9746406 0.9653022 0.003727046 0.02195584 0.02148852    mrmr
25  50 0.9939729 0.9746406 0.9653022 0.003727046 0.02195584 0.02148852    mrmr
26  75 0.9933236 0.9742675 0.9627090 0.004333648 0.02157585 0.02211227    mrmr
27 100 0.9933236 0.9742675 0.9627090 0.004333648 0.02157585 0.02211227    mrmr
28   5 0.9887040 0.9599451 0.9308505 0.005317968 0.03666072 0.03844549 mdfs.1D
29  10 0.9906830 0.9700560 0.9466745 0.005505444 0.02799906 0.03005919 mdfs.1D
30  15 0.9913536 0.9723194 0.9501081 0.004386496 0.02890808 0.02507787 mdfs.1D
31  20 0.9933280 0.9761623 0.9616924 0.005341063 0.03030464 0.03140701 mdfs.1D
32  30 0.9926786 0.9736783 0.9575776 0.004807890 0.02944285 0.02862076 mdfs.1D
33  40 0.9933409 0.9778450 0.9623404 0.004334037 0.02264081 0.02245097 mdfs.1D
34  50 0.9926742 0.9717835 0.9585942 0.003655541 0.02008900 0.01819328 mdfs.1D
35  75 0.9926144 0.9738867 0.9578773 0.003885302 0.02171403 0.02105812 mdfs.1D
36 100 0.9932637 0.9742438 0.9616138 0.004538666 0.02215096 0.02495123 mdfs.1D
```

#### Getting information about biomarkers from the nine biological databases (GO, KEGG, React, WP, TF, MIRNA, HPA, CORUM, and HPO)

- combination of a set of biomarkers: ***union***
```r
info.gene <- get.info.top.gene(gene.top, condition.methods = 'union')
info.gene
         term source               term.ID                                                                                                                                      term.name
1       ADCY8  GO:MF            GO:0008294          calcium- and calmodulin-responsive adenylate cyclase activity
2       ADRB2  CORUM            CORUM:3830          ADRB2 homodimer complex
3       ADRB2  CORUM             CORUM:668          BKCA-beta2AR-AKAP79 signaling complex
4       ADRB2  CORUM             CORUM:672          BKCA-beta2AR complex
5       ADRB2  CORUM             CORUM:687          CFTR-NHERF-beta(2)AR signaling complex
6       ADRB2  GO:MF            GO:0004941          beta2-adrenergic receptor activity
7     ANGPTL7  GO:BP            GO:0036331          avascular cornea development in camera-type eye
8     ANGPTL7  GO:BP            GO:1901346          negative regulation of vasculature development involved in avascular cornea development in camera-type eye
9         CAT  GO:BP            GO:0061691          detoxification of hydrogen peroxide
10        CAT  GO:BP            GO:0061692          cellular detoxification of hydrogen peroxide
11        CAT  GO:CC            GO:0062151          catalase complex
12        CAT     HP            HP:0012517          Reduced catalase level
13        CAT     HP            HP:0040113          Old-aged sensorineural hearing impairment
14       CAV1  CORUM            CORUM:2462          Caveolin-1 homodimer complex
15       CAV1  CORUM            CORUM:5714          NOS3-CAV1 complex
16       CAV1  CORUM             CORUM:550          NOS3-CAV1-NOSTRIN complex
17       CAV1  CORUM            CORUM:5862          CAV1-VDAC1-ESR1 complex
18       CAV1  GO:BP            GO:1900085          negative regulation of peptidyl-tyrosine autophosphorylation
19       CAV1  GO:BP            GO:1903609          negative regulation of inward rectifier potassium channel activity
20       CAV1  GO:CC            GO:0002095          caveolar macromolecular signaling complex
21       CAV1  GO:MF            GO:0070320          inward rectifier potassium channel inhibitor activity
22       CAV1     HP            HP:0005320          Lack of facial subcutaneous fat
23     CHRNA2  CORUM            CORUM:6440          CHRNA2-CHRNB4 complex
24     CYP1A2  GO:BP            GO:0009403          toxin biosynthetic process
25     CYP1A2   KEGG            KEGG:00232          Caffeine metabolism
26     CYP1A2   REAC     REAC:R-HSA-211957          Aromatic amines can be N-hydroxylated or N-dealkylated by CYP1A2
27     CYP1A2   REAC    REAC:R-HSA-9018681          Biosynthesis of protectins
28     CYP1A2     TF             TF:M06228          Factor: ZNF543; motif: KGGWATRTGGGA
29     CYP1A2     WP             WP:WP2646          Lidocaine metabolism
30     CYP1A2     WP             WP:WP3633          Caffeine and theobromine metabolism
31     CYP1A2     WP              WP:WP694          Arylamine metabolism
32     CYP1A2     WP             WP:WP2542          Sulindac metabolic pathway
33     CYP1A2     WP              WP:WP699          Aflatoxin B1 metabolism
34      EPAS1  CORUM            CORUM:7274          ARNTL-EPAS1 complex
... etc.
```


## Example 3 - creating and adding a new feature filter to the default list of basic FS methods
```r
data <- read.csv2('exampleData.csv')
decisions <- data$class
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
- ***x*** tabular input data, a numeric type, where columns are features, and the rows are observations;
- ***y*** decision variable as a binary vector, y length equals the total number of observations;
- ***params*** list of hyperparameters for a new FS method if such is required, and in their absence, the params argument of the above-mentioned
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
- selected feature filters: ***methods = c("fs.utest", "fs.mcfs", "fs.mrmr", "fs.mdfs.1D", "fs.relieff")***;
- U-test and MDFS parameter: ***adjust = "fdr"***;
- MRMR parameter: ***feature.number = 150;***;
- MCFS parameter: ***cutoff.method = "kmeans"***;
- model validation technique: ***3-fold cross-validation, repeated 5 times;***
- the cut-off value of the Spearman correlation coefficient: ***0.75;***
```r
result2 <- ensembleFS(x = data,
                     y = decisions,
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
graph.result(result2$model, "mcc")
graph.result(result2$stability, "stability")
```
![Fig.4](https://github.com/biocsuwb/Images/blob/main/MCC&ASM.png?raw=true)
Fig.4 The average values of the Matthews correlation coefficient (MCC) vs N top features (5, 10, 15, 20, ..., 50, 75, 100) for various features filters and the ASM similarity measure between m = 15 feature subsets vs N top features.

## Example 4 - identification of relevant features from clinical data

#### Loading data
```r
data <- read.csv2(https://raw.githubusercontent.com/biocsuwb/EnsembleFS-package/main/data/exampleData_METABRIC_clinical.csv')
# data <- read.csv2('exampleData_METABRIC_clinical.csv')
decisions <- data$class
data$class <- NULL
```
#### Building and testing ML model on top N features with each of the selected feature filters
- selected feature filters: ***methods = c("fs.utest", "fs.mcfs", "fs.mrmr", "fs.mdfs.1D")***;
- U-test and MDFS parameter: ***adjust = "holm"***;
- MRMR parameter: ***feature.number = 10;***;
- MCFS parameter: ***cutoff.method = "kmeans"***;
- model validation technique: ***rsampling, test set 30%, repeated 10 times;***
- the cut off value of the Spearman correlation coefficient: ***0.75;***
```r
result <- ensembleFS(x = data,
                     y = decisions,
                     methods = c("fs.utest", "fs.mcfs", "fs.mrmr", "fs.mdfs.1D"),
                     method.cv = "kfoldcv",
                     params.cv = list(k = 3, niter = 10),
                     level.cor = 0.75,
                     params = list(adjust = "holm", cutoff.method = "kmeans", feature.number = 10, alpha = 0.05),
                     asm = c("fs.utest", "fs.mcfs", "fs.mrmr", "fs.mdfs.1D"),
                     model = c("fs.utest", "fs.mcfs", "fs.mrmr", "fs.mdfs.1D")
                     )

#### Showing the list of top features.
- maximal number of top N biomarkers for each of individual feature filters: ***number.gene = 10***
- how many times a biomarker has minimal occurred in m feature subsets: ***level.freq = 15***
```r
feature.top <- get.top.gene(list.imp.var.cv = result$selected.feature, level.freq = 5, number.gene = 10)
print(feature.top)

result <- ensembleFS(x = data,
                     y = decisions,
                     methods = c("fs.utest", "fs.mcfs", "fs.mdfs.1D"),
                     method.cv = "kfoldcv",
                     params.cv = list(k = 3, niter = 10),
                     level.cor = 0.75,
                     params = list(adjust = "holm", cutoff.method = "kmeans", alpha = 0.05),
                     asm = c("fs.utest", "fs.mcfs", "fs.mdfs.1D"),
                     model = c("fs.utest", "fs.mcfs", "fs.mdfs.1D")
)

feature.top <- get.top.gene(list.imp.var.cv = result$selected.feature, level.freq = 15, number.gene = 10)
print(feature.top)

$utest
[1] AGE_AT_DIAGNOSIS   BREAST_SURGERY    GRADE         NPI             TUMOR_SIZE      
[6] TUMOR_STAGE        HER2_STATUS       COHORT        ONCOTREE_CODE   


$mcfs
[1] AGE_AT_DIAGNOSIS   NPI    TUMOR_SIZE     BREAST_SURGERY   COHORT          

$mdfs.1D
[1] AGE_AT_DIAGNOSIS   NPI              TUMOR_SIZE       TUMOR_STAGE       COHORT          
[6] GRADE              BREAST_SURGERY   THREEGENE        CLAUDIN_SUBTYPE 

```

