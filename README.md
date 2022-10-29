# EnsembleFS: Ensemble feature selection methods for analysis of molecular data


### Install the development version from GitHub:

```r
# install.packages("devtools")
devtools::install_github("biocsuwb/EnsembleFS-package")
```
#### Note: system require Java (>= 7)
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
