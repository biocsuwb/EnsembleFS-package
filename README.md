# ensembleFS: Ensemble feature selection methods for analysis of molecular data

**Features**
* Automatic optimal predictor ensembling via cross-validation with one line of code.
*
*
*
*
*
*

### Install the development version from GitHub:

```r
# install.packages("devtools")
devtools::install_github("biocsuwb/ensembleFS")
```
## Examples 

ensembleFS makes it trivial to run many algorithms and use the best one or an ensemble.

```r

data <- read.csv2('examplesDataTest.csv')
class <- data$class
data$class <- NULL

list.methods()

result <- ensembleFS(x = data,
                     y = class,
                     methods = c("fs.utest", "fs.mcfs", "fs.mrmr", "fs.mdfs.1D", "fs.mdfs.2D"),
                     method.cv = "kfoldcv",
                     params.cv = list(k = 3, niter = 5),
                     level.cor = 1,
                     params = list(adjust = "holm", feature.number = 10, alpha = 0.05),
                     asm = c("fs.utest", "fs.mcfs", "fs.mrmr", "fs.mdfs.1D", "fs.mdfs.2D"),
                     model = c("fs.utest", "fs.mcfs", "fs.mrmr", "fs.mdfs.1D", "fs.mdfs.2D"))
                     

graph.result(result$stability, "stability")
graph.result(result$model, "auc")

gene.top <- get.top.gene(result$selected.feature, 15 , 20)
info.gene <- gene.info.top.gene(gene.top, condition.methods = 'union')


var.utest <- fs.utest(x = data, y = class, params = list(adjust = "holm", alpha = 0.05))

var.mcfs <- fs.mcfs(x = data, y = class)


list.index.cross <- cross.val(x = data,
                              y = class,
                              method = 'kfoldcv',
                              params.cv = list(niter = 10, k = 3)
                              
list.selected.var <- feature.selection.cv(x = data,
                                          y = class,
                                          method = 'fs.mdfs.2D',
                                          list.index.cross = list.index.cross,
                                          params = list(adjust = 'holm', alpha = 0.05)
 
asm <- stabilty.selection(list.selected.var, list.index.cross, 100)

model.result <- build.model.crossval(x = data,
                                     y = class,
                                     list.selected.var = list.selected.var,
                                     list.index.cross = list.index.cross,
                                     nvar = 100)



```
