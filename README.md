
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BioM2: Biologically-informed multi-stage machine learning for phenotype prediction using omics data <img src="https://github.com/jkkomm/img/blob/main/BioM2.png" align="right" height="180" />

<!-- Add buttons here -->
  
  ![GitHub release (latest by date including pre-releases)](https://img.shields.io/github/v/release/navendu-pottekkat/awesome-readme?include_prereleases)



## Motivation
Identifying reproducible and interpretable biological patterns from high-dimensional omics data is a critical factor in understanding the risk mechanism of complex disease. As such, explainable machine learning can offer biological insight in addition to personalized risk scoring.

## Core Workflow
![BioM2](https://github.com/jkkomm/img/blob/main/Architect.jpg)
## Deliverables
We have implemented a biologically informed multi-stage machine learning framework termed __BioM2__ specifically for phenotype prediction using omics-scale data based on prior biological information including gene ontological (GO) and/or KEGG pathways.   

**Features of BioM2 in a nutshell**:   
  
1. Phenotype prediction using whole-genome DNA methylation data and genome-wide gene expression data.
2. Intrinsic ranking of outcome-associated functional patterns.    
3. Functional patterns can be used for biological stratification on an individual subject basis.   
4. Modularisation of functional patterns and subsequent network analysis of these modules.   
5. Various choice of conventional machine learning models that can be integrated within the BioM2 framework.   

## :writing_hand: Authors

Shunjie Zhang  ----  (E-mail: 202220148490@mail.scut.edu.cn)

# Tutorial

## Installation
BioM2 has been uploaded to CRAN .
```
install.packages('BioM2')
```
The latest release can be installed using the code provided below.
```

install.packages("devtools")
library(devtools)
install_github("BioTransAI/BioM2")
```
BioM2 is built on the mlr3 package. To use additional learners, please install the mlr3extralearners package.
```
remotes::install_github("mlr-org/mlr3extralearners@*release")
```
## Sample Pairing (MatchIt)

MatchIt package can be employed to identify the most optimal pairs of samples by matching the experimental and control samples based on covariates such as the age and gender of the samples,which allows for increased robustness to the choice of model.
```
library(MatchIt)

#Meta include 'Label' and other covariates of interest
head(meta[,1:5])
                   Label Age Female       Cd8t       Cd4t
201172200060_R08C01     1  26      0 0.06150491 0.19702276
201172200055_R06C01     1  43      1 0.08574886 0.02965455
201810150072_R03C01     1  55      0 0.16006435 0.16938068
201800470034_R07C01     0  27      1 0.06554267 0.14906774
201172200012_R08C01     0  62      1 0.07471174 0.09312758
201800470020_R08C01     0  40      0 0.08413080 0.11795726

#Sample imbalance
 table(meta$Label)

  0   1 
 68 361


#Set minority class samples to 'TRUE'
meta$Label=as.logical(meta$Label == 0)


#Paired samples based on all covariates(setting method and pairing ratio)
meta_match <- matchit(Label~., data =meta, method="nearest", ratio=1)
summary(meta_match)

Summary of Balance for All Data:
             Means Treated Means Control Std. Mean Diff. Var. Ratio eCDF Mean
distance            0.2994        0.1320          0.8504     2.8508    0.2727
Age                47.0588       39.5623          0.5322     1.3622    0.1535
Female              0.6618        0.6122          0.1048          .    0.0496
Cd8t                0.0819        0.0740          0.1888     0.9874    0.0656
...
population_10        0.9117       -0.1822          0.2190     5.7112    0.0476
smoking             1.0109        2.5165         -0.3554     0.7408    0.0910



Summary of Balance for Matched Data:
             Means Treated Means Control Std. Mean Diff. Var. Ratio eCDF Mean
distance            0.2994        0.2642          0.1789     1.6232    0.0141
Age                47.0588       44.2206          0.2015     1.5343    0.0817
Female              0.6618        0.7059         -0.0933          .    0.0441
Cd8t                0.0819        0.0929         -0.2656     0.7571    0.0579
...
population_10        0.9117        0.3284          0.1168     6.4489    0.0413
smoking             1.0109        0.8027          0.0492     1.0853    0.0275



Sample Sizes:
          Control Treated
All           361      68
Matched        68      68
Unmatched     293       0
Discarded       0       0


#Extraction of paired samples
match <- match.data(meta_match)

```

## Data requirements

BioM2 requires:
- Genome-wide DNA methylation data / Bulk RNA-seq data
- Feature annotation data
- Pathway annotation data  (Gene ontological and KEGG pathways are used in this tutorial)

```
<<< Data requirements >>>
  
  -----------------------------------------------------------------------------
  ## [Genome-wide DNA methylation data]  (data.frame)
  # First column name must be 'label', and the rest are the features (e.g., CpGs).
  
  label cg21870274 cg09499020 cg16535257 cg00168193
0     0.0057     0.0002    -0.0313     0.0002
0    -0.0317    -0.0444    -0.0578    -0.0160
1    -0.0341    -0.0541    -0.0056    -0.0230
1     0.0811    -0.0029     0.0049     0.0274
1    -0.0187     0.0475     0.1168     0.0169
0    -0.0158     0.0032    -0.0173     0.0133

--------------------------------------------------------------------------------
  ## [Feature annotation data]  (data.frame)
  # The data frame must contain the two column names 'ID' and 'entrezID' .
  
  ID entrezID symbol
cg00000029     5934   RBL2
cg00000109    64778 FNDC3B
cg00000155   221927  BRAT1
cg00000221   162282 ANKFN1
cg00000236     7419  VDAC3
cg00000289       87  ACTN1

-----------------------------------------------------------------------------
  ## [Pathway annotation data] (list)
  # The name of each subset of the list is the ID of the pathway, and each subset contains a vector of gene entrezIDs.
  
  List of 15719
$ GO:0000002: chr [1:31] "142" "291" "1763" "1890" ...
$ GO:0000003: chr [1:1513] "2" "18" "49" "51" ...
$ GO:0000012: chr [1:12] "142" "1161" "2074" "3981" ...
$ GO:0000017: chr [1:2] "6523" "6524"
$ GO:0000018: chr [1:131] "60" "86" "142" "604" ...
$ GO:0000019: chr [1:7] "2068" "4292" "4361" "7014" ...
$ GO:0000022: chr [1:9] "4926" "6795" "9055" "9212" ...
$ GO:0000023: chr [1:3] "2548" "2595" "8972"

-----------------------------------------------------------------------------
  
  ```

## Phenotype prediction
To predict the phenotype, follow these steps.
```
library(mlr3verse)
library(caret)
library(parallel)
library(BioM2)

result=BioM2 (  TrainData = data , TestData = NULL ,                               ## If you only have one dataset
                pathlistDB = pathlistDB ,                                          ## ==>> [Pathway annotation data]
                FeatureAnno = FeatureAnno ,                                        ## ==>> [Feature annotation data]
                classifier = 'liblinear' , nfolds = 5 ,                            ## Choose your learner( use "lrns()" ) , currently only cross-validation is supported
                Inner_CV = FALSE , inner_folds=10 ,                                ## Whether to use nested resampling
                Stage1_FeartureSelection_Method = "cor", cutoff=0,                 ## Stage-1 feature selection method and cutoff
                Stage2_FeartureSelection_Method = "RemoveHighcor",cutoff2 = 0.80,  ## Stage-2 feature selection method and cutoff
                Add_FeartureSelection_Method = "wilcox.test", Unmapped_num = 300,    ## Unmapped feature selection method and cutoff
                classifier2=NULL,                                                  ## Learner for stage 2 prediction(if classifier2==NULL,then it is the same as the learner in stage 1.)
                target='predict',                                                  ## Phenotype prediction
                cores = 5                                                          ## Parallel support
)
...(More Detail)

[1] "-----------------------------------------------------------"
[1] "------------========<<<<  Completed!  >>>>======-----------"
[1] "-----------------------------------------------------------"
[1] "{|>>>=====Learner: liblinear---Performance Metric---==>> AUC:0.953 ACC:0.876 PCCs:0.785 ======<<<|}"
resampling_id learner_name       AUC       ACC      PCCs
1             1    liblinear 0.9642857 0.9285714 0.8309358
2             2    liblinear 0.9387755 0.7857143 0.7114370
3             3    liblinear 0.9132653 0.8214286 0.7304734
4             4    liblinear 0.9940828 0.9615385 0.9101721
5             5    liblinear 0.9526627 0.8846154 0.7415218
Time difference of 10.99379 mins
[1] "######——————  Well Done！！！——————######"

> str(result)
$ Prediction :List of 5
..$ Resample No. 1:'data.frame':	28 obs. of  2 variables:
  .. ..$ sample    : chr [1:28] "201172200006_R04C01" "201172200008_R08C01" "201172200011_R04C01" "201172200012_R07C01" ...
.. ..$ prediction: num [1:28] 0.354 0.781 0.723 0.544 0.401 ...
..$ Resample No. 2:'data.frame':	27 obs. of  2 variables:
  .. ..$ sample    : chr [1:27] "201172200003_R06C01" "201172200004_R06C01" "201172200005_R06C01" "201172200011_R05C01" ...
.. ..$ prediction: num [1:27] 0.97 0.9942 0.2156 0.6554 0.0587 ...
..$ Resample No. 3:'data.frame':	27 obs. of  2 variables:
  .. ..$ sample    : chr [1:27] "201172200001_R07C01" "201172200006_R06C01" "201172200006_R08C01" "201172200019_R07C01" ...
.. ..$ prediction: num [1:27] 0.024 0.9319 0.9827 0.0119 0.9646 ...
..$ Resample No. 4:'data.frame':	27 obs. of  2 variables:
  .. ..$ sample    : chr [1:27] "201172200006_R05C01" "201172200012_R08C01" "201172200013_R04C01" "201172200050_R05C01" ...
.. ..$ prediction: num [1:27] 0.925 0.197 0.914 0.944 0.421 ...
..$ Resample No. 5:'data.frame':	27 obs. of  2 variables:
  .. ..$ sample    : chr [1:27] "201172200001_R06C01" "201172200004_R04C01" "201172200005_R07C01" "201172200016_R04C01" ...
.. ..$ prediction: num [1:27] 0.0504 0.4289 0.1212 0.9876 0.0164 ...
$ Metric     :'data.frame':	5 obs. of  5 variables:
  ..$ resampling_id: int [1:5] 1 2 3 4 5
..$ learner_name : chr [1:5] "liblinear" "liblinear" "liblinear" "liblinear" ...
..$ AUC          : num [1:5] 0.964 0.938 0.913 0.994 0.952
..$ ACC          : num [1:5] 0.928 0.785 0.821 0.961 0.884
..$ PCCs         : num [1:5] 0.830 0.711 0.730 0.910 0.741
$ TotalMetric: Named num [1:3] 0.953 0.876 0.785
..- attr(*, "names")= chr [1:3] "AUC" "ACC" "PCCs"

```

## Parameters optimization
Hyperparameter selection for BioM2 using HyBioM2() to improve prediction performance
```
library(mlr3verse)
library(parallel)
library(caret)
library(BioM2)

#Selection of stage-1 basemodels
classifier1=c('liblinear','svm')


#stage-1 feature_selection
stage1_cutoff=c(0.3,0.5)


#Number of unmapped features
Unmapped_num=c(5,10)


#stage-2 feature_selection
stage2_cutoff=c(0.9,0.8)

#Selection of stage-2 basemodels(The default is the same as the stage-1)
classifier2=NULL

#A data frame contains hyperparameter results
result=HyBioM2(TrainData=TrainData,pathlistDB=pathlistDB,FeatureAnno=FeatureAnno,resampling=NULL,nfolds=2,classifier=classifier1,
           PathwaySizeUp=200,PathwaySizeDown=150,MinfeatureNum_pathways=10,
           Add_FeartureSelection_Method='wilcox.test',Unmapped_num=Unmapped_num,
           Inner_CV=F,inner_folds=10,
           Stage1_FeartureSelection_Method='cor',stage1_cutoff=stage1_cutoff,
           Stage2_FeartureSelection_Method='RemoveHighcor',stage2_cutoff=stage2_cutoff,
           classifier2=NULL,cores=20,verbose=TRUE)


#View the optimal hyperparameter combination
head(result[order(result$AUC,decreasing = T),c(1,3:6,8)])

#   stage1_learner stage1_cutoff stage2_cutoff Unmapped_num       AUC       PCC
#15            svm         0.001            10            5 0.7373737 0.4020263
#2       liblinear         0.010             5           10 0.7171717 0.3806561
#4       liblinear         0.010            10           10 0.7171717 0.3815476
#6       liblinear         0.001             5           10 0.7171717 0.3798767
#8       liblinear         0.001            10           10 0.7171717 0.3809311
#16            svm         0.001            10           10 0.7121212 0.4094221




```
##  Biological interpretability
To explore the potential impact of biological pathways on the disease/phenotype, set the parameter (target='pathways').
Show the association between each biological pathway used for prediction and the phenotype.

```
library(mlr3verse)
library(caret)
library(parallel)
library(BioM2)

result=BioM2 (  TrainData = data , TestData = NULL ,                               
                pathlistDB = pathlistDB ,                                          
                FeatureAnno = FeatureAnno ,                                        
                classifier = 'liblinear' , nfolds = 5 ,                            
                Inner_CV = FALSE , inner_folds=10 ,                                
                Stage1_FeartureSelection_Method = "cor", cutoff=0,                 
                Stage2_FeartureSelection_Method = "RemoveHighcor",cutoff2 = 0.80,    
                classifier2=NULL,                                                  
                target='pathways',            ##==>>  [ target = 'pathways']
                cores = 5                                                          
)

[1] "-----------------------------------------------------------"
[1] "------------========<<<<  Completed!  >>>>======-----------"
[1] "-----------------------------------------------------------"
id       cor      p.value adjust_p.value
254  GO:0001708 0.6421436 4.839108e-17   1.440668e-13
165  GO:0072073 0.6349709 1.102139e-16   3.287145e-13
1141 GO:0072009 0.6263595 3.534259e-16   1.051309e-13
692  GO:0035265 0.6217819 4.435684e-16   1.323840e-12
557  GO:0003341 0.6499000 7.614604e-16   2.264562e-12
term
254                cell fate specification
165          kidney epithelium development
1141        nephron epithelium development
692                           organ growth
557                        cilium movement
Time difference of 8.649476 mins
[1] "######——————  Well Done！！！——————######"

> str(result)
List of 2
$ PathwaysMatrix: num [1:136, 1:2974] 0 1 1 0 0 0 1 0 1 1 ...
..- attr(*, "dimnames")=List of 2
.. ..$ : NULL
.. ..$ : chr [1:2974] "label" "GO:0000002" "GO:0000018" "GO:0000038" ...
$ PathwaysResult:'data.frame':	2973 obs. of  5 variables:
  ..$ id            : chr [1:2973] "GO:0001708" "GO:0072073" "GO:0072009" "GO:0035265" ...
..$ cor           : num [1:2973] 0.64 0.634 0.626 0.625 0.621 ...
..$ p.value       : num [1:2973] 4.83e-17 1.10e-16 3.53e-16 4.43e-16 7.61e-16 ...
..$ adjust_p.value: num [1:2973] 1.44e-13 3.28e-13 1.05e-12 1.32e-12 2.26e-12 ...
..$ term          : chr [1:2973] "cell fate specification" "kidney epithelium development" "nephron epithelium development" "organ growth" ...

```

## Pathway Modules

A pathway matrix can be obtained by using BioM2(, target = 'pathways'). The WGCNA method aggregates pathways with similar expression patterns into a module, and uses biological semantic information to assist in screening modules with high biological interpretability, and compares these biological pathway modules association with phenotype.

### FindParaModule（）：Using Biological Semantic Information to Assist in Selecting Optimal Parameters
```
library(mlr3verse)
library(caret)
library(parallel)
library(BioM2)


result=BioM2 (  TrainData = data , TestData = NULL ,                               
                pathlistDB = pathlistDB ,                                          
                FeatureAnno = FeatureAnno ,                                        
                classifier = 'liblinear' , nfolds = 5 ,                            
                Inner_CV = FALSE , inner_folds=10 ,                                
                Stage1_FeartureSelection_Method = "cor", cutoff=0,                 
                Stage2_FeartureSelection_Method = "RemoveHighcor",cutoff2 = 0.80,  
                classifier2=NULL,                                                  
                target='pathways',                ##==>>  [ target = 'pathways']
                cores = 5                                                          
)
Matrix=result$PathwaysMatrix

library(WGCNA)

Para=FindParaModule(pathways_matrix = Matrix, control_label=0,minModuleSize =seq(5,20,5), mergeCutHeight=seq(0.1,0.4,0.1),power=NULL)

> str(Para)
List of 2
$ TotalResult  :'data.frame':	33 obs. of  6 variables:
  ..$ mergeCutHeight      : num [1:33] 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.1 0.15 0.2 ...
..$ Number_clusters     : int [1:33] 160 156 152 144 122 100 65 115 114 113 ...
..$ Mean_number_pathways: num [1:33] 18.4 18.8 19.3 20.4 21.3 ...
..$ Mean_Fraction       : num [1:33] 65 64.7 64.4 64.3 64.3 ...
..$ Sd_Fraction         : num [1:33] 18.4 18.5 18.2 18.5 18.2 ...
..$ minModuleSize       : num [1:33] 10 10 10 10 10 10 10  ...
$ BestParameter: Named num [1:3] 8 10 0.4
..- attr(*, "names")= chr [1:3] "power" "ModuleSize" "mergeCutHeight"

```
### PathwaysModule（）: Identifying Illness-relevant  Modules with High Biological Interpretability
The optimal parameters can be provided by FindParaModule() or chosen by the user. 
The modules relevant to the illness can then be obtained, with high biological interpretability.
```
library(WGCNA)

Modules=PathwaysModule(pathways_matrix = Matrix , control_label = 0, minModuleSize = 10, mergeCutHeight = 0.4, cutoff = 70,power=NULL)

> str(Modules)
List of 4
$ ModuleResult      :'data.frame':	2939 obs. of  2 variables:
  ..$ ID     : chr [1:2939] "GO:0000002" "GO:0000018" "GO:0000038" "GO:0000041" ...
..$ cluster: num [1:2939] 0 1 62 30 35 27 50 1 50 50 ...
$ RAW_PathwaysModule:'data.frame':	66 obs. of  5 variables:
  ..$ module       : chr [1:66] "ME0" "ME1" "ME10" "ME11" ...
..$ Num_pathways : int [1:66] 200 1351 41 38 38 37 37 35 34 28 ...
..$ Fraction     : num [1:66] 36.5 44.3 73.2 71.1 63.2 ...
..$ adjust_pvalue: num [1:66] 1.57e-06 4.85e-09 9.47e-05 2.03e-03 4.32e-06 ...
..$ cor          : num [1:66] 0.504 0.574 0.392 0.311 0.458 ...
$ DE_PathwaysModule :'data.frame':	9 obs. of  5 variables:
  ..$ module       : chr [1:9] "ME28" "ME52" "ME40" "ME65" ...
..$ Num_pathways : int [1:9] 16 8 11 6 41 17 14 15 7
..$ Fraction     : num [1:9] 81.2 75 90.9 100 73.2 ...
..$ adjust_pvalue: num [1:9] 6.81e-07 1.57e-06 1.92e-05 1.97e-05 9.47e-05 ...
..$ cor          : num [1:9] 0.418 0.454 0.4 0.435 0.392 ...
$ Matrix            : num [1:136, 1:2974] 0 1 1 0 0 0 1 0 1 1 ...
..- attr(*, "dimnames")=List of 2
.. ..$ : NULL
.. ..$ : chr [1:2974] "label" "GO:0000002" "GO:0000018" "GO:0000038" ...

```
### ShowModule（）: Display the Term of the Pathway in Each Pathway-level Module
```
ModulesInner = ShowModule(Modules,c(25,27,34))

> str(ModulesInner)
List of 3
$ ME25:'data.frame':	16 obs. of  4 variables:
  ..$ GO        : chr [1:16] "GO:0006023" "GO:0006024" "GO:0006029" "GO:0015012" ...
..$ Name      : chr [1:16] "aminoglycan biosynthetic process" "glycosaminoglycan biosynthetic process" "proteoglycan metabolic process" "heparan sulfate proteoglycan biosynthetic process" ...
..$ Ancestor  : chr [1:16] "organic substance metabolic process" "organic substance metabolic process" "organic substance metabolic process" "organic substance metabolic process" ...
..$ AncestorGO: chr [1:16] "GO:0071704" "GO:0071704" "GO:0071704" "GO:0071704" ...

$ ME27:'data.frame':	11 obs. of  4 variables:
  ..$ GO        : chr [1:11] "GO:0007173" "GO:0007176" "GO:0038127" "GO:0042058" ...
..$ Name      : chr [1:11] "epidermal growth factor receptor signaling pathway" "regulation of epidermal growth factor-activated receptor activity" "ERBB signaling pathway" "regulation of epidermal growth factor receptor signaling pathway" ...
..$ Ancestor  : chr [1:11] "regulation of cellular process" "regulation of cellular process" "regulation of cellular process" "regulation of cellular process" ...
..$ AncestorGO: chr [1:11] "GO:0050794" "GO:0050794" "GO:0050794" "GO:0050794" ...

$ ME34:'data.frame':	8 obs. of  4 variables:
  ..$ GO        : chr [1:8] "GO:0003351" "GO:0007288" "GO:0030317" "GO:0060294" ...
..$ Name      : chr [1:8] "epithelial cilium movement involved in extracellular fluid movement" "sperm axoneme assembly" "flagellated sperm motility" "cilium movement involved in cell motility" ...
..$ Ancestor  : chr [1:8] "movement of cell or subcellular component" "movement of cell or subcellular component" "movement of cell or subcellular component" "movement of cell or subcellular component" ...
..$ AncestorGO: chr [1:8] "GO:0006928" "GO:0006928" "GO:0006928" "GO:0006928" ...

```



## Visualization
- PlotPathFearture() : Visualisation of significant pathway-level features
- PlotPathInner() : Visualisation Original features that make up the pathway
- PlotPathNet() : Network diagram of pathway-level features
- VisMultiModule() : Visualisation of the results of the analysis of the pathway modules
- PlotCorModule() : Correlalogram for Illness-relevant Modules



### PlotPathFearture ()
**Visualisation of significant pathway-level features**
```
#Load the required R packages
library(ggplot2)
library(viridis)


PlotPathFearture(BioM2_pathways_obj=result , pathlistDB = pathlistDB)
```
![barplot](https://github.com/jkkomm/img/blob/main/barplot.png)




### PlotPathInner ()
**Visualisation Original features that make up the pathway**
```
#Load the required R packages
library(CMplot)

#Select the top 10 most significant pathways
PathNames=result$PathwaysResult$id[1:10]

PlotPathInner(data=TrainData,pathlistDB=pathlistDB,FeatureAnno=FeatureAnno,
              PathNames=PathNames)
```
![CMantan](https://github.com/jkkomm/img/blob/main/CMantan.png)




### PlotPathNet ()
**Network diagram of pathway-level features**
```
#Load the required R packages
library(igraph)
library(ggnetwork)
library(ggplot2)


#Select the top 10 most significant pathways
PathNames=result$PathwaysResult$id[1:10]

PlotPathNet(data=TrainData,BioM2_pathways_obj=result,
            FeatureAnno = FeatureAnno,pathlistDB=pathlistDB,
            PathNames=PathNames)

```
![network](https://github.com/jkkomm/img/blob/main/network.jpg)



### VisMultiModule

**VisMultiModule ( ,BioM2_pathways_obj ) ：**
Each point represents a pathway, and each pathway belongs to a biological category. The higher the point, the more significant the difference between the pathway and the phenotype
```
#Load the required R packages
library(ggpubr)
library(ggthemes)
library(CMplot)
library(ggplot2)
library(RColorBrewer)
library(webshot)
library(wordcloud2)
library(jiebaR)
library("htmlwidgets")


VisMultiModule(BioM2_pathways_obj = result)
```
![PathwaysResult](https://github.com/jkkomm/img/blob/main/CManhan2.png)



**VisMultiModule ( , FindParaModule_obj ) ：**
Visualize the process of selecting optimal parameters based on biological terms.
```
Matrix=result$PathwaysMatrix
Para=FindParaModule(pathways_matrix = Matrix, minModuleSize = c(6,7,8), mergeCutHeight=c(0.2,0.25,0.3,0.35,0.4,0.45,0.5),power=NULL)

VisMultiModule(FindParaModule_obj=Para)

```
![FindParaModule](https://github.com/jkkomm/img/blob/main/Figure%205a(1).jpg)



**VisMultiModule ( , PathwaysModule_obj ) ：**
Each point represents a path, and points of the same color belong to the same illness-relevant module
```
Matrix=result$PathwaysMatrix
Modules=PathwaysModule(pathways_matrix = Matrix , control_label = 0, minModuleSize = 6, mergeCutHeight = 0.3, cutoff = 70,power=NULL)

VisMultiModule(PathwaysModule_obj=Modules)
```
![DE_PathwaysModule](https://github.com/jkkomm/img/blob/main/Figure%205b(1).jpg)

Violin plot showing statistics for the pathway modules
```
# xxx -> module index
VisMultiModule(PathwaysModule_obj=Modules,volin=TRUE,control_label=0,module= xxx )

img=list()

 ##(Setting up a visual module index)
modules=c(14,15,28,4)

for(i in 1:length(xxx)){
  pic=VisMultiModule(PathwaysModule_obj=Modules,volin=TRUE,control_label=0,module=modules[i])
  img[[i]]=pic
}
d=1
img[[1]]+theme(plot.margin = unit(c(d,d,d,d), "cm"))+
  img[[2]]+theme(plot.margin = unit(c(d,d,d,d), "cm"))+
  img[[3]]+theme(plot.margin = unit(c(d,d,d,d), "cm"))+
  img[[4]]+theme(plot.margin = unit(c(d,d,d,d), "cm"))
```
![Violin](https://github.com/jkkomm/img/blob/main/Figure%208a(1).jpg)

**VisMultiModule ( , ShowModule_obj )  ：**
Summarize the biological information of the pathways in the module with a wordcloud.
```
# xxx -> module index
ShowModule(Modules,xxx)

VisMultiModule(ShowModule_obj=ModulesInner)


 ##(Setting up a visual module index)
modules=c(14,15,28,4)

output=paste0('Module',modules,'_WordCloud.png')
for(i in 1:length(modules)){
  ModulesInner = ShowModule(Modules,modules[i])
  my_graph=VisMultiModule(ShowModule_obj=ModulesInner)
  saveWidget(my_graph,'tmp.html',selfcontained = F)
  webshot('tmp.html',output[i])
}
library(magick)
ff=list.files(pattern = 'png')
gg=do.call(cbind,lapply(1:4,function(x) image_read(ff[x])))
library(patchwork)
library(ggthemes)
ff2=gsub('_WordCloud.png','',ff)
ff2=gsub('ME','Module',ff2)
ba='#FFFFFF'
d=2
image_ggplot(gg[[1]])+labs(title = ff2[1])+
  theme(text = element_text(family = "serif", size = 10, color = "black",face='bold'),
        panel.border = element_rect(color = "black", fill = NA, size = 2),
        plot.margin = unit(c(d, d, d, d), "mm"),
        panel.background = element_rect(fill = ba, color = ba),
        plot.background = element_rect(fill = ba, color = ba))+
  image_ggplot(gg[[2]])+labs(title = ff2[2])+
  theme(text = element_text(family = "serif", size = 10, color = "black",face='bold'),
        panel.border = element_rect(color = "black", fill = NA, size = 2),
        plot.margin = unit(c(d, d, d, d), "mm"),
        panel.background = element_rect(fill = ba, color = ba),
        plot.background = element_rect(fill = ba, color = ba))+
  image_ggplot(gg[[3]])+labs(title = ff2[3])+
  theme(text = element_text(family = "serif", size = 10, color = "black",face='bold'),
        panel.border = element_rect(color = "black", fill = NA, size = 2),
        plot.margin = unit(c(d, d, d, d), "mm"),
        panel.background = element_rect(fill = ba, color = ba),
        plot.background = element_rect(fill = ba, color = ba))+
  image_ggplot(gg[[4]])+labs(title = ff2[4])+
  theme(text = element_text(family = "serif", size = 10, color = "black",face='bold'),
        panel.border = element_rect(color = "black", fill = NA, size = 2),
        plot.margin = unit(c(d, d, d, d), "mm"),
        panel.background = element_rect(fill = ba, color = ba),
        plot.background = element_rect(fill = ba, color = ba))





```
![SM25](https://github.com/jkkomm/img/blob/main/WordCloud.png)

### PlotCorModule()
**Correlalogram for illness-relevant modules**
```
PlotCorModule(PathwaysModule_obj=Modules)

```
![Cor](https://github.com/jkkomm/img/blob/main/Figure%206a(1).jpg)
# Citation
- NIPS ML4H submission: Chen, J. and Schwarz, E., 2017. BioMM: Biologically-informed Multi-stage Machine learning for identification of epigenetic fingerprints. arXiv preprint arXiv:1712.00336.
- Chen, Junfang, et al. "Association of a reproducible epigenetic risk profile for schizophrenia with brain methylation and function." JAMA psychiatry 77.6 (2020): 628-636.





